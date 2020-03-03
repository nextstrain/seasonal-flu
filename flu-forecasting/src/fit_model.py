"""Fit a model for the given data using the requested predictors and evaluate the model by time series cross-validation.
"""
import argparse
import csv
import json
import numpy as np
import pandas as pd
from scipy.optimize import minimize
import sys
import time

from forecast.fitness_model import get_train_validate_timepoints
from forecast.metrics import add_pseudocounts_to_frequencies, negative_information_gain
from forecast.metrics import mean_absolute_error, sum_of_squared_errors, root_mean_square_error
from weighted_distances import get_distances_by_sample_names, get_distance_matrix_by_sample_names

MAX_PROJECTED_FREQUENCY = 1e3
FREQUENCY_TOLERANCE = 1e-3

np.random.seed(314159)


def sum_of_differences(observed, estimated, y_diff, **kwargs):
    """
    Calculates the sum of squared errors for observed and estimated values.

    Parameters
    ----------
    observed : numpy.ndarray
        observed values

    estimated : numpy.ndarray
        estimated values

    y_diff : numpy.ndarray
        differences between observed and estimated values

    Returns
    -------
    float :
        sum of differences between estimated and observed future values
    """
    return np.sum(y_diff)


class ExponentialGrowthModel(object):
    def __init__(self, predictors, delta_time, l1_lambda, cost_function):
        """Construct an empty exponential growth model instance.

        Parameters
        ----------
        predictors : list
            a list of predictors to estimate coefficients for

        delta_time : float
            number of years into the future to project frequencies

        l1_lambda : float
            hyperparameter to scale L1 regularization penalty for non-zero coefficients

        cost_function : callable
            function returning the error to be minimized between observed and estimated values

        Returns
        -------
        ExponentialGrowthModel
        """
        self.predictors = predictors
        self.delta_time = delta_time
        self.l1_lambda = l1_lambda
        self.cost_function = cost_function

    def calculate_mean_stds(self, X, predictors):
        """Calculate mean standard deviations of predictors by timepoints prior to
        fitting.

        Parameters
        ----------
        X : pandas.DataFrame
            standardized tip attributes by timepoint

        predictors : list
            names of predictors for which mean standard deviations should be calculated

        Returns
        -------
        ndarray :
            mean standard deviation per predictor across all timepoints

        """
        # Note that the pandas standard deviation method ignores missing data
        # whereas numpy requires the use of specific NaN-aware functions (nanstd).
        return X.loc[:, ["timepoint"] + predictors].groupby("timepoint").std().mean().values

    def standardize_predictors(self, predictors, mean_stds, initial_frequencies):
        """Standardize the values for the given predictors by centering on the mean of
        each predictor and scaling by the mean standard deviation provided.

        Parameters
        ----------
        predictors : ndarray
            matrix of values per sample (rows) and predictor (columns)

        mean_stds : ndarray
            mean standard deviations of predictors across all training
            timepoints

        initial_frequencies : ndarray
            initial frequencies of samples corresponding to each row of the
            given predictors

        Returns
        -------
        ndarray :
            standardized predictor values

        """
        means = np.average(predictors, weights=initial_frequencies, axis=0)
        variances = np.average((predictors - means) ** 2, weights=initial_frequencies, axis=0)
        stds = np.sqrt(variances)

        nonzero_stds = np.where(stds)[0]

        if len(nonzero_stds) == 0:
            return predictors

        standardized_predictors = predictors
        standardized_predictors[:, nonzero_stds] = (predictors[:, nonzero_stds] - means[nonzero_stds]) / stds[nonzero_stds]

        return standardized_predictors

    def get_fitnesses(self, coefficients, predictors):
        """Apply the coefficients to the predictors and sum them to get strain
        fitnesses.

        Parameters
        ----------
        coefficients : ndarray or list
            coefficients for given predictors

        predictors : ndarray
            predictor values per sample (n x p matrix for p predictors and n samples)

        Returns
        -------
        ndarray :
            fitnesses per sample
        """
        return np.sum(predictors * coefficients, axis=-1)

    def project_frequencies(self, initial_frequencies, fitnesses, delta_time):
        """Project the given initial frequencies into the future by the given delta time
        based on the given fitnesses.

        Returns the projected frequencies normalized to sum to 1.

        Parameters
        ----------
        initial_frequencies : ndarray
            floating point frequencies for all samples in a timepoint

        fitnesses : ndarray
            floating point fitnesses for all samples in same order as given frequencies

        delta_time : float
            number of years to project into the future

        Returns
        -------
        ndarray :
            projected and normalized frequencies
        """
        # Exponentiate the fitnesses and multiply them by strain frequencies.
        projected_frequencies = initial_frequencies * np.exp(fitnesses * self.delta_time)

        # Replace infinite values a very large number that can still be summed
        # across all timepoints. This addresses the case of buffer overflows in
        # exponentiation which can produce both of these problematic values.
        projected_frequencies[np.isinf(projected_frequencies)] = MAX_PROJECTED_FREQUENCY

        # Sum the projected frequencies.
        total_projected_frequencies = projected_frequencies.sum()

        # Normalize the projected frequencies.
        projected_frequencies = projected_frequencies / total_projected_frequencies

        # Confirm that projected frequencies sum to 1.
        assert np.isclose(projected_frequencies.sum(), np.ones(1), atol=FREQUENCY_TOLERANCE)

        # Confirm that all projected frequencies are proper numbers.
        assert np.isnan(projected_frequencies).sum() == 0

        return projected_frequencies

    def _fit(self, coefficients, X, y, use_l1_penalty=True):
        """Calculate the error between observed and estimated values for the given
        parameters and data.

        Parameters
        ----------
        coefficients : ndarray
            coefficients for each of the model's predictors

        X : pandas.DataFrame
            standardized tip attributes by timepoint

        y : pandas.DataFrame
            final clade frequencies at delta time in the future from each
            timepoint in the given tip attributes table

        Returns
        -------
        float :
            error between estimated values using the given coefficients and
            input data and the observed values
        """
        # Estimate final frequencies.
        y_hat = self.predict(X, coefficients)

        # Merge estimated and observed frequencies. The left join enables
        # tracking of clades that die in the future and are therefore not
        # observed in the future frequencies data frame.
        frequencies = y_hat.merge(
            y,
            how="left",
            on=["timepoint", "clade_membership"],
            suffixes=["_estimated", "_observed"]
        )
        frequencies["frequency_observed"] = frequencies["frequency_observed"].fillna(0.0)

        # Calculate initial frequencies for use by cost function.
        initial_frequencies = X.groupby([
            "timepoint",
            "clade_membership"
        ])["frequency"].sum().reset_index()

        # Annotate future frequencies with initial frequencies.
        frequencies = frequencies.merge(
            initial_frequencies,
            how="inner",
            on=["timepoint", "clade_membership"]
        )

        # Calculate the error between the observed and estimated frequencies.
        error = self.cost_function(
            frequencies["frequency_observed"],
            frequencies["frequency_estimated"],
            initial=frequencies["frequency"]
        )

        if use_l1_penalty:
            l1_penalty = self.l1_lambda * np.abs(coefficients).sum()
        else:
            l1_penalty = 0.0

        return error + l1_penalty

    def fit(self, X, y):
        """Fit a model to the given input data, producing beta coefficients for each of
        the model's predictors.

        Coefficients are stored in the `coef_` attribute, after the pattern of
        scikit-learn models.

        Parameters
        ----------
        X : pandas.DataFrame
            standardized tip attributes by timepoint

        y : pandas.DataFrame
            final clade frequencies at delta time in the future from each
            timepoint in the given tip attributes table

        Returns
        -------
        float :
            model training error

        """
        # Calculate mean standard deviations of predictors by timepoints prior
        # to fitting.
        self.mean_stds_ = self.calculate_mean_stds(X, self.predictors)

        # Find coefficients that minimize the model's cost function.
        if hasattr(self, "coef_"):
            # Use the previous coefficients +/- a small random offset (+/- 0.05)
            # to prevent getting stuck in local minima.
            initial_coefficients = self.coef_ + (0.1 * np.random.random(len(self.predictors)) - 0.05)
        else:
            # If no previous coefficients exist, sample random values between -0.5 and 0.5.
            initial_coefficients = np.random.random(len(self.predictors)) - 0.5

        results = minimize(
            self._fit,
            initial_coefficients,
            args=(X, y),
            method="Nelder-Mead",
            options={"disp": False}
        )
        self.coef_ = results.x

        training_error = self.score(X, y)

        return training_error

    def predict(self, X, coefficients=None, mean_stds=None):
        """Calculate the estimate final frequencies of all clades in the given tip
        attributes data frame using previously calculated beta coefficients.

        Parameters
        ----------
        X : pandas.DataFrame
            standardized tip attributes by timepoint

        coefficients : ndarray
            optional coefficients to use for each of the model's predictors
            instead of the model's currently defined coefficients

        mean_stds : ndarray
            optional mean standard deviations of predictors across all training
            timepoints

        Returns
        -------
        pandas.DataFrame
            estimated final clade frequencies at delta time in the future for
            each clade from each timepoint in the given tip attributes table

        """
        # Use model coefficients, if none are provided.
        if coefficients is None:
            coefficients = self.coef_

        if mean_stds is None:
            mean_stds = self.mean_stds_

        estimated_frequencies = []
        for timepoint, timepoint_df in X.groupby("timepoint"):
            # Select frequencies from timepoint.
            initial_frequencies = timepoint_df["frequency"].values

            # Select predictors from the timepoint.
            predictors = timepoint_df.loc[:, self.predictors].values

            # Standardize predictors by timepoint centering by means at
            # timepoint and mean standard deviation provided.
            standardized_predictors = self.standardize_predictors(predictors, mean_stds, initial_frequencies)

            # Calculate fitnesses.
            fitnesses = self.get_fitnesses(coefficients, standardized_predictors)

            # Project frequencies.
            projected_frequencies = self.project_frequencies(
                initial_frequencies,
                fitnesses,
                self.delta_time
            )

            # Sum the estimated frequencies by clade.
            projected_timepoint_df = timepoint_df[["timepoint", "clade_membership"]].copy()
            projected_timepoint_df["frequency"] = projected_frequencies
            projected_clade_frequencies = projected_timepoint_df.groupby([
                "timepoint",
                "clade_membership"
            ])["frequency"].sum().reset_index()

            estimated_frequencies.append(projected_clade_frequencies)

        # Collect all estimated frequencies by timepoint.
        estimated_frequencies = pd.concat(estimated_frequencies)
        return estimated_frequencies

    def score(self, X, y):
        """Calculate model error between the estimated final clade frequencies for the
        given tip attributes, `X`, and the observed final clade frequencies in
        `y`.

        Parameters
        ----------
        X : pandas.DataFrame
            standardized tip attributes by timepoint

        y : pandas.DataFrame
            final clade frequencies at delta time in the future from each
            timepoint in the given tip attributes table

        Returns
        -------
        float :
            model error
        """
        return self._fit(self.coef_, X, y, use_l1_penalty=False)


class DistanceExponentialGrowthModel(ExponentialGrowthModel):
    def __init__(self, predictors, delta_time, l1_lambda, cost_function, distances):
        super().__init__(predictors, delta_time, l1_lambda, cost_function)
        self.distances = distances

    def _fit(self, coefficients, X, y, use_l1_penalty=True):
        """Calculate the error between observed and estimated values for the given
        parameters and data.

        Parameters
        ----------
        coefficients : ndarray
            coefficients for each of the model's predictors

        X : pandas.DataFrame
            standardized tip attributes by timepoint

        y : pandas.DataFrame
            final weighted distances at delta time in the future from each
            timepoint in the given tip attributes table

        Returns
        -------
        float :
            error between estimated values using the given coefficients and
            input data and the observed values
        """
        import cv2

        # Estimate target values.
        y_hat = self.predict(X, coefficients)

        # Calculate EMD for each timepoint in the estimated values and sum that
        # distance across all timepoints.
        error = 0.0
        count = 0
        for timepoint, timepoint_df in y_hat.groupby("timepoint"):
            samples_a = timepoint_df["strain"]
            sample_a_initial_frequencies = timepoint_df["frequency"].values.astype(np.float32)
            sample_a_frequencies = timepoint_df["projected_frequency"].values.astype(np.float32)

            future_timepoint_df = y[y["timepoint"] == timepoint]
            assert future_timepoint_df.shape[0] > 0

            samples_b = future_timepoint_df["strain"]
            sample_b_frequencies = future_timepoint_df["frequency"].values.astype(np.float32)

            distance_matrix = get_distance_matrix_by_sample_names(
                samples_a,
                samples_b,
                self.distances
            ).astype(np.float32)

            # Estimate the distance between the model's estimated future and the
            # observed future populations.
            model_emd, _, self.model_flow = cv2.EMD(
                sample_a_frequencies,
                sample_b_frequencies,
                cv2.DIST_USER,
                cost=distance_matrix
            )

            error += model_emd
            count += 1

        error = error / float(count)

        if use_l1_penalty:
            l1_penalty = self.l1_lambda * np.abs(coefficients).sum()
        else:
            l1_penalty = 0.0

        return error + l1_penalty

    def _fit_distance(self, coefficients, X, y, use_l1_penalty=True):
        """Calculate the error between observed and estimated values for the given
        parameters and data.

        Parameters
        ----------
        coefficients : ndarray
            coefficients for each of the model's predictors

        X : pandas.DataFrame
            standardized tip attributes by timepoint

        y : pandas.DataFrame
            final weighted distances at delta time in the future from each
            timepoint in the given tip attributes table

        Returns
        -------
        float :
            error between estimated values using the given coefficients and
            input data and the observed values
        """
        # Estimate target values.
        y_hat = self.predict(X, coefficients)

        # Calculate weighted distance to the future for each timepoint in the
        # estimated values and sum that distance across all timepoints.
        error = 0.0
        null_error = 0.0
        count = 0
        for timepoint, timepoint_df in y_hat.groupby("timepoint"):
            samples_a = timepoint_df["strain"]
            sample_a_initial_frequencies = timepoint_df["frequency"].values
            sample_a_frequencies = timepoint_df["projected_frequency"].values
            sample_a_weighted_distance_to_future = timepoint_df["weighted_distance_to_future"].values

            future_timepoint_df = y[y["timepoint"] == timepoint]
            assert future_timepoint_df.shape[0] > 0

            samples_b = future_timepoint_df["strain"]
            sample_b_frequencies = future_timepoint_df["frequency"].values
            sample_b_weighted_distance_to_present = future_timepoint_df["weighted_distance_to_present"].values

            d_t_u = (sample_a_initial_frequencies * sample_a_weighted_distance_to_future).sum()
            d_u_hat_u = (sample_a_frequencies * sample_a_weighted_distance_to_future).sum()
            d_u_u = (sample_b_frequencies * sample_b_weighted_distance_to_present).sum()

            null_error += d_t_u

            error += (d_u_hat_u - d_u_u) / d_t_u
            count += 1

        null_error = null_error / float(count)
        error = error / float(count)

        if use_l1_penalty:
            l1_penalty = self.l1_lambda * np.abs(coefficients).sum()
        else:
            l1_penalty = 0.0

        return error + l1_penalty

    def predict(self, X, coefficients=None, mean_stds=None):
        """Calculate the estimated final weighted distance between tips at each
        timepoint and at that timepoint plus delta months in the future.

        Parameters
        ----------
        X : pandas.DataFrame
            standardized tip attributes by timepoint

        coefficients : ndarray
            optional coefficients to use for each of the model's predictors
            instead of the model's currently defined coefficients

        mean_stds : ndarray
            optional mean standard deviations of predictors across all training
            timepoints

        Returns
        -------
        pandas.DataFrame
            estimated weighted distances at delta time in the future for
            each tip from each timepoint in the given tip attributes table

        """
        # Use model coefficients, if none are provided.
        if coefficients is None:
            coefficients = self.coef_
            model_is_fit = True
        else:
            model_is_fit = False

        if mean_stds is None:
            mean_stds = self.mean_stds_

        estimated_targets = []
        for timepoint, timepoint_df in X.groupby("timepoint"):
            # Select frequencies from timepoint.
            initial_frequencies = timepoint_df["frequency"].values

            # Select predictors from the timepoint.
            predictors = timepoint_df.loc[:, self.predictors].values

            # Standardize predictors by timepoint centering by means at
            # timepoint and mean standard deviation provided.
            mean_stds = timepoint_df.loc[:, self.predictors].std().values
            standardized_predictors = self.standardize_predictors(predictors, mean_stds, initial_frequencies)

            # Calculate fitnesses.
            fitnesses = self.get_fitnesses(coefficients, standardized_predictors)

            # Project frequencies.
            projected_frequencies = self.project_frequencies(
                initial_frequencies,
                fitnesses,
                self.delta_time
            )

            # Calculate observed distance between current tips and the future
            # using projected frequencies and weighted distances to the future.
            columns_to_extract = ["timepoint", "strain", "frequency"]
            optional_columns = ["weighted_distance_to_present", "weighted_distance_to_future"]
            for column in optional_columns:
                if column in timepoint_df.columns:
                    columns_to_extract.append(column)

            projected_timepoint_df = timepoint_df[columns_to_extract].copy()
            projected_timepoint_df["fitness"] = fitnesses
            projected_timepoint_df["projected_frequency"] = projected_frequencies

            if model_is_fit:
                # Calculate estimate distance between current tips and future tips
                # based on projections of current tips.
                estimated_weighted_distance_to_future = []
                for current_tip, current_tip_frequency in projected_timepoint_df.loc[:, ["strain", "frequency"]].values:
                    weighted_distance_to_future = 0.0
                    for other_tip, other_tip_projected_frequency in projected_timepoint_df.loc[:, ["strain", "projected_frequency"]].values:
                        weighted_distance_to_future += other_tip_projected_frequency * self.distances[current_tip][other_tip]

                    estimated_weighted_distance_to_future.append(weighted_distance_to_future)

                projected_timepoint_df["y"] = np.array(estimated_weighted_distance_to_future)
            else:
                projected_timepoint_df["y"] = np.nan

            estimated_targets.append(projected_timepoint_df)

        # Collect all estimated targets by timepoint.
        estimated_targets = pd.concat(estimated_targets, ignore_index=True)
        return estimated_targets


def cross_validate(model_class, model_kwargs, data, targets, train_validate_timepoints, coefficients=None, group_by="clade_membership",
                   include_attributes=False):
    """Calculate cross-validation scores for the given data and targets across the
    given train/validate timepoints.

    Parameters
    ----------
    model : ExponentialGrowthModel
        an instance of a model with defined hyperparameters including a list of
        predictors to use for fitting

    data : pandas.DataFrame
        standardized input attributes to use for model fitting

    targets : pandas.DataFrame
        observed outputs to fit the model to

    train_validate_timepoints : list
        a list of dictionaries of lists indexed by "train" and "validate" keys
        and containing timepoints to use for model training and validation,
        respectively

    coefficients : ndarray
        an optional array of fixed coefficients for the given model's predictors
        to use when calculating cross-validation error for specific models
        (e.g., naive forecasts)

    group_by : string
        column of the tip attributes by which they should be grouped to
        calculate the total number of samples in the model (e.g., group by clade
        or strain)

    include_attributes : boolean
        specifies whether tip attribute data used to train/validate models
        should be included in the output per training window

    Returns
    -------
    list
        a list of dictionaries containing cross-validation results with scores,
        training and validation results, and beta coefficients per timepoint

    """
    results = []
    differences_of_model_and_naive_errors = []
    previous_coefficients = None

    for timepoints in train_validate_timepoints:
        model = model_class(**model_kwargs)

        if previous_coefficients is not None:
            model.coef_ = previous_coefficients

        # Get training and validation timepoints.
        training_timepoints = pd.to_datetime(timepoints["train"])
        validation_timepoint = pd.to_datetime(timepoints["validate"])

        # Get training data by timepoints.
        training_X = data[data["timepoint"].isin(training_timepoints)].copy()
        training_y = targets[targets["timepoint"].isin(training_timepoints)].copy()

        # Fit a model to the training data.
        if coefficients is None:
            start_time = time.time()
            training_error = model.fit(training_X, training_y)
            end_time = time.time()
            previous_coefficients = model.coef_
            null_training_error = model._fit(np.zeros_like(model.coef_), training_X, training_y)
        else:
            start_time = end_time = time.time()
            model.coef_ = coefficients
            model.mean_stds_ = model.calculate_mean_stds(training_X, model.predictors)
            training_error = model.score(training_X, training_y)
            null_training_error = training_error

        # Get validation data by timepoints.
        validation_X = data[data["timepoint"] == validation_timepoint].copy()
        validation_y = targets[targets["timepoint"] == validation_timepoint].copy()

        # Calculate the model score for the validation data.
        validation_error = model.score(validation_X, validation_y)
        null_validation_error = model._fit(np.zeros_like(model.coef_), validation_X, validation_y)
        differences_of_model_and_naive_errors.append(validation_error - null_validation_error)
        print(
            "%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%s\t%.2f\t%.2f" % (
                training_timepoints[-1].strftime("%Y-%m"),
                validation_timepoint.strftime("%Y-%m"),
                training_error,
                null_training_error,
                validation_error,
                null_validation_error,
                model.coef_,
                (np.array(differences_of_model_and_naive_errors) < 0).sum() / float(len(differences_of_model_and_naive_errors)),
                end_time - start_time
            ),
            flush=True
        )

        # Get the estimated frequencies for training and validation sets to export.
        training_y_hat = model.predict(training_X)
        validation_y_hat = model.predict(validation_X)

        # Convert timestamps to a serializable format.
        for df in [training_X, training_y, training_y_hat, validation_X, validation_y, validation_y_hat]:
            for column in ["timepoint", "future_timepoint"]:
                if column in df.columns:
                    df[column] = df[column].dt.strftime("%Y-%m-%d")

        # Store training results, beta coefficients, and validation results.
        result = {
            "predictors": model.predictors,
            "training_data": {
                "y": training_y.to_dict(orient="records"),
                "y_hat": training_y_hat.to_dict(orient="records")
            },
            "training_n": training_X[group_by].unique().shape[0],
            "training_error": training_error,
            "coefficients": model.coef_.tolist(),
            "mean_stds": model.mean_stds_.tolist(),
            "validation_data": {
                "y": validation_y.to_dict(orient="records"),
                "y_hat": validation_y_hat.to_dict(orient="records")
            },
            "validation_n": validation_X[group_by].unique().shape[0],
            "validation_error": validation_error,
            "null_validation_error": null_validation_error,
            "last_training_timepoint": training_timepoints[-1].strftime("%Y-%m-%d"),
            "validation_timepoint": validation_timepoint.strftime("%Y-%m-%d")
        }

        # Include tip attributes, if requested.
        if include_attributes:
            result["training_data"]["X"] = training_X.to_dict(orient="records")
            result["validation_data"]["X"] = validation_X.to_dict(orient="records")

        results.append(result)

    # Return results for all validation timepoints.
    print("Mean difference between model and naive: %.4f" % (sum(differences_of_model_and_naive_errors) / len(differences_of_model_and_naive_errors)), flush=True)
    print("Proportion of timepoints when model < naive: %.2f" % ((np.array(differences_of_model_and_naive_errors) < 0).sum() / float(len(differences_of_model_and_naive_errors))))
    return results


def test(model_class, model_kwargs, data, targets, timepoints, coefficients=None, group_by="clade_membership",
                   include_attributes=False):
    """Calculate test scores for the given data and targets across the given
    timepoints.

    Parameters
    ----------
    model : ExponentialGrowthModel
        an instance of a model with defined hyperparameters including a list of
        predictors to use for fitting

    data : pandas.DataFrame
        standardized input attributes to use for model fitting

    targets : pandas.DataFrame
        observed outputs to test the model with

    timepoints : list
        a list of timepoint strings in YYYY-MM-DD format

    coefficients : ndarray
        an array of fixed coefficients for the given model's predictors

    group_by : string
        column of the tip attributes by which they should be grouped to
        calculate the total number of samples in the model (e.g., group by clade
        or strain)

    include_attributes : boolean
        specifies whether tip attribute data used to test models should be
        included in the output per timepoint

    Returns
    -------
    list
        a list of dictionaries containing test results with scores per timepoint

    """
    results = []
    differences_of_model_and_naive_errors = []

    for timepoint in timepoints:
        model = model_class(**model_kwargs)
        model.coef_ = coefficients
        model.mean_stds_ = np.zeros_like(coefficients)

        # Get training and validation timepoints.
        test_timepoint = pd.to_datetime(timepoint)

        # Get test data by timepoints.
        test_X = data[data["timepoint"] == test_timepoint].copy()
        test_y = targets[targets["timepoint"] == test_timepoint].copy()

        # Calculate the model score for the validation data.
        test_error = model.score(test_X, test_y)
        null_test_error = model._fit(np.zeros_like(model.coef_), test_X, test_y)
        differences_of_model_and_naive_errors.append(test_error - null_test_error)
        print(
            "%s\t%.2f\t%.2f\t%s\t%.2f" % (
                test_timepoint.strftime("%Y-%m"),
                test_error,
                null_test_error,
                model.coef_,
                (np.array(differences_of_model_and_naive_errors) < 0).sum() / float(len(differences_of_model_and_naive_errors))
            ),
            flush=True
        )

        # Get the estimated frequencies for test sets to export.
        test_y_hat = model.predict(test_X)

        # Convert timestamps to a serializable format.
        for df in [test_X, test_y, test_y_hat]:
            for column in ["timepoint", "future_timepoint"]:
                if column in df.columns:
                    df[column] = df[column].dt.strftime("%Y-%m-%d")

        # Store test results and beta coefficients.
        result = {
            "predictors": model.predictors,
            "coefficients": model.coef_.tolist(),
            "mean_stds": model.mean_stds_.tolist(),
            "validation_data": {
                "y": test_y.to_dict(orient="records"),
                "y_hat": test_y_hat.to_dict(orient="records")
            },
            "validation_n": test_X[group_by].unique().shape[0],
            "validation_error": test_error,
            "null_validation_error": null_test_error,
            "validation_timepoint": test_timepoint.strftime("%Y-%m-%d")
        }

        # Include tip attributes, if requested.
        if include_attributes:
            result["validation_data"]["X"] = test_X.to_dict(orient="records")

        results.append(result)

    # Return results for all validation timepoints.
    print("Mean difference between model and naive: %.4f" % (sum(differences_of_model_and_naive_errors) / len(differences_of_model_and_naive_errors)), flush=True)
    print("Proportion of timepoints when model < naive: %.2f" % ((np.array(differences_of_model_and_naive_errors) < 0).sum() / float(len(differences_of_model_and_naive_errors))))
    return results


def summarize_scores(scores, include_scores=False):
    """Summarize model errors across timepoints.

    Parameters
    ----------
    scores : list
        a list of cross-validation results including training errors,
        cross-validation errors, and beta coefficients OR a list of test errors

    include_scores : boolean
        specifies whether cross-validation scores should be included in the
        output per timepoint

    Returns
    -------
    dict :
        a dictionary of all cross-validation results plus summary statistics for
        training, cross-validation, and beta coefficients OR test results

    """
    summary = {
        "predictors": scores[0]["predictors"]
    }

    if include_scores:
        summary["scores"] = scores

    validation_errors = [score["validation_error"] for score in scores]
    summary["cv_error_mean"] = np.mean(validation_errors)
    summary["cv_error_std"] = np.std(validation_errors)

    coefficients = np.array([
        np.array(score["coefficients"])
        for score in scores
    ])
    summary["coefficients_mean"] = coefficients.mean(axis=0).tolist()
    summary["coefficients_std"] = coefficients.std(axis=0).tolist()

    mean_stds = np.array([
        np.array(score["mean_stds"])
        for score in scores
    ])
    summary["mean_stds_mean"] = mean_stds.mean(axis=0).tolist()
    summary["mean_stds_std"] = mean_stds.std(axis=0).tolist()

    return summary


def get_errors_by_timepoint(scores):
    """Convert cross-validation errors into a data frame by timepoint and predictors.

    Parameters
    ----------
    scores : list
        a list of cross-validation results including training errors,
        cross-validation errors, and beta coefficients

    Returns
    -------
    pandas.DataFrame
    """
    predictors = "-".join(scores[0]["predictors"])
    errors_by_time = []
    for score in scores:
        errors_by_time.append({
            "predictors": predictors,
            "validation_timepoint": pd.to_datetime(score["validation_timepoint"]),
            "validation_error": score["validation_error"],
            "null_validation_error": score["null_validation_error"],
            "validation_n": score["validation_n"]
        })

    return pd.DataFrame(errors_by_time)


def get_coefficients_by_timepoint(scores):
    """Convert model coefficients into a data frame by timepoint and predictors.

    Parameters
    ----------
    scores : list
        a list of cross-validation results including training errors,
        cross-validation errors, and beta coefficients

    Returns
    -------
    pandas.DataFrame
    """
    predictors = "-".join(scores[0]["predictors"])
    coefficients_by_time = []
    for score in scores:
        for predictor, coefficient in zip(score["predictors"], score["coefficients"]):
            coefficients_by_time.append({
                "predictors": predictors,
                "predictor": predictor,
                "coefficient": coefficient,
                "validation_timepoint": pd.to_datetime(score["validation_timepoint"])
            })

    return pd.DataFrame(coefficients_by_time)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tip-attributes", required=True, help="tab-delimited file describing tip attributes at all timepoints with standardized predictors")
    parser.add_argument("--output", required=True, help="JSON representing the model fit with training and cross-validation results, beta coefficients for predictors, and summary statistics")
    parser.add_argument("--predictors", nargs="+", help="tip attribute columns to use as predictors of final clade frequencies; optional if a fixed model is provided")
    parser.add_argument("--delta-months", type=int, help="number of months to project clade frequencies into the future")
    parser.add_argument("--target", required=True, choices=["clades", "distances"], help="target for models to fit")
    parser.add_argument("--final-clade-frequencies", help="tab-delimited file of clades per timepoint and their corresponding tips and tip frequencies at the given delta time in the future")
    parser.add_argument("--distances", help="tab-delimited file of distances between pairs of samples")
    parser.add_argument("--training-window", type=int, default=4, help="number of years required for model training")
    parser.add_argument("--l1-lambda", type=float, default=0.0, help="L1 regularization lambda")
    parser.add_argument("--cost-function", default="sse", choices=["sse", "rmse", "mae", "information_gain", "diffsum"], help="name of the function that returns the error between observed and estimated values")
    parser.add_argument("--pseudocount", type=float, help="pseudocount numerator to adjust all frequencies by, enabling some information theoretic metrics like information gain")
    parser.add_argument("--include-attributes", action="store_true", help="include attribute data used to train/validate models in the cross-validation output")
    parser.add_argument("--include-scores", action="store_true", help="include score data resulting from cross-validation output")
    parser.add_argument("--errors-by-timepoint", help="optional data frame of cross-validation errors by validation timepoint")
    parser.add_argument("--coefficients-by-timepoint", help="optional data frame of coefficients by validation timepoint")
    parser.add_argument("--fixed-model", help="optional model JSON to use as a fixed model for calculation of test error in forecasts")

    args = parser.parse_args()

    cost_functions_by_name = {
        "sse": sum_of_squared_errors,
        "rmse": root_mean_square_error,
        "mae": mean_absolute_error,
        "information_gain": negative_information_gain,
        "diffsum": sum_of_differences
    }

    # Load standardized tip attributes subsetting to tip name, clade, frequency,
    # and requested predictors.
    tips = pd.read_csv(
        args.tip_attributes,
        sep="\t",
        parse_dates=["timepoint"]
    )

    if args.target == "clades":
        # Load final clade tip frequencies.
        final_clade_tip_frequencies = pd.read_csv(
            args.final_clade_frequencies,
            sep="\t",
            parse_dates=["initial_timepoint", "final_timepoint"]
        )

        # If a pseudocount numerator has been provided, update the given tip
        # frequencies both from current and future timepoints.
        if args.pseudocount is not None and args.pseudocount > 0.0:
            tips = add_pseudocounts_to_frequencies(tips, args.pseudocount)
            print("Sum of tip frequencies by timepoint: ",
                  tips.groupby("timepoint")["frequency"].sum())
            final_clade_tip_frequencies = add_pseudocounts_to_frequencies(
                final_clade_tip_frequencies,
                args.pseudocount,
                timepoint_column="initial_timepoint"
            )
            print("Sum of tip frequencies by timepoint: ",
                  final_clade_tip_frequencies.groupby("initial_timepoint")["frequency"].sum())

        # Aggregate final clade frequencies.
        final_clade_frequencies = final_clade_tip_frequencies.groupby([
            "initial_timepoint",
            "clade_membership"
        ])["frequency"].sum().reset_index()

        # Rename initial timepoint column for comparison with tip attribute data.
        targets = final_clade_frequencies.rename(
            columns={"initial_timepoint": "timepoint"}
        )
        model_class = ExponentialGrowthModel
        model_kwargs = {}
        group_by_attribute = "clade_membership"
    elif args.target == "distances":
        # Scale each tip's weighted distance to future populations by one minus
        # the tip's current frequency. This ensures that lower frequency tips do
        # not considered closer to the future.
        tips["y"] = tips["weighted_distance_to_future"]

        # Get strain frequency per timepoint and subtract delta time from
        # timepoint to align strain frequencies with the previous timepoint and
        # make them appropriate as targets for the model.
        targets = tips.loc[:, ["strain", "timepoint", "frequency", "weighted_distance_to_present", "weighted_distance_to_future", "y"]].copy()
        targets["future_timepoint"] = targets["timepoint"]

        model_class = DistanceExponentialGrowthModel

        with open(args.distances, "r") as fh:
            print("Read distances", flush=True)
            reader = csv.DictReader(fh, delimiter="\t")
            print("Get distances by sample names", flush=True)
            distances_by_sample_names = get_distances_by_sample_names(reader)
            print("Data loaded", flush=True)

        model_kwargs = {"distances": distances_by_sample_names}
        group_by_attribute = "strain"

    # Identify all available timepoints from tip attributes.
    timepoints = tips["timepoint"].dt.strftime("%Y-%m-%d").unique()

    # If a fixed model is provided, calculate test errors. Otherwise, calculate
    # cross-validation errors.
    if args.fixed_model is not None:
        # Load model details and extract mean coefficients.
        with open(args.fixed_model, "r") as fh:
            model_json = json.load(fh)

        coefficients = np.array(model_json["coefficients_mean"])
        delta_months = model_json["delta_months"]
        delta_time = delta_months / 12.0
        l1_lambda = model_json["l1_lambda"]
        training_window = model_json["training_window"]
        cost_function_name = model_json["cost_function"]
        cost_function = cost_functions_by_name[cost_function_name]
        model_kwargs.update({
            "predictors": model_json["predictors"],
            "delta_time": delta_time,
            "l1_lambda": l1_lambda,
            "cost_function": cost_function
        })

        # Find the latest timepoint we can project from based on the given delta
        # months.
        latest_timepoint = pd.to_datetime(timepoints[-1]) - pd.DateOffset(months=delta_months)
        test_timepoints = [
            timepoint
            for timepoint in timepoints
            if pd.to_datetime(timepoint) <= latest_timepoint
        ]

        # Calculate test errors/scores for the given coefficients and data at
        # the identified test timepoints.
        targets["timepoint"] = targets["timepoint"] - pd.DateOffset(months=delta_months)
        scores = test(
            model_class,
            model_kwargs,
            tips,
            targets,
            test_timepoints,
            coefficients,
            group_by=group_by_attribute,
            include_attributes=args.include_attributes
        )
    else:
        # First, confirm that all predictors are defined in the given tip
        # attributes.
        if not all([predictor in tips.columns for predictor in args.predictors]):
            print("ERROR: Not all predictors could be found in the given tip attributes table.", file=sys.stderr)
            sys.exit(1)

        # Select the cost function.
        cost_function_name = args.cost_function
        cost_function = cost_functions_by_name[cost_function_name]

        # Identify train/validate splits from timepoints.
        training_window = args.training_window
        train_validate_timepoints = get_train_validate_timepoints(
            timepoints,
            args.delta_months,
            training_window
        )

        # For each train/validate split, fit a model to the training data, and
        # evaluate the model with the validation data, storing the training results,
        # beta parameters, and validation results.
        delta_months = args.delta_months
        delta_time = delta_months / 12.0
        l1_lambda = args.l1_lambda
        model_kwargs.update({
            "predictors": args.predictors,
            "delta_time": delta_time,
            "l1_lambda": l1_lambda,
            "cost_function": cost_function
        })

        # If this is a naive model, set the coefficients to zero so cross-validation
        # can run under naive model conditions.
        if "naive" in args.predictors:
            coefficients = np.zeros(len(args.predictors))
        else:
            coefficients = None

        targets["timepoint"] = targets["timepoint"] - pd.DateOffset(months=delta_months)
        scores = cross_validate(
            model_class,
            model_kwargs,
            tips,
            targets,
            train_validate_timepoints,
            coefficients,
            group_by=group_by_attribute,
            include_attributes=args.include_attributes
        )

    # Summarize model errors including in-sample errors by AIC, out-of-sample
    # errors by cross-validation, and beta parameters across timepoints.
    model_results = summarize_scores(scores, args.include_scores)

    # Annotate parameters used to produce models.
    model_results["cost_function"] = cost_function_name
    model_results["l1_lambda"] = l1_lambda
    model_results["delta_months"] = delta_months
    model_results["training_window"] = training_window
    model_results["pseudocount"] = args.pseudocount

    # Save model fitting hyperparameters, raw results, and summary of results to
    # JSON.
    with open(args.output, "w") as fh:
        json.dump(model_results, fh, indent=1)

    # Save errors by timepoint, if requested.
    if args.errors_by_timepoint:
        errors_by_timepoint_df = get_errors_by_timepoint(scores)
        errors_by_timepoint_df.to_csv(args.errors_by_timepoint, sep="\t", header=True, index=False)

    # Save coefficients by timepoint, if requested.
    if args.coefficients_by_timepoint:
        coefficients_by_timepoint_df = get_coefficients_by_timepoint(scores)
        coefficients_by_timepoint_df.to_csv(args.coefficients_by_timepoint, sep="\t", header=True, index=False)
