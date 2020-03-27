import argparse
from augur.frequency_estimators import logit_transform, TreeKdeFrequencies
from augur.lbi import select_nodes_in_season
from augur.utils import write_json
from collections import defaultdict
import numpy as np
import os
import pandas as pd
from scipy.interpolate import interp1d
from scipy.stats import linregress, pearsonr
import sys

from .fitness_predictors import fitness_predictors


min_tips = 10
pc=1e-2
regularization = 1e-3
default_predictors = ['lb', 'ep', 'ne_star']


def float_to_datestring(time):
    """Convert a floating point date to a date string
    """
    year = int(time)
    month = int(((time - year) * 12) + 1)
    day = 1
    return "-".join(map(str, (year, month, day)))

def timestamp_to_float(time):
    """Convert a pandas timestamp to a floating point date.
    """
    return time.year + ((time.month - 1) / 12.0)

def process_predictor_args(predictors, params=None, sds=None):
    """Returns a predictor data structure for the given lists of predictors, params,
    and standard deviations.

    When no parameters or deviations are provided, the predictors are a simple
    list. When parameters and deviations are provided, the predictor are a
    dictionary indexed by predictor name with values corresponding to each
    predictor's param and global standard deviation.

    >>> process_predictor_args(None, None, None)
    >>> process_predictor_args(['ep'])
    ['ep']
    >>> process_predictor_args(['ep'], None, None)
    ['ep']
    >>> process_predictor_args(['ep'], [1], [5])
    {'ep': [1, 5]}
    """
    if predictors is None:
        processed_predictors = None
    elif params is None or sds is None:
        processed_predictors = predictors
    else:
        merged_params = map(list, zip(params, sds))
        processed_predictors = dict(zip(predictors, merged_params))

    return processed_predictors


def make_pivots(start, stop, pivots_per_year=12, precision=2):
    """Makes an array of pivots (i.e., timepoints) between the given start and stop
    by the given pivots per year. The generated pivots are floating point values
    that are then rounded to the requested decimal precision.

    >>> list(make_pivots(2000.0, 2001.0, 5))
    [2000.0, 2000.25, 2000.5, 2000.75, 2001.0]
    """
    # Calculate number of pivots (i.e., months) in the requested interval.
    number_of_pivots = np.ceil((stop - start) * pivots_per_year)

    # Build an evenly-spaced closed interval (including the start and stop
    # points) based on the calculated number of pivots.
    return np.around(
        np.linspace(start, stop, number_of_pivots),
        precision
    )


def matthews_correlation_coefficient(tp, tn, fp, fn):
    """Return Matthews correlation coefficient for values from a confusion matrix.
    Implementation is based on the definition from wikipedia:

    https://en.wikipedia.org/wiki/Matthews_correlation_coefficient
    """
    numerator = (tp * tn) - (fp * fn)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if denominator == 0:
            denominator = 1

    return float(numerator) / denominator


def get_matthews_correlation_coefficient_for_data_frame(freq_df, return_confusion_matrix=False):
        """Calculate Matthew's correlation coefficient from a given pandas data frame
        with columns for initial, observed, and predicted frequencies.
        """
        observed_growth = (freq_df["observed_freq"] > freq_df["initial_freq"])
        predicted_growth = (freq_df["predicted_freq"] > freq_df["initial_freq"])
        true_positives = ((observed_growth) & (predicted_growth)).sum()
        false_positives= ((~observed_growth) & (predicted_growth)).sum()

        observed_decline = (freq_df["observed_freq"] <= freq_df["initial_freq"])
        predicted_decline = (freq_df["predicted_freq"] <= freq_df["initial_freq"])
        true_negatives = ((observed_decline) & (predicted_decline)).sum()
        false_negatives = ((~observed_decline) & (predicted_decline)).sum()

        mcc = matthews_correlation_coefficient(
            true_positives,
            true_negatives,
            false_positives,
            false_negatives
        )

        if return_confusion_matrix:
            confusion_matrix = {
                "tp": true_positives,
                "tn": true_negatives,
                "fp": false_positives,
                "fn": false_negatives
            }

            return mcc, confusion_matrix
        else:
            return mcc


def sum_of_squared_errors(observed_freq, predicted_freq):
    """
    Calculates the sum of squared errors for observed and predicted frequencies.

    Args:
        observed_freq (numpy.ndarray): observed frequencies
        predicted_freq (numpy.ndarray): predicted frequencies

    Returns:
        float: sum of squared errors between observed and predicted frequencies
    """
    return np.sum((observed_freq - predicted_freq) ** 2)


def mean_absolute_error(observed_freq, predicted_freq):
    """
    Calculates the mean absolute error between observed and predicted frequencies.

    Args:
        observed_freq (numpy.ndarray): observed frequencies
        predicted_freq (numpy.ndarray): predicted frequencies

    Returns:
        float: mean absolute error between observed and predicted frequencies
    """
    return np.mean(np.abs(observed_freq - predicted_freq))


def project_clade_frequencies_by_delta_from_time(tree, model, time, delta, delta_steps_per_year=12):
    """
    Project clade frequencies from a given time to the future by a given delta.
    """
    # Calculate the steps between the projection date and delta time into the
    # future. First, find the frequency pivot that is closest to the requested
    # projection date.
    max_date = model.timepoints[np.searchsorted(model.timepoints, time)]
    future_date = max_date + delta

    # Then, calculate a fixed number of steps between that pivot and delta time
    # into the future.
    projected_pivots = np.linspace(max_date, future_date, int(delta_steps_per_year * delta))
    deltas = projected_pivots - max_date

    # Identify tip predictors and frequencies at the current time point.
    all_pred = model.predictor_arrays[max_date]
    all_freqs = model.freq_arrays[max_date]

    # For each requested delta, project current tip frequencies using the model
    # and calculate the corresponding projected clade frequencies.
    projected_clade_frequencies = defaultdict(list)

    for delta in deltas:
        # Project all tip frequencies.
        pred_freq = model.projection(model.model_params, all_pred, all_freqs, delta)

        # Normalize projected frequencies.
        pred_freq = pred_freq / pred_freq.sum()

        # Store projected frequencies by clade id.
        for i, tip in enumerate(model.tips):
            projected_clade_frequencies[tip.name].append(pred_freq[i])

        # Calculate projected frequencies for internal nodes and store by clade it.
        for node in tree.find_clades(order="postorder"):
            if not node.is_terminal():
                projected_clade_frequencies[node.name].append(pred_freq[node.tips].sum())

    projected_frequencies = {
        "params": {
            "max_date": max_date
        },
        "data": {
            "pivots": projected_pivots.tolist(),
            "frequencies": projected_clade_frequencies
        }
    }

    return projected_frequencies


def get_train_validate_timepoints(timepoints, delta_time, training_window):
    """Return all possible train-validate timepoints from the given complete list of
    timepoints, a delta time to project forward by, and the required number of
    timepoints to include in each training window.

    Parameters
    ----------
    timepoints : list
        Date/time strings to use for model training and validation

    delta_time : int
        Number of months into the future that the model will project

    training_window : int
        Number of years to include in each training window

    Returns
    -------
    list
        List of dictionaries containing all possible train-validate timepoints indexed by "train" and "validate" keys
    """
    # Convert list of date/time strings into pandas datetimes.
    timepoints = pd.to_datetime(timepoints)

    # Convert delta time and training window to pandas offsets.
    delta_time = pd.DateOffset(months=delta_time)
    training_window = pd.DateOffset(years=training_window)

    # Filter timepoints to those with enough future years to train and project
    # from. This means timepoints must not extend beyond the last training
    # timepoint plus its delta time in the future and the delta time in the
    # future for the validation interval.
    #
    # If the timepoints range from October 2005 to October 2015, the training
    # window is 4 years, and the delta time is 1 year, then the last validation
    # timepoint is October 2014 and the last training timepoint is October
    # 2013. This allows the model to train from October 2013 to October 2014 and
    # then validate from October 2014 to 2015. Assuming there is enough data in
    # the tree up to October 2005, the earliest point we can project from is
    # then October 2008 such that the previous 4 years are included in that
    # projection.
    is_valid_projection_timepoint = (timepoints + training_window + delta_time + delta_time) <= timepoints[-1]
    projection_timepoints = timepoints[is_valid_projection_timepoint].copy()

    # Split valid timepoint index values into all possible train/test sets.
    train_validate_timepoints = []
    for start_timepoint in projection_timepoints:
        end_timepoint = start_timepoint + training_window
        train_timepoints = timepoints[
            (timepoints >= start_timepoint) &
            (timepoints <= end_timepoint)
        ]
        validate_timepoint = train_timepoints[-1] + delta_time

        # Store train/validate timepoints as datetime strings to enable
        # downstream use by other tools.
        train_validate_timepoints.append({
            "train": train_timepoints.strftime("%Y-%m-%d").tolist(),
            "validate": validate_timepoint.strftime("%Y-%m-%d")
        })

    return train_validate_timepoints


class fitness_model(object):

    def __init__(self, tree, frequencies, predictor_input, cross_validate=False, censor_frequencies=True,
                 pivot_spacing=1.0 / 12, verbose=0, enforce_positive_predictors=True, predictor_kwargs=None,
                 cost_function=sum_of_squared_errors, delta_time=1.0, end_date=None, step_size=0.5, min_tip_freq=1e-3,
                 min_training_window=4.0, **kwargs):
        """

        Args:
            tree (Bio.Phylo): an annotated tree for which a fitness model is to be determined
            frequencies (TreeKdeFrequencies): a frequency estimator and its parameters
            predictor_input: a list of predictors to fit or dict of predictors to coefficients / std deviations
            censor_frequencies (bool): whether frequencies should censor future data or not
            pivot_spacing:
            verbose:
            enforce_positive_predictors:
            predictor_kwargs:
            cost_function (callable): a function that takes observed and predicted frequencies and returns a single error value
            **kwargs:
        """
        self.tree = tree
        self.frequencies = frequencies
        self.cross_validate = cross_validate
        self.censor_frequencies = censor_frequencies
        self.pivot_spacing = pivot_spacing
        self.verbose = verbose
        self.enforce_positive_predictors = enforce_positive_predictors
        self.estimate_coefficients = True
        self.min_freq = kwargs.get("min_freq", 0.1)
        self.max_freq = kwargs.get("max_freq", 0.99)
        self.min_tip_freq = min_tip_freq
        self.cost_function = cost_function
        self.end_date = end_date
        self.min_training_window = min_training_window

        if predictor_kwargs is None:
            self.predictor_kwargs = {}
        else:
            self.predictor_kwargs = predictor_kwargs

        self.time_window = kwargs.get("time_window", 6.0 / 12.0)

        if isinstance(predictor_input, dict):
            predictor_names = list(predictor_input.keys())
            self.estimate_coefficients = False
        else:
            predictor_names = predictor_input
        if "estimate_fitness_model" in kwargs:
            if kwargs["estimate_fitness_model"]:
                self.estimate_coefficients = True

        # Reestimate frequencies if they have not already been estimated or if internal nodes were excluded.
        if not hasattr(self.frequencies, "frequencies") or not self.frequencies.include_internal_nodes:
            sys.stderr.write("Recalculating frequencies\n")
            frequency_params = self.frequencies.get_params()
            frequency_params["include_internal_nodes"] = True
            self.frequencies = TreeKdeFrequencies(**frequency_params)
            self.frequencies.estimate(self.tree)

        # Pivots should be defined by frequencies.
        self.pivots = self.frequencies.pivots

        # final timepoint is end of interval and is only projected forward, not tested
        self.time_interval = (self.frequencies.start_date, self.frequencies.end_date)
        self.timepoint_step_size = step_size # amount of time between timepoints chosen for fitting
        self.delta_time = delta_time         # amount of time projected forward to do fitting

        # Convert delta time and step size to month integers from floating points.
        months_step_size = int(12 * step_size)
        delta_time_months = pd.DateOffset(months=int(delta_time * 12))

        # Calculate timepoints with a date/time-aware range method.
        # Exclude the first timepoint when there are no viruses sampled yet.
        # Somewhat counterintuitively, this requires the closed attribute to be "right", below.
        timepoints = pd.date_range(
            float_to_datestring(self.frequencies.start_date),
            float_to_datestring(self.frequencies.end_date),
            freq="%sMS" % months_step_size,
            closed="right"
        )

        # Convert date range to floats for consistency with downstream functions.
        timepoints = [timestamp_to_float(timepoint) for timepoint in timepoints]
        self.timepoints = np.array(timepoints)

        # If an end date was provided, exclude all timepoints after that date.
        if self.end_date:
            original_number_of_timepoints = len(self.timepoints)
            self.timepoints = [time for time in self.timepoints
                               if time < self.end_date]
            filtered_number_of_timepoints = len(self.timepoints)
            sys.stderr.write("Fit model up to %s (filtered from %s to %s timepoints)\n" % (self.end_date, original_number_of_timepoints, filtered_number_of_timepoints))

        # Determine the set of timepoints the model can project from. This will
        # filter out timepoints for which we need to estimate frequencies
        # because they are an endpoint of a projection but which cannot be used
        # themselves to project forward by delta time.
        self.projection_timepoints = np.array([
            timepoint
            for timepoint in self.timepoints
            if timepoint + self.delta_time <= self.timepoints[-1]
        ])

        self.predictors = predictor_names

        self.model_params = np.zeros(len(self.predictors))
        if isinstance(predictor_input, dict):
            self.model_params = np.array([predictor_input[k][0] for k in predictor_names])

        self.to_standardize = np.array([p!='dfreq' for p in self.predictors])
        if isinstance(predictor_input, dict):
            self.global_sds = np.array([predictor_input[k][1] for k in predictor_names])
        else:
            self.global_sds = np.zeros(len(self.predictors))

        self.fp = fitness_predictors(predictor_names = predictor_names, **kwargs)

        # Map node names to parents.
        self.node_parents = {}
        for clade in self.tree.find_clades(order='level'):
            for child in clade:
                self.node_parents[child] = clade

    def prep_nodes(self):
        """Assigns data from the tree to top-level fitness model attributes.

        TODO: consider moving this code directly into the `predict`
        method since it is only ever called there.
        """
        self.nodes = [node for node in self.tree.find_clades(order="postorder")]
        self.tips = [node for node in self.nodes if node.is_terminal()]
        self.rootnode = self.tree.root
        self.rootnode.pivots = self.pivots

        # Create a list of tip indices under node.tips that map to self.tips
        # list.
        tip_index_region_specific = 0
        for node in self.nodes:
            tmp_tips = []
            if node.is_terminal():
                tmp_tips.append((tip_index_region_specific, node.attr["num_date"]))
                tip_index_region_specific += 1

            for child in node.clades:
                tmp_tips.extend(child.tips)

            # Sort tips by corresponding date.
            node.tips = np.array([x for x in sorted(tmp_tips, key = lambda x: x[1])])

        # Erase the dates from the tip lists and cast to int such that they can
        # be used for indexing. These operations must happen after all nodes
        # have been processed and sorted.
        for node in self.nodes:
            if len(node.tips.shape) == 2:
                node.tips = np.array(node.tips[:, 0], dtype=int)
            else:
                node.tips = np.array([], dtype=int)

    def calc_node_frequencies(self):
        '''
        goes over all nodes and calculates frequencies at timepoints based on previously calculated frequency trajectories
        '''
        region = "global"
        frequencies = self.frequencies.frequencies

        # Annotate frequencies on nodes using all available data regardless of tip frequency censoring status.
        for node in self.nodes:
            node.freq = {
                region: frequencies[node.name]
            }
            node.logit_freq = {
                region: logit_transform(frequencies[node.name], 1e-4)
            }

        for node in self.nodes:
            interpolation = interp1d(self.rootnode.pivots, node.freq[region], kind='linear', bounds_error=True)
            node.timepoint_freqs = {}
            for time in self.timepoints:
                node.timepoint_freqs[time] = np.asscalar(interpolation(time))

        # Estimate frequencies for tips at specific timepoints.
        # Censor future tips from estimations unless these data are explicitly allowed.
        # freq_arrays list *all* tips for each initial timepoint
        self.freq_arrays={}
        frequency_parameters = self.frequencies.get_params()
        for time in self.timepoints:
            tmp_freqs = []

            if self.censor_frequencies and not frequency_parameters["censored"]:
                # Censor frequencies by omitting all tips sampled after the current timepoint.
                if self.verbose:
                    sys.stderr.write("Calculating censored frequencies for %s\n" % time)

                frequency_parameters["max_date"] = time
                frequency_estimator = TreeKdeFrequencies(**frequency_parameters)
                frequencies = frequency_estimator.estimate(self.tree)

                # Determine the frequency of each tip at the given timepoint.
                total_freq = 0.0
                for tip in self.tips:
                    interpolation = interp1d(
                        self.pivots,
                        frequencies[tip.name],
                        kind="linear",
                        bounds_error=True
                    )
                    censored_frequency = np.asscalar(interpolation(time))
                    total_freq += censored_frequency

                    if not hasattr(tip, "censored_freqs"):
                        tip.censored_freqs = {}

                    tip.censored_freqs[time] = censored_frequency
                    tmp_freqs.append(censored_frequency)

                # Normalize tip frequencies interpolated from pivots at timepoints.
                if total_freq > 0:
                    tmp_freqs = [tmp_freq / total_freq for tmp_freq in tmp_freqs]

                    for tip in self.tips:
                        tip.censored_freqs[time] /= total_freq
            else:
                for tip in self.tips:
                    tip.censored_freqs = tip.timepoint_freqs
                    tmp_freqs.append(tip.timepoint_freqs[time])

            self.freq_arrays[time] = np.array(tmp_freqs)

        # Assign observed final frequencies to clades based on the sum of their respective tip frequencies.
        # This allows clade frequencies to be censored at future timepoints if tip frequencies are too.
        for node in self.nodes:
            node.observed_final_freqs = {}
            for time in self.projection_timepoints:
                node.observed_final_freqs[time] = self.freq_arrays[time + self.delta_time][node.tips].sum(axis=0)

    def calc_predictors(self, timepoint):
        for pred in self.predictors:
            # calculate the predictors for all nodes of the tree and save as node.attr
            if pred != 'dfreq':
                self.fp.setup_predictor(self.tree, pred, timepoint, **self.predictor_kwargs)

    def calc_time_censored_tree_frequencies(self):
        print("fitting time censored tree frequencies")
        # this doesn't interfere with the previous freq estimates via difference in region: global_censored vs global
        region = "global_censored"
        frequency_parameters = self.frequencies.get_params()
        freq_cutoff = 25.0
        pivots_fit = 6
        freq_window = 1.0
        for node in self.nodes:
            node.fit_frequencies = {}
            node.freq_slope = {}
        for time in self.timepoints:
            time_interval = [time - freq_window, time]
            node_filter_func = lambda node: node.attr['num_date'] >= time_interval[0] and node.attr['num_date'] < time_interval[1]

            # Recalculate  frequencies for the given time interval and its corresponding pivots.
            frequency_parameters["start_date"] = time_interval[0]
            frequency_parameters["end_date"] = time_interval[1]
            frequency_parameters["max_date"] = time_interval[1]
            frequency_estimator = TreeKdeFrequencies(**frequency_parameters)
            frequencies = frequency_estimator.estimate(self.tree)

            # Annotate censored frequencies on nodes.
            # TODO: replace node-based annotation with dicts indexed by node name.
            for node in self.nodes:
                node.freq = {
                    region: frequencies[node.name]
                }
                node.logit_freq = {
                    region: logit_transform(frequencies[node.name], 1e-4)
                }

            for node in self.nodes:
                if node.logit_freq[region] is not None:
                    node.fit_frequencies[time] = np.minimum(freq_cutoff, np.maximum(-freq_cutoff,node.logit_freq[region]))
                else:
                    node.fit_frequencies[time] = self.node_parents[node].fit_frequencies[time]
                try:
                    slope, intercept, rval, pval, stderr = linregress(self.pivots[pivots_fit:-1], node.fit_frequencies[time][pivots_fit:-1])
                    node.freq_slope[time] = slope
                except:
                    import ipdb; ipdb.set_trace()

        # reset pivots in tree to global pivots
        self.rootnode.pivots = self.pivots


    def calc_all_predictors(self, estimate_frequencies = True):
        if estimate_frequencies and 'dfreq' in [x for x in self.predictors]:
            self.calc_time_censored_tree_frequencies()
        # predictor_arrays list *all* tips for each timepoint
        self.predictor_arrays={}
        for node in self.nodes:
            node.predictors = {}
        for time in self.timepoints:
            if self.verbose: print("calculating predictors for time", time)
            select_nodes_in_season(self.tree, time, self.time_window)
            self.calc_predictors(time)

            for node in self.tips:
                if 'dfreq' in [x for x in self.predictors]: node.dfreq = node.freq_slope[time]

                predictors_at_time = []
                for pred in self.predictors:
                    if hasattr(node, pred):
                        predictors_at_time.append(getattr(node, pred))
                    else:
                        predictors_at_time.append(node.attr[pred])

                node.predictors[time] = np.array(predictors_at_time)

            tmp_preds = []
            for tip in self.tips:
                tmp_preds.append(tip.predictors[time])
            self.predictor_arrays[time]=np.array(tmp_preds, dtype=float)

    def standardize_predictors(self):
        self.predictor_means = {}
        self.predictor_sds = {}
        if self.verbose: print("standardizing predictors")
        for time in self.timepoints:
            values = self.predictor_arrays[time]
            weights = self.freq_arrays[time]

            # Do not calculate weighted summary statistics when all frequencies at the current timepoint are zero.
            if weights.sum() == 0:
                weights = None

            means = np.average(values, weights=weights, axis=0)
            variances = np.average((values-means)**2, weights=weights, axis=0)
            sds = np.sqrt(variances)
            self.predictor_means[time] = means
            self.predictor_sds[time] = sds

        if self.estimate_coefficients:
            self.global_sds = np.mean(list(self.predictor_sds.values()), axis=0)

        for time in self.timepoints:
            for node in self.tips:
                if node.predictors[time] is not None and (self.global_sds == 0).sum() == 0:
                    node.predictors[time] = (node.predictors[time]-self.predictor_means[time]) / self.global_sds
            self.predictor_arrays[time][:,self.to_standardize] -= self.predictor_means[time][self.to_standardize]

            if (self.global_sds == 0).sum() == 0:
                self.predictor_arrays[time][:,self.to_standardize] /= self.global_sds[self.to_standardize]

    def select_clades_for_fitting(self):
        # for each time point, select clades that are within the specified frequency window
        # keep track in the dict fit_clades that maps timepoint to clade list
        self.fit_clades = {}
        for time in self.projection_timepoints:
            self.fit_clades[time] = []
            for node in self.nodes:
                # Only select clades for fitting if their censored frequencies are within the specified thresholds.
                node_freq = self.freq_arrays[time][node.tips].sum(axis=0)

                if self.min_freq <= node_freq <= self.max_freq:
                    # Exclude subclades whose frequency is identical to their parent clade.
                    parent_node_freq = self.freq_arrays[time][self.node_parents[node].tips].sum(axis=0)
                    if node_freq < parent_node_freq:
                        self.fit_clades[time].append(node)

    def select_nonoverlapping_clades_for_fitting(self):
        """For each timepoint, identify clades that originated at or after the previous
        timepoint and whose frequencies at the current timepoint are within a
        specified range.
        """
        self.fit_clades = {}

        for timepoint in self.projection_timepoints:
            previous_timepoint = timepoint - self.timepoint_step_size
            total_freq = []
            candidate_clades = []
            self.fit_clades[timepoint] = []

            tips_in_window = {}
            for node in self.tree.find_clades():
                # Select all tips that are likely sampled from this timepoint
                # based on a minimum frequency.
                if node.is_terminal() and node.censored_freqs[timepoint] > self.min_tip_freq:
                    tips_in_window[node.name] = node

            if self.verbose:
                sys.stderr.write("Tips in window sum to frequency of %.2f\n" % sum([node.censored_freqs[timepoint]
                                                                                    for node in tips_in_window.values()]))

            # Order tips by the date of their parents.
            # This preferentially adds parents that represent more tips.
            tip_names = [key for key, value in sorted(tips_in_window.items(), key=lambda item: item[1].parent.attr["num_date"])]

            for tip_name in tip_names:
                if tip_name in tips_in_window:
                    tip = tips_in_window[tip_name]

                    # Test whether the tip's parent contains any other tips in this window.
                    parent = tip.parent
                    while parent:
                        tips_added_by_parent = 0
                        for child in parent.get_terminals():
                            if child.name in tips_in_window:
                                # Note that this parent accounts for other tips in this window.
                                tips_added_by_parent += 1
                                del tips_in_window[child.name]

                        # Stop adding parents, if doing so provides no new tips.
                        # Save the current parent, otherwise.
                        if tips_added_by_parent == 0:
                            break
                        else:
                            candidate_clades.append(parent)

                        # Stop adding parents, if we have moved outside of the current window.
                        # Otherwise, keep traversing upward.
                        if parent.attr["num_date"] < previous_timepoint:
                            break
                        else:
                            parent = parent.parent

            # Reset attribute on nodes.
            for clade in self.tree.find_clades():
                if "clade_group" in clade.attr:
                    del clade.attr["clade_group"]

            clade_group = 0
            clade_group_freq = []
            nested = 0
            for clade in self.tree.find_clades():
                if clade.parent and "clade_group" in clade.parent.attr:
                    clade.attr["clade_group"] = clade.parent.attr["clade_group"]
                    nested += 1
                elif clade in candidate_clades:
                    clade_group += 1
                    clade.attr["clade_group"] = clade_group

                    # Filter clades by minimum and maximum summed censored
                    # frequency of their tips at the current timepoint.
                    node_freq = self.freq_arrays[timepoint][clade.tips].sum(axis=0)
                    clade_group_freq.append(node_freq)
                    if self.min_freq <= node_freq <= self.max_freq:
                        total_freq.append(node_freq)
                        self.fit_clades[timepoint].append(clade)

            if self.verbose:
                sys.stderr.write("%s: %s clades totalling %.2f (%s nested, %s clade groups at %.2f)\n" % (timepoint, len(self.fit_clades[timepoint]), sum(total_freq), nested, clade_group, sum(clade_group_freq)))
                sys.stderr.write("%s\n" % sorted(clade_group_freq))

    def clade_fit(self, params, test=False):
        # walk through initial/final timepoint pairs
        # tested that the sum of frequencies of tips within a clade is equal to the direct clade frequency
        timepoint_errors = []
        self.pred_vs_true = []
        pred_vs_true_values = []

        if test:
            timepoints = self.test_timepoints
            print("Testing fit parameters")
        else:
            timepoints = self.train_timepoints

        for time in timepoints:
            # Project all tip frequencies forward by the specific delta time.
            all_pred_freq = self.projection(params, self.predictor_arrays[time], self.freq_arrays[time], self.delta_time)
            assert all_pred_freq.shape == self.freq_arrays[time].shape

            # normalization factor for predicted tip frequencies
            total_pred_freq = np.sum(all_pred_freq)

            # project clades forward according to strain makeup
            clade_errors = []
            tmp_pred_vs_true = []
            for clade in self.fit_clades[time]:
                # The observed final frequency is calculated for each clade from all available data.
                obs_final_freq = clade.observed_final_freqs[time]

                # The initial frequency is calculated from the sum of each clade's censored tip frequencies.
                initial_freq = self.freq_arrays[time][clade.tips].sum(axis=0)

                # The predicted final frequency is also calculated from each clade's censored tip frequencies modified
                # by the fitness and model parameters.
                pred_final_freq = np.sum(all_pred_freq[clade.tips]) / total_pred_freq

                tmp_pred_vs_true.append((initial_freq, obs_final_freq, pred_final_freq))
                pred_vs_true_values.append((time, time + self.delta_time, clade.name, len(clade.tips), initial_freq, obs_final_freq, pred_final_freq))

            self.pred_vs_true.append(np.array(tmp_pred_vs_true))

        # Prepare a data frame with all initial, observed, and predicted frequencies by time and clade.
        self.pred_vs_true_df = pd.DataFrame(
            pred_vs_true_values,
            columns=("timepoint", "projected_timepoint", "clade", "clade_size", "initial_freq", "observed_freq", "predicted_freq")
        )

        training_error = self.cost_function(
            self.pred_vs_true_df["observed_freq"],
            self.pred_vs_true_df["predicted_freq"]
        )
        if np.isnan(training_error) or np.isinf(training_error):
            training_error = 1e10

        self.last_fit = training_error
        if self.verbose>2: print(params, self.last_fit)
        penalty = regularization*np.sum(params**2)
        if self.enforce_positive_predictors:
            for param in params:
                if param < 0:
                    penalty += 1
        return training_error + penalty

    def weighted_af(self, seqs, weights):
        af = np.zeros((4, seqs.shape[1]))
        for ni, nuc in enumerate('ACGT'):
            af[ni] += (weights*(seqs==nuc).T).sum(axis=1)/weights.sum()
        return af

    def af_fit(self, params):
        # TODO: fix me for continuous prediction
        seasonal_errors = []
        self.pred_vs_true = []
        for s,t in self.fit_test_season_pairs:
            weights = np.exp(self.fitness(params, self.predictor_arrays[s][self.tree.root.season_tips[s],:]))
            pred_af = self.weighted_af(self.seqs[s],weights)
            #seasonal_errors.append(np.mean(np.sum((pred_af-self.af[t])**2, axis=0), axis=0))
            future_diameter = 0.5*np.sum(np.sum(self.af[t]*(1-self.af[t]), axis=0), axis=0)
            seasonal_errors.append(np.sum(np.sum(pred_af*(1-self.af[t]), axis=0), axis=0)-future_diameter)
            good_ind = self.af[s]*(1-self.af[s])>0.05
            self.pred_vs_true.append(np.array(zip(self.af[s][good_ind], self.af[t][good_ind], pred_af[good_ind])))

        mean_error = np.mean(seasonal_errors)
        if any(np.isnan(seasonal_errors)+np.isinf(seasonal_errors)):
            mean_error = 1e10
        self.last_fit = mean_error
        if self.verbose>2: print(params, self.last_fit)
        return mean_error + regularization*np.sum(params**2)

    def fitness(self, params, pred):
        return np.sum(params*pred, axis=-1)

    def projection(self, params, pred, freqs, delta):
        return freqs * np.exp(self.fitness(params, pred) * delta);

    def minimize_clade_error(self):
        from scipy.optimize import fmin as minimizer
        if self.verbose:
            print("initial function value:", self.clade_fit(self.model_params))
            print("initial parameters:", self.model_params)
        self.model_params = minimizer(self.clade_fit, self.model_params, disp = self.verbose>1)
        if self.verbose:
            print("final function value:", self.clade_fit(self.model_params))
            print("final parameters:", self.model_params, '\n')

    def prep_af(self):
        if not hasattr(self,'variable_nuc'):
            self.determine_variable_positions()
        fit_aln = np.zeros((len(self.tips), len(self.variable_nuc)), dtype='S1')
        for i in range(len(self.tips)):
            tip = self.tips[i]
            fit_aln[i] = np.fromstring(tip.seq, 'S1')[self.variable_nuc]
        self.seqs = fit_aln
        self.af = {}
        for time in self.timepoints:
            self.af[time] = self.weighted_af(self.seqs, self.freq_arrays[time])

    def minimize_af_error(self):
        from scipy.optimize import fmin as minimizer
        if self.verbose:
            print("initial function value:", self.af_fit(self.model_params))
            print("initial parameters:", self.model_params)
        self.model_params = minimizer(self.af_fit, self.model_params, disp = self.verbose>1)
        if self.verbose:
            print("final function value:", self.af_fit(self.model_params))
            print("final parameters:", self.model_params, '\n')


    def learn_parameters(self, niter = 10, fit_func = "clade"):
        if fit_func=='clade':
            minimize_error=self.minimize_clade_error
            fit_func=self.clade_fit
        elif fit_func=="af":
            minimize_error=self.minimize_af_error
            fit_func=self.af_fit
        else:
            print("fit function", fit_func,"does not exist")
            raise NotImplementedError

        print("fitting parameters of the fitness model\n")

        params_stack = []

        if self.verbose:
            print("null parameters")
        self.model_params = 0*np.ones(len(self.predictors))  # initial values
        minimize_error()
        params_stack.append((self.last_fit, self.model_params))

        for ii in range(niter):
            if self.verbose:
                print("iteration:", ii+1)
            self.model_params = np.random.rand(len(self.predictors)) #0*np.ones(len(self.predictors))  # initial values
            minimize_error()
            params_stack.append((self.last_fit, self.model_params))

        self.model_params = params_stack[np.argmin([x[0] for x in params_stack])][1]
        fit_func(self.model_params)
        if self.verbose:
            print("best after",niter,"iterations\nfunction value:", self.last_fit)
            print("fit parameters:")
            for pred, val in zip(self.predictors, self.model_params):
                print(pred,':', val)


    def assign_fitness(self):
        if self.verbose: print("calculating predictors for the final timepoint")
        final_timepoint = self.timepoints[-1]

        for node in self.tips:
            if node.predictors[final_timepoint] is not None:
                node.fitness = self.fitness(self.model_params, node.predictors[final_timepoint])
            else:
                node.fitness = 0.0

            node.attr["fitness"] = node.fitness

    def assign_predicted_frequency(self, delta=1.0):
        total_freq = 0
        timepoint = self.timepoints[-1]
        for node in self.tree.get_terminals():
            pred = self.predictor_arrays[timepoint][node.tips]
            freqs = self.freq_arrays[timepoint][node.tips]
            node.predicted_freq = self.projection(self.model_params, pred, freqs, delta)[0]
            total_freq += node.predicted_freq

        for node in self.tree.get_terminals():
            node.predicted_freq /= total_freq
            node.attr["predicted_freq"] = node.predicted_freq

    def split_timepoints(self):
        """Split timepoints into train/test sets for cross-validation.
        """
        # Split valid timepoint index values into all possible train/test sets.
        train_test_splits = []
        total_timepoints = len(self.projection_timepoints)

        for i in range(total_timepoints):
            # Skip training intervals that are smaller than the minimum.
            if self.projection_timepoints[i] - self.projection_timepoints[0] < self.min_training_window:
                continue

            train = np.arange(i + 1)

            # Since timepoints occur more frequently than the prediction interval (delta t),
            # the test index for any given set of training indices must be at least delta t
            # into the future from the last training index to prevent training the model
            # on data overlapping with the test time intervals.
            found_test_index = False
            for test_index in range(i + 1, total_timepoints):
                if self.projection_timepoints[train[-1]] + self.delta_time <= self.projection_timepoints[test_index]:
                    found_test_index = True
                    break

            # If we found a valid test index, add to the existing train/test splits.
            # Otherwise, we have no more timepoints to consider.
            if found_test_index:
                test = np.array([test_index])
                train_test_splits.append([train, test])
            else:
                break

        assert len(train_test_splits) > 0, "No train/test timepoints found; consider changing the value of `min_training_window` from %s." % self.min_training_window
        return train_test_splits

    def predict(self, niter = 10, estimate_frequencies = True):
        self.prep_nodes()
        self.calc_node_frequencies()
        self.calc_all_predictors(estimate_frequencies = estimate_frequencies)
        self.standardize_predictors()
        self.select_nonoverlapping_clades_for_fitting()
        #self.select_clades_for_fitting()

        validation_df = None
        if self.estimate_coefficients:
            # If cross-validation was requested, identify train/test timepoints.
            # Otherwise, just learn parameters for all timepoints.
            if self.cross_validate:
                train_test_splits = self.split_timepoints()
                results = []
                for train, test in train_test_splits:
                    self.train_timepoints = self.projection_timepoints[train]
                    self.test_timepoints = self.projection_timepoints[test]

                    clades_for_training = sum([len(self.fit_clades[timepoint])
                                               for timepoint in self.train_timepoints])
                    clades_for_testing = sum([len(self.fit_clades[timepoint])
                                              for timepoint in self.test_timepoints])

                    if clades_for_training == 0 or clades_for_testing == 0:
                        sys.stderr.write(f"WARNING: No clades for training ({self.train_timepoints}) or testing ({self.test_timepoints}), skipping\n")
                        continue

                    print("train: %s, test: %s" % (self.train_timepoints, self.test_timepoints))
                    self.learn_parameters(niter = niter, fit_func = "clade")
                    training_error = self.clade_fit(self.model_params, test=False)
                    training_correlation = self.get_correlation()
                    training_mcc = get_matthews_correlation_coefficient_for_data_frame(self.pred_vs_true_df)

                    testing_error = self.clade_fit(self.model_params, test=True)
                    testing_correlation = self.get_correlation()
                    testing_mcc, testing_matrix = get_matthews_correlation_coefficient_for_data_frame(self.pred_vs_true_df, return_confusion_matrix=True)

                    result = {
                        "last_training_timepoint": self.projection_timepoints[train[-1]],
                        "test_timepoint": self.projection_timepoints[test[0]],
                        "training_windows": len(train),
                        "training_accuracy": training_mcc,
                        "training_correlation": training_correlation[-1][0],
                        "testing_accuracy": testing_mcc,
                        "testing_correlation": testing_correlation[-1][0],
                        "mae": mean_absolute_error(self.pred_vs_true_df["observed_freq"], self.pred_vs_true_df["predicted_freq"])
                    }
                    for predictor, parameter in zip(self.predictors, self.model_params):
                        result["parameter-%s" % predictor] = parameter

                    result.update(testing_matrix)
                    results.append(result)
                    print("train error: %s, test error: %s" % (training_error, testing_error))

                validation_df = pd.DataFrame(results)
                validation_df["predictors"] = "-".join(self.predictors)
            else:
                self.train_timepoints = self.projection_timepoints
                self.learn_parameters(niter = niter, fit_func = "clade")
        else:
            self.train_timepoints = self.projection_timepoints

        self.assign_fitness()
        self.assign_predicted_frequency()

        return validation_df

    def get_correlation(self):
        rho_null = pearsonr(self.pred_vs_true_df["initial_freq"], self.pred_vs_true_df["observed_freq"])
        rho_raw = pearsonr(self.pred_vs_true_df["observed_freq"], self.pred_vs_true_df["predicted_freq"])
        rho_rel = pearsonr((self.pred_vs_true_df["observed_freq"] / self.pred_vs_true_df["initial_freq"]),
                           (self.pred_vs_true_df["predicted_freq"] / self.pred_vs_true_df["initial_freq"]))

        return rho_null, rho_raw, rho_rel

    def validate_prediction(self, plot=False, test=False):
        test = test and self.cross_validate
        abs_clade_error = self.clade_fit(self.model_params, test=test)

        if plot:
            import matplotlib.pyplot as plt

            fig, axs = plt.subplots(1,4, figsize=(10,5))
            for time, pred_vs_true in zip(self.projection_timepoints, self.pred_vs_true):
                # 0: initial, 1: observed, 2: predicted
                axs[0].scatter(pred_vs_true[:,1], pred_vs_true[:,2])
                axs[1].scatter(pred_vs_true[:,1]/pred_vs_true[:,0],
                               pred_vs_true[:,2]/pred_vs_true[:,0], c=pred_vs_true[0])
                for s, o, p  in pred_vs_true:
                    axs[2].arrow(s, s, o-s, p-s)
                axs[3].scatter(pred_vs_true[:,0],
                               (pred_vs_true[:,2]+0.01)/(pred_vs_true[:,1]+0.01))

            axs[0].set_ylabel('predicted')
            axs[0].set_xlabel('observed')
            axs[1].set_ylabel('predicted/initial')
            axs[1].set_xlabel('observed/initial')
            axs[1].set_yscale('linear')
            axs[1].set_xscale('linear')
            axs[2].set_ylabel('predicted')
            axs[2].set_xlabel('observed')
            axs[2].set_ylim(-0.1, 1.1)
            axs[2].set_xlim(-0.1, 1.1)
            axs[3].set_ylabel('predicted / observed')
            axs[3].set_xlabel('initial')
            axs[3].set_yscale('log')

        print("Abs clade error:", abs_clade_error)

        rho_null, rho_raw, rho_rel = self.get_correlation()
        print("Pearson's R, null:", rho_null)
        print("Pearson's R, raw:", rho_raw)
        print("Pearson's R, rel:", rho_rel)

        trajectory_mcc, confusion_matrix = get_matthews_correlation_coefficient_for_data_frame(
            self.pred_vs_true_df,
            return_confusion_matrix=True
        )
        correct_growth = confusion_matrix["tp"]
        correct_decline = confusion_matrix["tn"]
        total_growth = float(correct_growth + confusion_matrix["fn"])
        total_decline = float(correct_decline + confusion_matrix["fp"])

        print("Correct at predicting growth: %s (%s / %s)" % ((correct_growth / total_growth), correct_growth, total_growth))
        print("Correct at predicting decline: %s (%s / %s)" % ((correct_decline / total_decline), correct_decline, total_decline))
        print("Correct classification:",  (correct_growth+correct_decline) / (total_growth+total_decline))
        print("Matthew's correlation coefficient: %s" % trajectory_mcc)
        print("Params:")
        print(list(zip(self.predictors, np.around(self.model_params, 4))))

        pred_data = []
        for time, pred_vs_true in zip(self.projection_timepoints, self.pred_vs_true):
            for entry in pred_vs_true:
                pred_data.append(np.append(entry, time))
        pred_vs_true_df = pd.DataFrame(pred_data, columns=['initial', 'obs', 'pred', 'time'])

        output_dir = "data"
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        pred_vs_true_df.to_csv(os.path.join(output_dir, "prediction_pairs.tsv"), sep="\t", index=False)

    def validate_trajectories(self):
        '''
        Project clade trajectories based on fitted fitness model and compare to observed trajectories
        '''
        self.trajectory_data = []
        series = 0
        for time in self.projection_timepoints:
            all_pred = self.predictor_arrays[time]
            all_freqs = self.freq_arrays[time]
            for clade in self.fit_clades[time]:
                initial_freq = clade.timepoint_freqs[time]
                pred = all_pred[clade.tips]
                freqs = all_freqs[clade.tips]
                interpolation = interp1d(self.rootnode.pivots, clade.freq['global'], kind='linear', bounds_error=True)
                for delta in np.arange(-0.5, 1.1, 0.1):
                    if time + delta >= self.rootnode.pivots[0] and time + delta <= self.rootnode.pivots[-1]:
                        obs_freq = np.asscalar(interpolation(time+delta))
                        pred_freq = obs_freq
                        if delta >= 0:
                            total_pred_freq = np.sum(self.projection(self.model_params, all_pred, all_freqs, delta))
                            pred_freq = np.sum(self.projection(self.model_params, pred, freqs, delta)) / total_pred_freq
                        self.trajectory_data.append([series, str(clade), time, time+delta, obs_freq, pred_freq])
                series += 1

        self.trajectory_data_df = pd.DataFrame(self.trajectory_data, columns=['series', 'clade', 'initial_time', 'time', 'obs', 'pred'])
        self.trajectory_data_df.to_csv("data/prediction_trajectories.tsv", sep="\t", index=False)

        import seaborn as sns
        import matplotlib.pyplot as plt
        cols = sns.color_palette(n_colors=6)
        fig, axs = plt.subplots(6,4, sharey=True)
        for tp, ax in zip(self.projection_timepoints, axs.flatten()):
            traj = self.trajectory_data_df[self.trajectory_data_df.initial_time == tp]
            clades = np.unique(traj['series'])
            for ci in clades:
                tmp = traj[traj['series']==ci]
                ax.plot(tmp['time'], tmp['obs'], ls='-', c=cols[ci%6])
                ax.plot(tmp['time'], tmp['pred'], ls='--', c=cols[ci%6])

    def to_data_frames(self):
        """Return data frames representing the fitness model's inputs for each tip and clade, respectively.

        Inputs include tip metadata, frequencies, and standardized predictor values.
        """
        tip_records = []
        clade_records = []

        # Include only timepoints used to fit the model itself (excluding the last timepoint).
        for timepoint in self.projection_timepoints:
            # Create a record for each clade fit by the model despite nesting of clades.
            for clade in self.fit_clades[timepoint]:
                # Store information for each tip in the current clade despite
                # redundancy of information. This enables refitting the model
                # with the same data later.
                total_clade_tips = 0
                for tip_index in clade.tips:
                    if self.freq_arrays[timepoint][tip_index] > 0.0:
                        total_clade_tips += 1
                        tip = self.tips[tip_index]
                        tip_record = {
                            "timepoint": timepoint,
                            "clade_name": clade.name,
                            "name": tip.name,
                            "num_date": tip.attr["num_date"],
                            "observed_frequency": tip.timepoint_freqs[timepoint],
                            "censored_frequency": self.freq_arrays[timepoint][tip_index]
                        }

                        # Store standardized predictor values using a column name
                        # prefix that enables downstream analyses to easily identify
                        # predictor columns.
                        for predictor_index, predictor in enumerate(self.predictors):
                            tip_record["predictor:%s" % predictor] = self.predictor_arrays[timepoint][tip_index][predictor_index]

                        tip_records.append(tip_record)

                clade_record = {
                    "timepoint": timepoint,
                    "clade_name": clade.name,
                    "initial_frequency": self.freq_arrays[timepoint][clade.tips].sum(axis=0),
                    "final_frequency": clade.observed_final_freqs[timepoint],
                    "total_tips": total_clade_tips
                }
                clade_records.append(clade_record)

        return pd.DataFrame(tip_records), pd.DataFrame(clade_records)

    @classmethod
    def from_json(cls, tree, frequencies, json_dict):
        """Return an instance of the fitness model defined by the given tree, frequencies, and model JSON.

        The JSON input is expected to be the output of the `to_json` method of the fitness model class.
        """
        predictors = {record["predictor"]: [round(record["param"], 2), round(record["global_sd"], 2)]
                      for record in json_dict["params"]}
        predictors_key = "-".join(sorted([record["predictor"] for record in json_dict["params"]]))
        predictor_kwargs = json_dict["predictor_kwargs"]

        model = cls(
            tree,
            frequencies,
            predictors,
            epitope_masks_fname="%s/builds/flu/metadata/ha_masks.tsv" % augur_path,
            epitope_mask_version="wolf",
            tolerance_mask_version="HA1",
            min_freq=0.1,
            predictor_kwargs=predictor_kwargs,
            delta_time=json_dict["delta_time"],
            step_size=json_dict["step_size"]
        )
        model.prep_nodes()

        predictor_arrays = {}
        for key in json_dict["predictor_arrays"]:
            predictor_arrays[float(key)] = np.array(json_dict["predictor_arrays"][key])

        model.predictor_arrays = predictor_arrays

        freq_arrays = {}
        for key in json_dict["freq_arrays"]:
            freq_arrays[float(key)] = np.array(json_dict["freq_arrays"][key])

        model.freq_arrays = freq_arrays
        model.select_clades_for_fitting()

        return model

    def to_json(self, filename):
        """Export fitness model parameters, data, and accuracy statistics to JSON.
        """
        # Convert predictor parameters to a data frame to easily export as
        # records.
        params_df = pd.DataFrame({
            "predictor": self.predictors,
            "param": self.model_params.tolist(),
            "global_sd": self.global_sds.tolist()
        })

        correlation_null, correlation_raw, correlation_rel = self.get_correlation()
        mcc = get_matthews_correlation_coefficient_for_data_frame(self.pred_vs_true_df)

        # Do not try to export titer data if it was provided to the model.
        predictor_kwargs = self.predictor_kwargs.copy()
        if "transform" in predictor_kwargs:
            predictor_kwargs["transform"] = str(predictor_kwargs["transform"])

        if "titers" in predictor_kwargs:
            del predictor_kwargs["titers"]

        data = {
            "params": params_df.to_dict(orient="records"),
            "predictor_kwargs": predictor_kwargs,
            "data": self.pred_vs_true_df.to_dict(orient="records"),
            "accuracy": {
                "clade_error": self.clade_fit(self.model_params),
                "correlation_rel": correlation_rel[0],
                "mcc": mcc
            },
            "delta_time": self.delta_time,
            "step_size": self.timepoint_step_size,
            "end_date": self.end_date
        }

        predictor_arrays = {}
        for key in self.predictor_arrays:
            predictor_arrays[key] = self.predictor_arrays[key].tolist()

        data["predictor_arrays"] = predictor_arrays

        freq_arrays = {}
        for key in self.freq_arrays:
            freq_arrays[key] = self.freq_arrays[key].tolist()

        data["freq_arrays"] = freq_arrays

        write_json(data, filename)


def main(params):
    import time
    from io_util import read_json
    from io_util import write_json
    from tree_util import json_to_dendropy, dendropy_to_json

    print("--- Start fitness model optimization at " + time.strftime("%H:%M:%S") + " ---")

    tree_fname='data/tree_refine.json'
    tree =  json_to_dendropy(read_json(tree_fname))
    fm = fitness_model(tree, predictors = params['predictors'], verbose=1)
    fm.predict(niter = params['niter'])
    out_fname = "data/tree_fitness.json"
    write_json(dendropy_to_json(tree.root), out_fname)
    return out_fname

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Optimize predictor coefficients')
    parser.add_argument('-n', '--niter', type = int, default=10, help='number of replicate optimizations')
    parser.add_argument("-t", "--test", help="run test", action="store_true")
    parser.add_argument('-p', '--predictors', default=default_predictors, help='predictors to optimize', nargs='+')
    params = parser.parse_args().__dict__
    if params['test']:
        fm = test(params)
    else:
        main(params)
