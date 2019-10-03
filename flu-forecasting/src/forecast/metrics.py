"""Metrics for model comparison
"""
import numpy as np


def sum_of_squared_errors(observed, estimated, **kwargs):
    """
    Calculates the sum of squared errors for observed and estimated values.

    Parameters
    ----------
    observed : numpy.ndarray
        observed values

    estimated : numpy.ndarray
        estimated values

    Returns
    -------
    float :
        sum of squared errors between observed and estimated values
    """
    return np.sum((observed - estimated) ** 2)


def root_mean_square_error(observed, estimated, **kwargs):
    """
    Calculates the root mean square error between observed and estimated values.

    Parameters
    ----------
    observed : numpy.ndarray
        observed values

    estimated : numpy.ndarray
        estimated values

    Returns
    -------
    float :
        root mean square error between observed and estimated values
    """
    return np.sqrt(np.mean((observed - estimated) ** 2))


def mean_absolute_error(observed, estimated, **kwargs):
    """
    Calculates the mean absolute error between observed and estimated values.

    Parameters
    ----------
    observed : numpy.ndarray
        observed values

    estimated : numpy.ndarray
        estimated values

    Returns
    -------
    float :
        mean absolute error between observed and estimated values
    """
    return np.mean(np.abs(observed - estimated))


def negative_information_gain(observed, estimated, initial, **kwargs):
    """Calculates the negative information gain between a given model's estimated
    values and a null model estimate based on the initial values. Data are
    assumed to be based on time series where the initial values are collected at
    time t and the observed and estimated values are from time t + 1.

    The value returned is the negative of information gain to enable use of this
    function for minimization algorithms.

    Parameters
    ----------
    observed : numpy.ndarray
        observed values

    estimated : numpy.ndarray
        estimated values

    initial : numpy.ndarray
        initial values

    Returns
    -------
    float :
        sum of squared errors between observed and estimated values

    """
    # Consider only samples where the observed value is non-zero.
    nonzero_indices = np.nonzero(observed)[0]

    # Calculate the negative information gain across all samples.
    return -1 * (observed[nonzero_indices] * np.log(estimated[nonzero_indices] / initial[nonzero_indices])).sum()


def add_pseudocounts_to_frequencies(df, pseudocount, timepoint_column="timepoint"):
    """Add pseudocounts to the frequencies in the given data frame.

    Parameters
    ----------
    df : pandas.DataFrame
        table of samples including annotations for `timepoint`,
        `clade_membership`, and `frequency` unless an alternate timepoint column
        is provided

    pseudocount : float
        numerator of pseudocount to add to the frequency of each sample scaled
        by the number of other samples in the same clade as the given sample

    timepoint_column : string
        name of the timepoint column in the given data frame

    Returns
    -------
    pandas.DataFrame :
        copy of the input data frame with a modified `frequency` column

    """
    assert all([
        column in df.columns
        for column in (timepoint_column, "clade_membership", "frequency")
    ])

    # Calculate clade sample counts per timepoint and clade.
    pseudocount_by_timepoint_and_clade = df.groupby([
        timepoint_column,
        "clade_membership"
    ])["frequency"].count().reset_index().rename(columns={"frequency": "samples"})

    # Calculate the scaled pseudocount value per timepoint and clade.
    pseudocount_by_timepoint_and_clade["total_pseudocount"] = pseudocount
    pseudocount_by_timepoint_and_clade["pseudocount"] = (
        pseudocount_by_timepoint_and_clade["total_pseudocount"] / pseudocount_by_timepoint_and_clade["samples"]
    )

    # Calculate the total frequency per timepoint after adding pseudocounts.
    total_frequency_by_timepoint = (
        pseudocount_by_timepoint_and_clade.groupby(timepoint_column)["total_pseudocount"].sum() + df.groupby(timepoint_column)["frequency"].sum()
    ).reset_index().rename(columns={0: "total_frequency"})

    # Annotate pseudocounts with total frequency per timepoint.
    pseudocount_by_timepoint_and_clade = pseudocount_by_timepoint_and_clade.merge(
        total_frequency_by_timepoint,
        on=[timepoint_column]
    )

    # Annotate original data frame with pseudocounts.
    df_with_pseudocounts = df.merge(
        pseudocount_by_timepoint_and_clade,
        on=[timepoint_column, "clade_membership"]
    )

    # Calculate new normalized frequencies with pseudocounts for all timepoints
    # that have non-zero total frequencies.
    timepoint_frequencies = df.groupby(timepoint_column)["frequency"].sum().reset_index()
    nonzero_timepoints = timepoint_frequencies[timepoint_frequencies["frequency"] > 0][timepoint_column].unique()
    nonzero_indices = df_with_pseudocounts[timepoint_column].isin(nonzero_timepoints)

    df_with_pseudocounts.loc[nonzero_indices, "frequency"] = (
        (df_with_pseudocounts.loc[nonzero_indices, "frequency"] + df_with_pseudocounts.loc[nonzero_indices, "pseudocount"]) / df_with_pseudocounts.loc[nonzero_indices, "total_frequency"]
    )

    df = df_with_pseudocounts.drop(columns=[
        "samples",
        "total_pseudocount",
        "pseudocount",
        "total_frequency"
    ])

    return df
