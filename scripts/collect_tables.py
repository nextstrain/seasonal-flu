"""Collect two or more data frames in tab-delimited format into a single table.

Optionally, add new columns based on transformations of existing
columns. Transformations are limited to functions defined in the numpy
namespace.
"""
import argparse
import numpy as np
import pandas as pd
import re
import sys


class BadTransformException(Exception):
    pass


def parse_transform(transform):
    """Parse a transform string into its corresponding new column name, Python function, and original column name.

    Return `None` for the Python function if the requested function is not valid.

    Parameters
    ----------
    transform : str
        transformation definition string (e.g., "log_lbi=log(lbi)")

    Returns
    -------
    str, callable, str
        new column name, transformation function, and original column name

    >>> parse_transform("log_lbi=log(lbi)")
    ('log_lbi', <ufunc 'log'>, 'lbi')
    >>> parse_transform("fake_col=fake(col)")
    Traceback (most recent call last):
        ...
    collect_tables.BadTransformException: the requested function was invalid

    >>> parse_transform("bad_transform")
    Traceback (most recent call last):
        ...
    collect_tables.BadTransformException: the requested transform was malformed

    """
    match = re.match(r"(?P<new_column>\w+)=(?P<function>\w+)\((?P<column>\w+)\)", transform)
    if match is None:
        raise BadTransformException("the requested transform was malformed")

    new_column, function_string, column = match.groups()
    function = getattr(np, function_string, None)
    if function is None:
        raise BadTransformException("the requested function was invalid")

    return new_column, function, column


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Collect two or more data frame tables",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tables", nargs="+", required=True, help="tab-delimited files with the same columns to be collected into a single file")
    parser.add_argument("--transforms", nargs="+", help="a list of new columns to create by transformation of existing columns (e.g., 'log_lbi=log(lbi)')")
    parser.add_argument("--output", required=True, help="tab-delimited output file collecting the given input tables")
    args = parser.parse_args()

    # Concatenate tip attributes across all timepoints.
    df = pd.concat([pd.read_table(table) for table in args.tables], ignore_index=True)

    # Apply transformations.
    if args.transforms:
        for transform in args.transforms:
            try:
                new_column, transform_function, column = parse_transform(transform)
                df[new_column] = transform_function(df[column])
            except BadTransformException as e:
                print(f"Error: Could not apply transformation '{transform}' because {e}", file=sys.stderr)
            except Exception as e:
                print(f"Error: Failed to apply transformation '{transform}' ({e})", file=sys.stderr)

    df.to_csv(args.output, sep="\t", index=False)
