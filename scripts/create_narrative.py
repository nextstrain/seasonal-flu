#!/usr/bin/env python3
import argparse
import datetime
import sys

from jinja2 import Environment, FileSystemLoader, select_autoescape
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--template", required=True, help="Template of Markdown narrative in Jinja2 format.")
    parser.add_argument("--earliest-date", help="Earliest date in format of YYYY-MM-DD to use for year/month filters in Nextstrain views. Defaults to 6 months before the given report date.")
    parser.add_argument("--date", help="Date in the format of YYYY-MM-DD to use for the report. Defaults to the current date.")
    parser.add_argument("--markdown-includes", nargs="+", help="One or more key/value pairs of variable name and path to a Markdown file to be included in the template environment in the format of 'variable_name=path_to_markdown_file.md'.")
    parser.add_argument("--output", required=True, help="Rendered Markdown narrative.")

    args = parser.parse_args()

    env = Environment(
        loader=FileSystemLoader("."),
        autoescape=select_autoescape()
    )

    template = env.get_template(args.template)

    if args.date:
        date = datetime.datetime.strptime(args.date, "%Y-%m-%d")
        date_string = args.date
    else:
        date = datetime.datetime.today()
        date_string = date.strftime("%Y-%m-%d")

    # Build year/month filters for Nextstrain views based on the current report
    # date and the earliest date requested. If no earliest date is given,
    # default to a fixed number of months in the past. This dynamic construction
    # of the date filters ensures that we always include the most recent month
    # in the filters.
    if args.earliest_date:
        earliest_date_string = args.earliest_date
    else:
        earliest_date = date - pd.DateOffset(months=6)
        earliest_date_string = earliest_date.strftime("%Y-%m-01")

    year_month_periods = pd.period_range(
        earliest_date_string,
        date_string,
        freq="M",
    )
    year_month_filter = ",".join(
        str(period)
        for period in year_month_periods
    )

    # Setup variables that the template always expects to find.
    variables = {
        "date": date_string,
        "long_month": date.strftime("%B"),
        "day": date.strftime("%d"),
        "year": date.strftime("%Y"),
        "year_month_filter": year_month_filter,
    }

    # Include the contents of additional Markdown files that have been generated
    # by the workflow.
    markdown_includes = args.markdown_includes
    if markdown_includes is None:
        markdown_includes = []

    for markdown_include in markdown_includes:
        variable_name, markdown_path = markdown_include.split("=")
        if variable_name in variables:
            print(
                f"Warning: Overwriting a variable of name '{variable_name}' that already exists in the template's environment.",
                file=sys.stderr
            )

        with open(markdown_path, "r", encoding="utf-8") as fh:
            variables[variable_name] = fh.read()

    with open(args.output, "w", encoding="utf-8") as oh:
        print(template.render(**variables), file=oh)
