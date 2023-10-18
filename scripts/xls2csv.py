"""Minimal script to convert Excel XLS format to CSV.
"""
import argparse
import csv
import xlrd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--xls", required=True, help="path to XLS file to convert")
    parser.add_argument("--output", required=True, help="path to CSV output")

    args = parser.parse_args()

    workbook = xlrd.open_workbook(args.xls)
    sheet = workbook.sheet_by_index(0)
    field_names = [field.value for field in sheet.row(0)]

    with open(args.output, "w", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile, dialect="unix")
        writer.writerow(field_names)

        for row_index in range(1, sheet.nrows):
            writer.writerow([field.value for field in sheet.row(row_index)])
