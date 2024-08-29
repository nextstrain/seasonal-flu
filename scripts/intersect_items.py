#!/usr/bin/env python3
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--items", nargs="+", required=True, help="one or more files containing a list of items")
    parser.add_argument("--output", required=True, help="list of items shared by all input files (the intersection)")

    args = parser.parse_args()

    with open(args.items[0], "r", encoding="utf-8") as fh:
        shared_items = {line.strip() for line in fh}

    for item_file in args.items[1:]:
        with open(item_file, "r", encoding="utf-8") as fh:
            items = {line.strip() for line in fh}

        shared_items = shared_items & items

    with open(args.output, "w", encoding="utf-8") as oh:
        for item in sorted(shared_items):
            print(item, file=oh)
