import argparse
from augur.utils import read_node_data, write_json


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--clades", required=True, help="clades to rename attribute for")
    parser.add_argument("--new-attribute", required=True, help="name of the new clade attribute for full clades")
    parser.add_argument("--output", required=True, help="clades with renamed attributes")
    args = parser.parse_args()

    data = read_node_data(args.clades)
    for strain in data["nodes"].keys():
        clade = data["nodes"][strain].pop("clade_membership")
        data["nodes"][strain] = {
            args.new_attribute: clade,
        }

    write_json(data, args.output)
