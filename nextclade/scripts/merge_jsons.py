import json, argparse

def get_clade_configs(name):
    return {
    "short_clades": {
        "name": "short_clade",
        "displayName": "Abbreviated clade name",
        "description": "For recent subclades with long names, the prefix describing their history is omitted."
    },
    "subclade": {
        "name": "subclade",
        "displayName": "Subclade",
        "description": "Experimental fine-grained subclade annotation."
    }}.get(name, {'name':name, "displayName":name})


if __name__=="__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--lineage", required=True, help="attribute info")
    parser.add_argument("--segment", required=True, help="attribute info")
    parser.add_argument("--reference", required=True, help="attribute info")
    parser.add_argument("--auspice-config", required=True, help="Auspice config JSON with coloring entry to have scale added to")
    parser.add_argument("--pathogen-jsons", nargs='+', required=True, help="name of the coloring field in the Auspice config JSON")
    parser.add_argument("--clades", nargs="+", required=True, help="list of values to assign colors to")
    parser.add_argument("--output-auspice", required=True, help="Auspice config JSON with scale added to the requested coloring")
    parser.add_argument("--output-pathogen", required=True, help="Auspice config JSON with scale added to the requested coloring")
    args = parser.parse_args()

    pathogen_json = {}
    for p in args.pathogen_jsons:
        with open(p) as fh:
            pathogen_json.update(json.load(fh))

    with open(args.auspice_config) as fh:
        auspice_json = json.load(fh)

    pathogen_json['attributes'] = {"name":{"value":args.lineage},
                                   "segment":{"value":args.segment},
                                   "reference":{"value":args.reference}}


    if len(args.clades):
        auspice_json['extensions']['nextclade']["clade_node_attrs"] =  [
            get_clade_configs(c) for c in args.clades
        ]

    with open(args.output_pathogen, 'w') as fh:
        json.dump(pathogen_json, fh, indent=2)

    with open(args.output_auspice, 'w') as fh:
        json.dump(auspice_json, fh, indent=2)

