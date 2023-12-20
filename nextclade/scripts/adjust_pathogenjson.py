import json, argparse

gene_to_segment = {
    'pb2': 1,
    'pb1': 2,
    'pa': 3,
    'ha': 4,
    'np': 5,
    'na': 6,
    'mp': 7,
    'ns': 8,
}

if __name__=="__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--lineage", required=True, help="attribute info")
    parser.add_argument("--segment", required=True, help="attribute info")
    parser.add_argument("--reference", required=True, help="attribute info")
    parser.add_argument("--reference-name", required=True, help="attribute info")
    parser.add_argument("--default-cds", required=True, help="attribute info")
    parser.add_argument("--pathogen-json",  required=True)
    parser.add_argument("--output-pathogen", required=True, help="Auspice config JSON with scale added to the requested coloring")
    args = parser.parse_args()

    pathogen_json = {}
    with open(args.pathogen_json) as fh:
        pathogen_json.update(json.load(fh))

    flu_type = {'h3n2':'A', 'h1n1pdm':'A', 'vic':'B', 'yam':'B'}[args.lineage]
    lineage_name = {'h3n2':'H3N2', 'h1n1pdm':'H1N1pdm', 'vic':'Victoria', 'yam':'Yamagata'}[args.lineage]

    pathogen_json['attributes'] = {"name": f"Influenza {flu_type} {lineage_name} {args.segment.upper()} (segment {gene_to_segment[args.segment]})",
                                   "segment": args.segment,
                                   "reference accession": args.reference,
                                   "reference name": args.reference_name}

    pathogen_json['cdsOrderPreference'] = {"ha": ["HA1", "HA2"], "na":["NA"]}.get(args.segment, [])
    if args.default_cds:
        pathogen_json['defaultCds'] = args.default_cds

    with open(args.output_pathogen, 'w') as fh:
        json.dump(pathogen_json, fh, indent=2)

