import json
import argparse

def get_clade_configs(name):
    return {
    "legacy-clade": {
        "name": "legacy-clade",
        "displayName": "Legacy clade",
        "description": "Previous clades system that is no longer actively used.",
        "hideInWeb": True,
    },
    "short-clade": {
        "name": "short-clade",
        "displayName": "Abbr.legacy clade name",
        "description": "Shortened names of legacy clades.",
        "skipAsReference": True
    },
    "subclade": {
        "name": "subclade",
        "displayName": "Subclade",
        "description": "Fine-grained subclade annotation.",
        "hideInWeb": True,
        "skipAsReference": True
    },
    "proposedSubclade": {
        "name": "proposedSubclade",
        "displayName": "Proposed subclade",
        "description": "Includes proposals of new subclades. These can change anytime.",
        "hideInWeb": True,
        "skipAsReference": True
    }}.get(name, {'name':name, "displayName":name, "description":""})

default_CDS = {
    "ha": "HA1",
    "na": "NA",
    "pb2": "PB2",
    "pb1": "PB1",
    "pa": "PA",
    "np": "NP",
    "mp": "M1",
    "ns": "NS1"
}



if __name__=="__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--lineage", required=True, help="attribute info")
    parser.add_argument("--segment", required=True, help="attribute info")
    parser.add_argument("--reference", required=True, help="attribute info")
    parser.add_argument("--reference-name", required=True, help="attribute info")
    parser.add_argument("--auspice-config", required=True, help="Auspice config JSON with coloring entry to have scale added to")
    parser.add_argument("--pathogen-jsons", nargs='+', required=True, help="name of the coloring field in the Auspice config JSON")
    parser.add_argument("--clades", nargs="+", required=False, help="list of values to assign colors to")
    parser.add_argument("--output-auspice", required=True, help="Auspice config JSON with scale added to the requested coloring")
    parser.add_argument("--output-pathogen", required=True, help="Auspice config JSON with scale added to the requested coloring")
    args = parser.parse_args()

    pathogen_json = {}
    for p in args.pathogen_jsons:
        with open(p) as fh:
            pathogen_json.update(json.load(fh))

    with open(args.auspice_config) as fh:
        auspice_json = json.load(fh)

    flu_type = {'h3n2':'A', 'h1n1pdm':'A', 'vic':'B', 'yam':'B'}[args.lineage]
    lineage_name = {'h3n2':'H3N2', 'h1n1pdm':'H1N1pdm', 'vic':'Victoria', 'yam':'Yamagata'}[args.lineage]

    pathogen_json['attributes'] = {"name": f"Influenza {flu_type} {lineage_name} {args.segment.upper()}",
                                   "segment": args.segment,
                                   "reference accession": args.reference,
                                   "reference name": args.reference_name}

    pathogen_json['cdsOrderPreference'] = {"ha": ["HA1", "HA2"], "na":["NA"]}.get(args.segment, [])
    pathogen_json['defaultCds'] = default_CDS[args.segment]

    if args.clades and len(args.clades)>0:
        auspice_json['extensions']['nextclade']["clade_node_attrs"] =  [
            get_clade_configs(c) for c in args.clades if c not in ['default']
        ]

    if args.segment in ['ha', 'na']:
        auspice_json['display_defaults'] = {
            "color_by": "clade_membership",
            "branch_label": "clade"
        }
        auspice_json['filters'].extend(['clade_membership', 'subclade'])
    else:
        auspice_json['display_defaults'] = {
            "color_by": "QC_status",
        }

    with open(args.output_pathogen, 'w') as fh:
        json.dump(pathogen_json, fh, indent=2)

    with open(args.output_auspice, 'w') as fh:
        json.dump(auspice_json, fh, indent=2)
