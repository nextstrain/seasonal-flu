import json, argparse



pathogen_data = {
  "schemaVersion": "3.0.0",
  "alignmentParams": {
    "excessBandwidth": 9,
    "terminalBandwidth": 100,
    "allowedMismatches": 4,
    "gapAlignmentSide": "right",
    "minSeedCover": 0.1
  },
  "compatibility": {
    "cli": "3.0.0-alpha.0",
    "web": "3.0.0-alpha.0"
  },
  "defaultCds": "HA1",
  "experimental": True,
  "qc": {
    "privateMutations": {
      "enabled": True,
      "typical": 5,
      "cutoff": 15,
      "weightLabeledSubstitutions": 2,
      "weightReversionSubstitutions": 1,
      "weightUnlabeledSubstitutions": 1
    },
    "missingData": {
      "enabled": False,
      "missingDataThreshold": 100,
      "scoreBias": 10
    },
    "snpClusters": {
      "enabled": False,
      "windowSize": 100,
      "clusterCutOff": 5,
      "scoreWeight": 50
    },
    "mixedSites": {
      "enabled": True,
      "mixedSitesThreshold": 4
    },
    "frameShifts": {
      "enabled": True
    },
    "stopCodons": {
      "enabled": True,
      "ignoredStopCodons": []
    }
  },
  "cdsOrderPreference": [
    "HA1",
    "HA2"
  ],
  "maintenance": {
    "website": [
      "https://nextstrain.org",
      "https://clades.nextstrain.org"
    ],
    "documentation": [
      "https://github.com/nextstrain/seasonal-flu"
    ],
    "source code": [
      "https://github.com/nextstrain/seasonal_flu"
    ],
    "issues": [
      "https://github.com/nextstrain/seasonal_flu/issues"
    ],
    "organizations": [
      "Nextstrain"
    ],
    "authors": [
      "Nextstrain team <https://nextstrain.org>"
    ]
  }
}

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="fix genome clade info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input-auspice-json', type=str, required=True, help="input auspice_json")
    parser.add_argument('--reference', type=str, required=True, help="reference")
    parser.add_argument('--build-name', type=str, required=True, help="nextclade build name")
    parser.add_argument('--reference-accession', type=str, required=True, help="reference accession")
    parser.add_argument('--output', type=str, metavar="JSON", required=True, help="output Auspice JSON")
    args = parser.parse_args()

    with open(args.input_auspice_json, 'r') as fh:
        data = json.load(fh)

    pathogen_data["attributes"] = {"reference accession": args.reference_accession, "reference name": args.reference, "name": args.build_name}

    data["meta"]["extensions"] = {'nextclade': {'pathogen': pathogen_data}}

    with open(args.output, 'w') as fh:
        json.dump(data, fh, indent=0)

