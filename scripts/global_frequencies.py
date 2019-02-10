import argparse, json
import numpy as np

population_sizes = {
    'africa':1.02,
    'europe': 0.74,
    'north_america': 0.54,
    'china': 1.36,
    'south_asia': 1.45,
    'japan_korea': 0.20,
    'oceania': 0.04,
    'south_america': 0.41,
    'southeast_asia': 0.62,
    'west_asia': 0.75
}

region_abbreviations = {
    'africa':'AF',
    'europe': 'EU',
    'north_america': 'NA',
    'china': 'CN',
    'south_asia': 'SAS',
    'japan_korea': 'JK',
    'oceania': 'OC',
    'south_america': 'SA',
    'southeast_asia': 'SEA',
    'west_asia': "WAS"
}

def format_frequencies(x):
    return [round(y,4) for y in x]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Combine regional frequencies into a global frequency estimate and export as json",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--region-frequencies', nargs='+', type=str, help="regions with frequency estimates")
    parser.add_argument('--regions', nargs='+', type=str, help="region names corresponding to estimated frequencies")
    parser.add_argument('--tree-frequencies', type=str, help="json file with tree frequencies")
    parser.add_argument('--output-augur', type=str,  help="name of file to frequencies too ")
    parser.add_argument('--output-auspice', type=str,  help="name of file to frequencies too ")

    args = parser.parse_args()

    frequencies = {'global':{}}
    for region, freq_file in zip(args.regions, args.region_frequencies):
        with open(freq_file) as fh:
            frequencies[region] = json.load(fh)

    all_mutations = sorted(filter(lambda x:('counts' not in x) and ('pivots' not in x),
                           set.union(*[set(frequencies[region].keys()) for region in frequencies])))
    pivots = frequencies[args.regions[0]]['pivots']
    frequencies['global']['pivots'] = format_frequencies(pivots)

    seasonal_profile = {}
    total_weights = {}
    for region in frequencies:
        all_genes = sorted(filter(lambda x:'counts' in x, frequencies[region].keys()))
        seasonal_profile[region] = {}
        for x in all_genes:
            gene = x.split(':')[0]
            tmp_counts = np.convolve(np.ones(6)/6.0, np.array(frequencies[region][x]), mode='same')
            seasonal_profile[region][gene] = np.maximum(0.1, tmp_counts/tmp_counts.max())
            if gene not in total_weights: total_weights[gene]=[]
            total_weights[gene].append(seasonal_profile[region][gene])

    for gene in total_weights:
        total_weights[gene] = np.sum(total_weights[gene], axis=0)

    for mutation in all_mutations:
        gene = mutation.split(':')[0]
        freqs = []
        weights = []
        for region in frequencies:
            if mutation in frequencies[region]:
                freqs.append(frequencies[region][mutation])
                weights.append(seasonal_profile[region][gene])

        frequencies['global'][mutation] = format_frequencies(
            np.sum(np.array(freqs)*np.array(weights), axis=0) / total_weights[gene]
        )

    with open(args.output_augur, 'wt') as fh:
        json.dump(frequencies, fh, indent=1)

    json_for_export = {'pivots':format_frequencies(pivots)}
    for region in frequencies:
        for mutation in frequencies[region]:
            key = '%s_%s'%(region_abbreviations.get(region, region), mutation)
            json_for_export[key] = frequencies[region][mutation]

    if args.tree_frequencies:
        with open(args.tree_frequencies) as fh:
            json_for_export.update(json.load(fh))

    with open(args.output_auspice, 'wt') as fh:
        json.dump(json_for_export, fh, indent=1)
