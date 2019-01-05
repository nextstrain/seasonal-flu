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
    parser.add_argument('--output', type=str,  help="names of file to save age_distribution histogram to ")

    args = parser.parse_args()

    frequencies = {'global':{}}
    for region, freq_file in zip(args.regions, args.region_frequencies):
        with open(freq_file) as fh:
            frequencies[region] = json.load(fh)

    all_mutations = sorted(filter(lambda x:('counts' not in x) and ('pivots' not in x),
                           set.union(*[set(frequencies[region].keys()) for region in frequencies])))
    pivots = frequencies[args.regions[0]]['pivots']

    seasonal_profile = {}
    for region in frequencies:
        all_counts = sorted(filter(lambda x:'counts' in x, frequencies[region].keys()))
        seasonal_profile[region] = {}
        for x in all_counts:
            gene = x.split(':')[0]
            tmp = np.array(frequencies[region][x])
            seasonal_profile[region][gene] = (tmp+0.05*tmp.max())/(tmp.mean()+0.05*tmp.max())


    for mutation in all_mutations:
        gene = mutation.split(':')[0]
        freqs = []
        weights = []
        for region in frequencies:
            if mutation in frequencies[region]:
                freqs.append(frequencies[region][mutation])
                weights.append(seasonal_profile[region][gene]*population_sizes[region])

        frequencies['global'][mutation] = format_frequencies(np.sum(np.array(freqs)*np.array(weights), axis=0)/np.sum(weights, axis=0))

    json_for_export = {}
    for region in frequencies:
        for mutation in frequencies[region]:
            key = '%s_%s'%(region_abbreviations.get(region, region), mutation)
            json_for_export[key] = frequencies[region][mutation]

    with open(args.output, 'wt') as fh:
        json.dump(json_for_export, fh)
