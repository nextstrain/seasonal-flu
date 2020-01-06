import argparse, json
import numpy as np
from flu_regions import region_properties

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

    all_mutations = sorted(filter(lambda x:('counts' not in x) and ('pivots' not in x) and (x != 'generated_by'),
                           set.union(*[set(frequencies[region].keys()) for region in frequencies])))
    pivots = frequencies[args.regions[0]]['pivots']
    frequencies['global']['pivots'] = format_frequencies(pivots)

    # determine the seasonal profiles and weights of different regions for later equitable global frequency calculation below
    seasonal_profile = {}
    total_weights = {}
    for region in frequencies:
        props = region_properties[region]
        all_gene_counts = sorted(filter(lambda x:'counts' in x, frequencies[region].keys()))
        seasonal_profile[region] = {}
        for x in all_gene_counts:
            gene = x.split(':')[0]
            # calculate a smoothed seasonal profile from the samples counts.
            tmp_counts = np.convolve(np.ones(6)/6.0, np.array(frequencies[region][x]), mode='same')
            # calculate temporal weight for each region as the product of region population size and seasonal pattern
            seasonal_profile[region][gene] = props['popsize']*np.ones_like(tmp_counts) #np.maximum(0.1, tmp_counts/tmp_counts.max())
            if gene not in total_weights: total_weights[gene]=[]
            total_weights[gene].append(seasonal_profile[region][gene])

    # sum the weights for each gene to obtain the total weight (denominator)
    for gene in total_weights:
        total_weights[gene] = np.sum(total_weights[gene], axis=0)

    # calculate the mutation frequencies by spatio-temporal weighing
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
        props = region_properties[region]
        for mutation in frequencies[region]:
            key = '%s_%s'%(props.get('abbr', region), mutation)
            json_for_export[key] = frequencies[region][mutation]

    if args.tree_frequencies:
        with open(args.tree_frequencies) as fh:
            json_for_export.update(json.load(fh))

    with open(args.output_auspice, 'wt') as fh:
        json.dump(json_for_export, fh, indent=1)
