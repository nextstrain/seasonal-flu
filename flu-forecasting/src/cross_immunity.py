"""Calculate cross-immunity for individual samples against past samples based on given pairwise distances and scaled by past frequencies.
"""
import argparse
from augur.utils import read_node_data
import json

from forecast.fitness_predictors import inverse_cross_immunity_amplitude, cross_immunity_cost


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Calculate cross-immunity",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--frequencies", required=True, help="JSON of frequencies per sample")
    parser.add_argument("--distances", required=True, help="JSON of distances between samples")
    parser.add_argument("--date-annotations", required=True, help="JSON of branch lengths and date annotations from augur refine for samples in the given tree")
    parser.add_argument("--distance-attributes", nargs="+", required=True, help="names of attributes to use from the given distances JSON")
    parser.add_argument("--immunity-attributes", nargs="+", required=True, help="names of attributes to use for the calculated cross-immunities")
    parser.add_argument("--decay-factors", nargs="+", required=True, type=float, help="list of decay factors (d_0) for each given immunity attribute")
    parser.add_argument("--years-to-wane", type=int, help="number of years after which immunity wanes completely")
    parser.add_argument("--output", required=True, help="cross-immunities calculated from the given distances and frequencies")
    args = parser.parse_args()

    # Load frequencies.
    with open(args.frequencies, "r") as fh:
        frequencies = json.load(fh)

    # Identify maximum frequency per sample.
    max_frequency_per_sample = {
        sample: float(max(sample_frequencies["frequencies"]))
        for sample, sample_frequencies in frequencies.items()
        if sample not in ["pivots", "generated_by"] and not sample.startswith("count")
    }
    current_timepoint = frequencies["pivots"][-1]

    # Load distances.
    with open(args.distances, "r") as fh:
        distances = json.load(fh)

    distances = distances["nodes"]

    # Load date annotations and annotate tree with them.
    date_annotations = read_node_data(args.date_annotations)
    date_by_node_name = {}
    for node, annotations in date_annotations["nodes"].items():
        date_by_node_name[node] = annotations["numdate"]

    """
  "A/Acre/15093/2010": {
   "ep": 9,
   "ne": 8,
   "rb": 3
  },
    """
    if args.years_to_wane is not None:
        print("Waning effect with max years of %i" % args.years_to_wane)
    else:
        print("No waning effect")

    # Calculate cross-immunity for distances defined by the given attributes.
    cross_immunities = {}
    for sample, sample_distances in distances.items():
        for distance_attribute, immunity_attribute, decay_factor in zip(args.distance_attributes, args.immunity_attributes, args.decay_factors):
            if distance_attribute not in sample_distances:
                continue

            if sample not in cross_immunities:
                cross_immunities[sample] = {}

            # Calculate cross-immunity cost from all distances to the current
            # sample. This negative value increases for samples that are
            # increasingly distant from previous samples.
            cross_immunity = 0.0
            for past_sample, distance in sample_distances[distance_attribute].items():
                # Calculate effect of waning immunity.
                if args.years_to_wane is not None:
                    waning_effect = max(1 - ((current_timepoint - date_by_node_name[past_sample]) / args.years_to_wane), 0)
                else:
                    waning_effect = 1.0

                # Calculate cost of cross-immunity with waning.
                if waning_effect > 0:
                    cross_immunity += waning_effect * max_frequency_per_sample[past_sample] * cross_immunity_cost(
                        distance,
                        decay_factor
                    )

            cross_immunities[sample][immunity_attribute] = -1 * cross_immunity

    # Export cross-immunities to JSON.
    with open(args.output, "w") as oh:
        json.dump({"nodes": cross_immunities}, oh, indent=1, sort_keys=True)
