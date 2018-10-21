import argparse
import numpy as np
from collections import defaultdict
from augur.utils import read_metadata, get_numerical_dates

vpm_dict = {
    2: 92,
    3: 62,
    6: 32,
    12: 18,
}

regions = [
    ('africa',            "",   1.02),
    ('europe',            "EU", 0.74),
    ('north_america',     "NA", 0.54),
    ('china',             "AS", 1.36),
    ('south_asia',        "AS", 1.45),
    ('japan_korea',       "AS", 0.20),
    ('oceania',           "OC", 0.04),
    ('south_america',     "SA", 0.41),
    ('southeast_asia',    "AS", 0.62),
    ('west_asia',         "AS", 0.75)
]

def populate_counts(seqs_to_count):
    sequence_count_total = defaultdict(int)
    sequence_count_region = defaultdict(int)
    for seq in seqs_to_count:
        sequence_count_total[(seq.attributes['date'].year,
                              seq.attributes['date'].month)]+=1
        sequence_count_region[(seq.attributes['region'],
                              seq.attributes['date'].year,
                              seq.attributes['date'].month)]+=1
    return (sequence_count_total, sequence_count_region)


def flu_subsampling(strains, viruses_per_month, titer_values=None):
    category = lambda x: (x['region'],
                          x['year'],
                          x['month'])

    #### DEFINE THE PRIORITY
    if titer_values is not None:
        HI_titer_count = TiterCollection.count_strains(titer_values)
    else:
        print("Couldn't load titer information - using random priorities")
        HI_titer_count = False
        def priority(seq):
            return np.random.random() + int(seq.name in reference_viruses[params.lineage])

    if HI_titer_count:
        def priority(seq):
            sname = seq.attributes['strain']
            if sname in HI_titer_count:
                pr = HI_titer_count[sname]
            else:
                pr = 0
            return (pr + len(seq.seq)*0.0001 - 0.01*np.sum([seq.seq.count(nuc) for nuc in 'NRWYMKSHBVD']) +
                    1e6*int(seq.name in reference_viruses[params.lineage]))


    region_threshold = int(np.ceil(1.0*viruses_per_month/len(regions)))
    def threshold(obj):
        """
        a higher order function which returns a fn which has access to
        some summary stats about the sequences (closure)
        """
        sequence_count_total, sequence_count_region = populate_counts(obj)
        def threshold_fn(x):
            #x is the collection key, in this case a tuple of (region, year, month)
            if sequence_count_total[(x[1], x[2])] < viruses_per_month:
                return viruses_per_month
            region_counts = sorted([sequence_count_region[(r[0], x[1], x[2])] for r in regions])
            if region_counts[0] > region_threshold:
                return region_threshold
            left_to_fill = viruses_per_month - len(regions)*region_counts[0]
            thres = region_counts[0]
            for ri, rc in zip(range(len(regions)-1, 0, -1), region_counts[1:]):
                if left_to_fill - ri*(rc-thres)>0:
                    left_to_fill-=ri*(rc-thres)
                    thres = rc
                else:
                    thres += left_to_fill/ri
                    break
            return max(1, int(thres))
        return threshold_fn

    return {
        "category": category,
        "priority": priority,
        "threshold": threshold
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Prepare fauna FASTA for analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('-v', '--viruses_per_month', type = int, default=15,
                        help='Subsample x viruses per country per month. Set to 0 to disable subsampling.')
    parser.add_argument('--metadata', nargs='+', help="FASTA file of virus sequences from fauna (e.g., zika.fasta)")
    parser.add_argument('--output', help="name of the file with selected strains")
    parser.add_argument('--verbose', action="store_true", help="turn on verbose reporting")

    parser.add_argument('-l', '--lineage', choices=['h3n2', 'h1n1pdm', 'vic', 'yam'], default='h3n2', type=str, help="single lineage to include (default: h3n2)")
    parser.add_argument('-r', '--resolution', choices=['2y', '3y', '6y', '12y'], default=['3y'], type = str,  help = "single resolution to include (default: 3y)")
    parser.add_argument('--ensure_all_segments', action="store_true", default=False,  help = "exclude all strains that don't have the full set of segments")
    parser.add_argument('-s', '--segments', default=['ha'], nargs='+', type = str,  help = "list of segments to include (default: ha)")
    parser.add_argument('--sampling', default = 'even', type=str,
                        help='sample evenly over regions (even) (default), or prioritize one region (region name), otherwise sample randomly')
    parser.add_argument('--time_interval', nargs=2, help="explicit time interval to use -- overrides resolutions"
                                                                     " expects YYYY-MM-DD YYYY-MM-DD")
    parser.add_argument('--strains', help="a text file containing a list of strains (one per line) to prepare without filtering or subsampling")
    parser.add_argument('--titers', help="tab-delimited file of titer strains and values from fauna (e.g., h3n2_hi_titers.tsv)")
    parser.add_argument('--complete_frequencies', action='store_true', help="compute mutation frequences from entire dataset")

    args = parser.parse_args()

    metadata = {}
    for segment, fname in zip(args.segments, args.metadata):
        tmp_meta, columns = read_metadata(fname)
        numerical_dates = get_numerical_dates(tmp_meta, fmt='%Y-%m-%d')
        for x in tmp_meta:
            tmp_meta[x]['num_date'] = np.mean(numerical_dates[x])
            tmp_meta[x]['year'] = int(tmp_meta[x]['num_date'])
            tmp_meta[x]['month'] = int((tmp_meta[x]['num_date']%1)*12)
        metadata[segment] = tmp_meta

    guide_segment = args.segments[0]
    strains_with_all_segments = set.intersection(*(set(metadata[x].keys()) for x in args.segments))
    subsampling = flu_subsampling({x:metadata[guide_segment][x] for x in strains_with_all_segments},
                                  args.viruses_per_month)

    virus_by_category = defaultdict(list)
    for v in strains_with_all_segments:
        virus_by_category[subsampling["category"](metadata[guide_segment][v])].append(v)

    selected_strains = []
    for cat, val in virus_by_category.items():
        selected_strains.extend(val[:args.viruses_per_month])

    with open(args.output, 'w') as ofile:
        ofile.write('\n'.join(selected_strains))
