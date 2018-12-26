import argparse
import sys
import os
from collections import defaultdict
from datetime import datetime, timedelta, date
import Bio
import Bio.SeqIO
import numpy as np
from treetime.utils import numeric_date
from augur.utils import read_metadata, get_numerical_dates

regions = ['africa', 'europe', 'north_america', 'china', 'south_asia', 'japan_korea', 'oceania', 'south_america', 'southeast_asia', 'west_asia']
subcats = regions

def read_strain_list(fname):
    """
    read strain names from a file assuming there is one strain name per line

    Parameters:
    -----------
    fname : str
        file name

    Returns:
    --------
    strain_list : list
        strain names

    """
    if os.path.isfile(fname):
        with open(fname, 'r') as fh:
            strain_list = [x.strip() for x in fh.readlines() if x[0]!='#']
    else:
        print("ERROR: file %s containing strain list not found"%fname)
        sys.exit(1)

    return strain_list


def count_titer_measurements(fname):
    """
    read how many titer measurements exist for each virus

    Parameters:
    -----------
    fname : str
        file name

    Returns:
    --------
    titer_count : defaultdict(int)
        dictionary with titer count for each strain
    """
    titer_count = defaultdict(int)
    if os.path.isfile(fname):
        with open(fname, 'r') as fh:
            for line in fh:
                titer_count[line.split()[0]] += 1
    else:
        print("ERROR: file %s containing strain list not found"%fname)
        sys.exit(1)

    return titer_count


def populate_categories(metadata):
    super_category = lambda x: (x['region'],
                                x['year'],
                                x['month'])

    category = lambda x: (x['region'],
                          x['year'],
                          x['month'])

    virus_by_category = defaultdict(list)
    virus_by_super_category = defaultdict(list)
    for v in metadata:
        virus_by_category[category(metadata[v])].append(v)
        virus_by_super_category[super_category(metadata[v])].append(v)

    return virus_by_super_category, virus_by_category


def flu_subsampling(metadata, viruses_per_month, time_interval, titer_fname=None):

    #### DEFINE THE PRIORITY
    if titer_fname:
        HI_titer_count = count_titer_measurements(titer_fname)
        def priority(strain):
            return HI_titer_count[strain] + np.random.random()
    else:
        print("No titer counts provided - using random priorities")
        def priority(strain):
            return np.random.random()

    subcat_threshold = int(np.ceil(1.0*viruses_per_month/len(subcats)))

    virus_by_super_category, virus_by_category = populate_categories(metadata)
    def threshold_fn(x):
        #x is the subsampling category, in this case a tuple of (region, year, month)

        # if there are not enough viruses by super category, take everything
        if len(virus_by_super_category[x[1:]]) < viruses_per_month:
            return viruses_per_month

        # otherwise, sort sub categories by strain count
        sub_counts = sorted([(r, virus_by_super_category[(r, x[1], x[2])]) for r in subcats],
                             key=lambda y:len(y[1]))

        # if all (the smallest) subcat has more strains than the threshold, return threshold
        if len(sub_counts[0][1]) > subcat_threshold:
            return subcat_threshold


        strains_selected = 0
        tmp_subcat_threshold = subcat_threshold
        for ri, (r, strains) in enumerate(sub_counts):
            current_threshold = int(np.ceil(1.0*(viruses_per_month-strains_selected)/(len(subcats)-ri)))
            if r==x[0]:
                return current_threshold
            else:
                strains_selected += min(len(strains), current_threshold)
        return subcat_threshold

    selected_strains = []
    for cat, val in virus_by_category.items():
        if cat_valid(cat, time_interval):
            val.sort(key=priority, reverse=True)
            selected_strains.extend(val[:threshold_fn(cat)])

    return selected_strains


def cat_valid(cat, time_interval):
    if cat[-2]<time_interval[1].year or cat[-2]>time_interval[0].year:
        return False
    if (cat[-2]==time_interval[1].year and cat[-1]<time_interval[1].month) or\
       (cat[-2]==time_interval[0].year and cat[-1]>time_interval[0].month):
        return False
    return True

def determine_time_interval(time_interval, resolution):
    # determine date range to include strains from
    if time_interval: # explicitly specified
        datetime_interval = sorted([datetime.strptime(x, '%Y-%m-%d').date() for x in args.time_interval], reverse=True)
    else: # derived from resolution arguments (explicit takes precedence)
        if resolution:
            years_back = int(resolution[:-1])
        else:
            years_back = 3
        datetime_interval = [datetime.today().date(), (datetime.today()  - timedelta(days=365.25 * years_back)).date()]
    return datetime_interval

def parse_metadata(segments, metadata_files):
    metadata = {}
    for segment, fname in zip(segments, metadata_files):
        tmp_meta, columns = read_metadata(fname)

        numerical_dates = get_numerical_dates(tmp_meta, fmt='%Y-%m-%d')
        for x in tmp_meta:
            tmp_meta[x]['num_date'] = np.mean(numerical_dates[x])
            tmp_meta[x]['year'] = int(tmp_meta[x]['num_date'])
            tmp_meta[x]['month'] = int((tmp_meta[x]['num_date']%1)*12)
        metadata[segment] = tmp_meta
    return metadata

def parse_sequences(segments, sequence_files):
    """Load sequence names into a dictionary of sets indexed by segment.
    """
    sequences = {}
    for segment, filename in zip(segments, sequence_files):
        sequence_set = Bio.SeqIO.parse(filename, "fasta")
        sequences[segment] = set()
        for seq in sequence_set:
            sequences[segment].add(seq.name)

    return sequences

def summary(strains, metadata, segments, keys):
    print("Summary of strain counts by category")
    for segment in segments:
        print("Categories for segment", segment)
        category_counts = {}
        for strain in strains:
            if all([key in metadata[segment][strain] for key in keys]):
                category = "-".join([str(metadata[segment][strain][key]) for key in keys])
                if category not in category_counts:
                    category_counts[category] = 1
                else:
                    category_counts[category] += 1
        sorted_categories = sorted(category_counts.keys())
        for category in sorted_categories:
            print(category, category_counts[category])
        print("total", len(strains))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Select strains for downstream analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('-v', '--viruses_per_month', type = int, default=15,
                        help='Subsample x viruses per country per month. Set to 0 to disable subsampling.')
    parser.add_argument('--sequences', nargs='+', help="FASTA file with viral sequences, one for each segment")
    parser.add_argument('--metadata', nargs='+', help="file with metadata associated with viral sequences, one for each segment")
    parser.add_argument('--output', help="name of the file to write selected strains to")
    parser.add_argument('--verbose', action="store_true", help="turn on verbose reporting")

    parser.add_argument('-l', '--lineage', choices=['h3n2', 'h1n1pdm', 'vic', 'yam'], default='h3n2', type=str, help="single lineage to include (default: h3n2)")
    parser.add_argument('-r', '--resolution',default='3y', type = str,  help = "single resolution to include (default: 3y)")
    parser.add_argument('-s', '--segments', default=['ha'], nargs='+', type = str,  help = "list of segments to include (default: ha)")
    parser.add_argument('--sampling', default = 'even', type=str,
                        help='sample evenly over regions (even) (default), or prioritize one region (region name), otherwise sample randomly')
    parser.add_argument('--time-interval', nargs=2, help="explicit time interval to use -- overrides resolutions"
                                                                     " expects YYYY-MM-DD YYYY-MM-DD")
    parser.add_argument('--titers', help="a text file titers. this will only read in how many titer measurements are available for a each virus"
                                          " and use this count as a priority for inclusion during subsampling.")
    parser.add_argument('--include', help="a text file containing strains (one per line) that will be included regardless of subsampling")
    parser.add_argument('--max-include-range', type=float, default=5, help="number of years prior to the lower date limit for reference strain inclusion")
    parser.add_argument('--exclude', help="a text file containing strains (one per line) that will be excluded")

    args = parser.parse_args()
    time_interval = determine_time_interval(args.time_interval, args.resolution)

    # derive additional lower inclusion date for "force-included strains"
    lower_reference_cutoff = date(year = time_interval[1].year - args.max_include_range, month=1, day=1)
    upper_reference_cutoff = time_interval[0]

    # read strains to exclude
    excluded_strains = read_strain_list(args.exclude) if args.exclude else []
    # read strains to include
    included_strains = read_strain_list(args.include) if args.include else []

    # read in sequence names to determine which sequences already passed upstream filters
    sequence_names_by_segment = parse_sequences(args.segments, args.sequences)

    # read in meta data, parse numeric dates
    metadata = parse_metadata(args.segments, args.metadata)

    # eliminate all metadata entries that do not have sequences
    filtered_metadata = {}
    for segment in metadata:
        filtered_metadata[segment] = {}
        for name in metadata[segment]:
            if name in sequence_names_by_segment[segment]:
                filtered_metadata[segment][name] = metadata[segment][name]

    # filter down to strains with sequences for all required segments
    guide_segment = args.segments[0]
    strains_with_all_segments = set.intersection(*(set(filtered_metadata[x].keys()) for x in args.segments))
    # exclude outlier strains
    strains_with_all_segments.difference_update(set(excluded_strains))
    # subsample by region, month, year
    selected_strains = flu_subsampling({x:filtered_metadata[guide_segment][x] for x in strains_with_all_segments},
                                  args.viruses_per_month, time_interval, titer_fname=args.titers)

    # add strains that need to be included
    for strain in included_strains:
        if strain in strains_with_all_segments and strain not in selected_strains:
            # Do not include strains sampled too far in the past or strains
            # sampled from the future relative to the requested build interval.
            if (filtered_metadata[guide_segment][strain]['year'] >= lower_reference_cutoff.year and
                filtered_metadata[guide_segment][strain]['num_date'] <= numeric_date(upper_reference_cutoff)):
                selected_strains.append(strain)

    # summary of selected strains by region
    summary(selected_strains, filtered_metadata, args.segments, ['region'])

    # write the list of selected strains to file
    with open(args.output, 'w') as ofile:
        ofile.write('\n'.join(selected_strains))
