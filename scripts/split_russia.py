from augur.utils import read_lat_longs


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description="assign Eastern parts of russia to Europe",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str, help="metadata file")
    parser.add_argument('--output', type=str, help="output file")
    parser.add_argument('--locations', type=str, help="existing location file to complete")
    args = parser.parse_args()

    coord = read_lat_longs(overrides=args.locations)
    with open(args.metadata) as fh:
        d = fh.readlines()

    entries = d[0].split('\t')
    region_index = entries.index('region')
    division_index = entries.index('division')

    with open(args.output, 'wt') as fh:
        fh.write(d[0])
        for line in d:
            if 'Russia' in line:
                e = line.strip().split('\t')
                try:
                    print(coord[('division', e[division_index].lower())])
                except:
                    print(e[division_index])
                if ('division', e[division_index].lower()) in coord and coord[('division', e[division_index].lower())]['longitude']<59:
                    fh.write(line.replace('West Asia', 'Europe'))
                else:
                    fh.write(line)
            else:
                fh.write(line)

                            
                
                

