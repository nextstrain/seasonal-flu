"""
This script runs the entire suite of seasonal flu builds via AWS batch
One batch run for 'live' builds and a separate batch run is created for 'who' builds
"""

import subprocess
import argparse
import os

def get_cpus(jobs):
    count = 1
    if jobs >= 72:
        count = 72
    elif jobs >= 36:
        count = 36
    elif jobs >= 16:
        count = 16
    elif jobs >= 8:
        count = 8
    elif jobs >= 4:
        count = 4
    elif jobs >= 2:
        count = 2
    return count

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run flu builds on aws', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--system', type = str, default = 'local', help='where to run, local or batch')
    parser.add_argument('-v', '--version', type = str, default = 'live', help='version to run, live or who')
    parser.add_argument('-l', '--lineages', nargs='+', type = str,  help ="flu lineages to include", default=['h3n2', 'h1n1pdm', 'vic', 'yam'])
    parser.add_argument('-r', '--resolutions', nargs='+', type = str,  help ="flu resolutions to include", default=['6m', '2y', '3y', '6y' , '12y'])
    parser.add_argument('-s', '--segments', nargs='+', type = str, help ="flu segments to include", default=['ha', 'na'])
    parser.add_argument('-c', '--centers', nargs='+', type = str, help ="who collaborating centers to include", default=['cdc', 'who', 'niid', 'vidrl', 'crick'])
    parser.add_argument('-p', '--passages', nargs='+', type = str, help ="passages to include", default=['cell', 'egg'])
    parser.add_argument('-a', '--assays', nargs='+', type = str, help ="titer assays to include", default=['hi', 'fra'])
    params = parser.parse_args()

    if not os.path.exists("logs"):
        os.makedirs("logs")
    else:
        rmc = 'rm -rf logs/*'
        subprocess.call(rmc, shell=True)

    if params.version == 'live':

        targets = []
        for lineage in params.lineages:
            for resolution in params.resolutions:
                for segment in params.segments:
                    targets.append('targets/flu_seasonal_%s_%s_%s'%(lineage, segment, resolution))

        cpus = get_cpus(len(targets))
        memory = 1800 * cpus
        if params.system == 'local':
            call = ['nextstrain', 'build', '.', '--jobs', '1']
        elif params.system == 'batch':
            call = ['nextstrain', 'build', '--aws-batch', '--detach', '--aws-batch-cpus', str(cpus), '--aws-batch-memory', str(memory), '.', '--jobs', str(cpus)]
        call.extend(targets)

        if targets:
            print(' '.join(call))
            pro = subprocess.call(call)


    if params.version == 'who':

        segment = 'ha'
        resolutions = [r for r in params.resolutions if r == '2y' or r == '6y']
        assays = [a for a in params.assays if a == 'hi']

        targets = []
        for lineage in params.lineages:
            for center in params.centers:
                for resolution in resolutions:
                    for passage in params.passages:
                        for assay in assays:
                            targets.append('targets/flu_%s_%s_%s_%s_%s_%s'%(center, lineage, segment, resolution, passage, assay))
        if 'h3n2' in params.lineages and 'fra' in params.assays:
            lineage = 'h3n2'
            assay = 'fra'
            for center in params.centers:
                for resolution in resolutions:
                    for passage in params.passages:
                        targets.append('targets/flu_%s_%s_%s_%s_%s_%s'%(center, lineage, segment, resolution, passage, assay))

        cpus = get_cpus(len(targets))
        memory = 1800 * cpus
        if params.system == 'local':
            call = ['nextstrain', 'build', '.', '--snakefile', 'Snakefile_WHO', '--jobs', '1']
        elif params.system == 'batch':
            call = ['nextstrain', 'build', '--aws-batch', '--detach', '--aws-batch-cpus', str(cpus), '--aws-batch-memory', str(memory), '.', '--snakefile', 'Snakefile_WHO', '--jobs', str(cpus)]
        call.extend(targets)

        if targets:
            print(' '.join(call))
            pro = subprocess.call(call)
