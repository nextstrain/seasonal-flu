"""
This script runs the entire suite of seasonal flu builds via AWS batch
One batch run is created per HA/NA combination, ie
h3n2_ha_2y + h3n2_na_2y is one build and
h3n2_ha_6y + h3n2_na_6y is another
"""

import subprocess
import argparse
import os

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run flu builds on aws', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--system', type = str, default = 'local', help='where to run, local or batch')
    parser.add_argument('-v', '--version', type = str, default = 'live', help='version to run, live or who or both')
    parser.add_argument('-l', '--lineages', nargs='+', type = str,  help ="flu lineages to include", default=['h3n2', 'h1n1pdm', 'vic', 'yam'])
    parser.add_argument('-r', '--resolutions', nargs='+', type = str,  help ="flu resolutions to include", default=['2y', '3y', '6y' , '12y'])
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

    if params.version == 'live' or params.version == 'both':
        for lineage in params.lineages:
            if params.system == 'local':
                call = ['nextstrain', 'build', '.', '--jobs', '1']
            elif params.system == 'batch':
                call = ['nextstrain', 'build', '--aws-batch', '--aws-batch-cpus', '8', '--aws-batch-memory', '15200', '.', '--jobs', '8']
            targets = []
            for resolution in params.resolutions:
                for segment in params.segments:
                    targets.append('targets/flu_seasonal_%s_%s_%s'%(lineage, segment, resolution))
            call.extend(targets)
            print(' '.join(call))
            log = open('logs/live_%s.txt'%(lineage), 'w')
            if params.system == 'local':
                pro = subprocess.call(call)
            if params.system == 'batch':
                pro = subprocess.Popen(call, stdout=log, stderr=log)

    if params.version == 'who' or params.version == 'both':
        for lineage in params.lineages:
            if params.system == 'local':
                call = ['nextstrain', 'build', '.', '--snakefile', 'Snakefile_WHO', '--jobs', '1']
            elif params.system == 'batch':
                call = ['nextstrain', 'build', '--aws-batch', '--aws-batch-cpus', '16', '--aws-batch-memory', '31000', '.', '--snakefile', 'Snakefile_WHO', '--jobs', '16']
            targets = []
            segment = 'ha'
            resolutions = [r for r in params.resolutions if r == '2y' or r == '6y']
            for center in params.centers:
                for resolution in resolutions:
                    for passage in params.passages:
                        assays = [assay for assay in params.assays if lineage == 'h3n2' or assay == 'hi']
                        for assay in assays:
                            targets.append('targets/flu_%s_%s_%s_%s_%s_%s'%(center, lineage, segment, resolution, passage, assay))
            call.extend(targets)
            print(' '.join(call))
            log = open('logs/who_%s.txt'%(lineage), 'w')
            if params.system == 'local':
                pro = subprocess.call(call)
            if params.system == 'batch':
                pro = subprocess.Popen(call, stdout=log, stderr=log)
