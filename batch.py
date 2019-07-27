"""
This script runs the entire suite of seasonal flu builds via AWS batch
One batch run for 'live' builds and a separate batch run is created for 'who' builds
"""

import subprocess
import argparse
import os

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run flu builds on aws', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--system', type = str, default = 'local', help='where to run, local or batch')
    parser.add_argument('-v', '--version', type = str, default = 'live', help='version to run, live or who or both')
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

    if params.version == 'live' or params.version == 'both':
        cpus = len(params.lineages) * len(params.resolutions) * len(params.segments)
        if cpus > 36:
            cpus = 36
        memory = 1800 * cpus
        if params.system == 'local':
            call = ['nextstrain', 'build', '.', '--jobs', '1']
        elif params.system == 'batch':
            call = ['nextstrain', 'build', '--aws-batch', '--aws-batch-cpus', str(cpus), '--aws-batch-memory', str(memory), '.', '--jobs', str(cpus)]
        targets = []
        for lineage in params.lineages:
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
        segment = 'ha'
        resolutions = [r for r in params.resolutions if r == '2y' or r == '6y']
        assay = 'hi'
        cpus = len(params.lineages) * len(params.centers) * len(resolutions) * len(params.passages)
        if cpus > 72:
            cpus = 72
        memory = 1800 * cpus
        if params.system == 'local':
            call = ['nextstrain', 'build', '.', '--snakefile', 'Snakefile_WHO', '--jobs', '1']
        elif params.system == 'batch':
            call = ['nextstrain', 'build', '--aws-batch', '--aws-batch-cpus', str(cpus), '--aws-batch-memory', str(memory), '.', '--snakefile', 'Snakefile_WHO', '--jobs', str(cpus)]
        targets = []
        for lineage in params.lineages:
            for center in params.centers:
                for resolution in resolutions:
                    for passage in params.passages:
                        targets.append('targets/flu_%s_%s_%s_%s_%s_%s'%(center, lineage, segment, resolution, passage, assay))
        if 'h3n2' in params.lineages and 'fra' in params.assays:
            assay = 'fra'
            for center in params.centers:
                for resolution in resolutions:
                    for passage in params.passages:
                        targets.append('targets/flu_%s_%s_%s_%s_%s_%s'%(center, lineage, segment, resolution, passage, assay))
        call.extend(targets)
        print(' '.join(call))
        log = open('logs/who_%s.txt'%(lineage), 'w')
        if params.system == 'local':
            pro = subprocess.call(call)
        if params.system == 'batch':
            pro = subprocess.Popen(call, stdout=log, stderr=log)
