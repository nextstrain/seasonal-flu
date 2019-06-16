"""
This script runs the entire suite of seasonal flu builds via AWS batch
One batch run is created per lineage
"""

import subprocess
import argparse
import os

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run flu builds on aws', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--system', type = str, default = 'local', help='where to run, local or batch')
    parser.add_argument('-l', '--lineages', nargs='+', type = str,  help ="flu lineages to include", default=['h3n2', 'h1n1pdm'])
    parser.add_argument('-r', '--resolutions', nargs='+', type = str,  help ="flu resolutions to include", default=['2y'])
    parser.add_argument('-s', '--segments', nargs='+', type = str, help ="flu segments to include", default=['ha', 'na', 'pb2', 'pb1', 'pa', 'np', 'mp', 'ns'])
    params = parser.parse_args()

    if not os.path.exists("logs"):
        os.makedirs("logs")
    else:
        rmc = 'rm -rf logs/*'
        subprocess.call(rmc, shell=True)

    for lineage in params.lineages:
        cpus = len(params.resolutions) * len(params.segments)
        memory = 1800 * cpus
        if params.system == 'local':
            call = ['nextstrain', 'build', '.', '--jobs', '1']
        elif params.system == 'batch':
            call = ['nextstrain', 'build', '--aws-batch', '--aws-batch-cpus', str(cpus), '--aws-batch-memory', str(memory), '.', '--jobs', str(cpus)]
        targets = []
        for resolution in params.resolutions:
            for segment in params.segments:
                targets.append('auspice/seattle_flu_seasonal_%s_%s_%s_global_tree.json'%(lineage, segment, resolution))
        call.extend(targets)
        print(' '.join(call))
        log = open('logs/%s.txt'%(lineage), 'w')
        if params.system == 'local':
            pro = subprocess.call(call)
        if params.system == 'batch':
            pro = subprocess.Popen(call, stdout=log, stderr=log)
