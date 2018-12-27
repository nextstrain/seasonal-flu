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

    parser = argparse.ArgumentParser(description='Run flu builds on aws')
    parser.add_argument('-l', '--lineages', nargs='+', type = str,  help ="flu lineages to include")
    parser.add_argument('-r', '--resolutions', nargs='+', type = str,  help ="flu resolutions to include")
    params = parser.parse_args()

    segments = ["ha", "na"]

    if not os.path.exists("logs"):
        os.makedirs("logs")

    if params.lineages is None:
        params.lineages = ['h3n2', 'h1n1pdm', 'vic', 'yam']

    if params.resolutions is None:
        params.resolutions = ['2y', '3y', '6y', '12y']

    for lineage in params.lineages:
        for resolution in params.resolutions:
            call = ['nextstrain', 'build', '--aws-batch', '.', '-j 2']
            targets = []
            for segment in segments:
                targets.append('auspice/flu_seasonal_%s_%s_%s_tree.json'%(lineage, segment, resolution))
                targets.append('auspice/flu_seasonal_%s_%s_%s_meta.json'%(lineage, segment, resolution))
                targets.append('auspice/flu_seasonal_%s_%s_%s_tip-frequencies.json'%(lineage, segment, resolution))
            call.extend(targets)
            print(' '.join(call))
            log = open('logs/%s_%s.txt'%(lineage, resolution), 'w')
            # run in background with subprocess.Popen
            pro = subprocess.Popen(call, stdout=log, stderr=log)
