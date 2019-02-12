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
                call = ['nextstrain', 'build', '.', '-j', '1']
            elif params.system == 'batch':
                call = ['nextstrain', 'build', '--aws-batch', '.', '-j', '4'] # 4 jobs to match xlarge nodes
            targets = []
            for resolution in params.resolutions:
                for segment in params.segments:
                    targets.append('auspice/flu_seasonal_%s_%s_%s_tree.json'%(lineage, segment, resolution))
                    targets.append('auspice/flu_seasonal_%s_%s_%s_meta.json'%(lineage, segment, resolution))
                    targets.append('auspice/flu_seasonal_%s_%s_%s_tip-frequencies.json'%(lineage, segment, resolution))
            call.extend(targets)
            print(' '.join(call))
            log = open('logs/live_flu_%s.txt'%(lineage), 'w')
            if params.system == 'local':
                pro = subprocess.call(call)
            if params.system == 'batch':
                pro = subprocess.Popen(call, stdout=log, stderr=log)

    if params.version == 'who' or params.version == 'both':
        for center in params.centers:
            for lineage in params.lineages:
                if params.system == 'local':
                    call = ['nextstrain', 'build', '.', '-s', 'Snakefile_WHO', '-j', '1']
                elif params.system == 'batch':
                    call = ['nextstrain', 'build', '--aws-batch', '.', '-s', 'Snakefile_WHO', '-j', '4'] # 4 jobs to match xlarge nodes
                targets = []
                segment = 'ha'
                resolutions = [r for r in params.resolutions if r == '2y' or r == '6y']
                for resolution in resolutions:
                    for passage in params.passages:
                        assays = [assay for assay in params.assays if lineage == 'h3n2' or assay == 'hi']
                        for assay in assays:
                            targets.append('auspice-who/flu_%s_%s_%s_%s_%s_%s_tree.json'%(center, lineage, segment, resolution, passage, assay))
                            targets.append('auspice-who/flu_%s_%s_%s_%s_%s_%s_meta.json'%(center, lineage, segment, resolution, passage, assay))
                            targets.append('auspice-who/flu_%s_%s_%s_%s_%s_%s_entropy.json'%(center, lineage, segment, resolution, passage, assay))
                            targets.append('auspice-who/flu_%s_%s_%s_%s_%s_%s_frequencies.json'%(center, lineage, segment, resolution, passage, assay))
                            targets.append('auspice-who/flu_%s_%s_%s_%s_%s_%s_sequences.json'%(center, lineage, segment, resolution, passage, assay))
                            targets.append('auspice-who/flu_%s_%s_%s_%s_%s_%s_titer-sub-model.json'%(center, lineage, segment, resolution, passage, assay))
                            targets.append('auspice-who/flu_%s_%s_%s_%s_%s_%s_titer-tree-model.json'%(center, lineage, segment, resolution, passage, assay))
                            targets.append('auspice-who/flu_%s_%s_%s_%s_%s_%s_titers.json'%(center, lineage, segment, resolution, passage, assay))
                call.extend(targets)
                print(' '.join(call))
                log = open('logs/who_flu_%s_%s.txt'%(center, lineage), 'w')
                if params.system == 'local':
                    pro = subprocess.call(call)
                if params.system == 'batch':
                    pro = subprocess.Popen(call, stdout=log, stderr=log)
