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
    parser.add_argument('-l', '--lineages', nargs='+', type = str,  help ="flu lineages to include", default=['h3n2', 'h1n1pdm', 'vic', 'yam'])
    parser.add_argument('-r', '--resolutions', nargs='+', type = str,  help ="flu resolutions to include", default=['2y', '6y'])
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

    for center in params.centers:
        for lineage in params.lineages:
            for resolution in params.resolutions:
                for passage in params.passages:
                    if lineage == "h3n2":
                        assays = params.assays
                        for assay in assays:
                            call = ['nextstrain', 'build', '--aws-batch', '.', '-s Snakefile_WHO', '-j 2']
                            targets = []
                            for segment in params.segments:
                                # auspice-who/flu_cdc_h3n2_ha_2y_cell_hi_tree.json, auspice-who/flu_cdc_h3n2_ha_2y_egg_hi_tree.json, auspice-who/flu_cdc_h3n2_ha_2y_cell_fra_tree.json, auspice-who/flu_cdc_h3n2_ha_2y_egg_fra_tree.json, auspice-who/flu_cdc_h3n2_ha_6y_cell_hi_tree.json, auspice-who/flu_cdc_h3n2_ha_6y_egg_hi_tree.json, auspice-who/flu_cdc_h3n2_ha_6y_cell_fra_tree.json, auspice-who/flu_cdc_h3n2_ha_6y_egg_fra_tree.json, auspice-who/flu_cdc_h3n2_na_2y_cell_hi_tree.json
                                targets.append('auspice-who/flu_%s_%s_%s_%s_%s_%s_tree.json'%(center, lineage, segment, resolution, passage, assay))
                                targets.append('auspice-who/flu_%s_%s_%s_%s_%s_%s_meta.json'%(center, lineage, segment, resolution, passage, assay))
                                targets.append('auspice-who/flu_%s_%s_%s_%s_%s_%s_entropy.json'%(center, lineage, segment, resolution, passage, assay))
                                # targets.append('auspice-who/flu_%s_%s_%s_%s_%s_%s_sequences.json'%(center, lineage, segment, resolution, passage, assay))
                                targets.append('auspice-who/flu_%s_%s_%s_%s_%s_%s_titer_subs_model.json'%(center, lineage, segment, resolution, passage, assay))
                                targets.append('auspice-who/flu_%s_%s_%s_%s_%s_%s_titer_tree_model.json'%(center, lineage, segment, resolution, passage, assay))
                                targets.append('auspice-who/flu_%s_%s_%s_%s_%s_%s_titers.json'%(center, lineage, segment, resolution, passage, assay))

                            call.extend(targets)
                            print(' '.join(call))
                            log = open('logs/%s_%s.txt'%(lineage, resolution), 'w')
                            # run in background with subprocess.Popen
                            pro = subprocess.Popen(call, stdout=log, stderr=log)
                    else:
                        assay = 'hi'
                        call = ['nextstrain', 'build', '--aws-batch', '.', '-s Snakefile_WHO', '-j 2']
                        targets = []
                        for segment in params.segments:
                            targets.append('auspice/flu_%s_%s_%s_%s_%s_%s_tree.json'%(center, lineage, segment, resolution, passage, assay))
                            targets.append('auspice-who/flu_%s_%s_%s_%s_%s_%s_meta.json'%(center, lineage, segment, resolution, passage, assay))
                            targets.append('auspice-who/flu_%s_%s_%s_%s_%s_%s_entropy.json'%(center, lineage, segment, resolution, passage, assay))
                            # targets.append('auspice-who/flu_%s_%s_%s_%s_%s_%s_sequences.json'%(center, lineage, segment, resolution, passage, assay))
                            targets.append('auspice-who/flu_%s_%s_%s_%s_%s_%s_titer_subs_model.json'%(center, lineage, segment, resolution, passage, assay))
                            targets.append('auspice-who/flu_%s_%s_%s_%s_%s_%s_titer_tree_model.json'%(center, lineage, segment, resolution, passage, assay))
                            targets.append('auspice-who/flu_%s_%s_%s_%s_%s_%s_titers.json'%(center, lineage, segment, resolution, passage, assay))

                        call.extend(targets)
                        print(' '.join(call))
                        log = open('logs/%s_%s.txt'%(lineage, resolution), 'w')
                        # run in background with subprocess.Popen
                        pro = subprocess.Popen(call, stdout=log, stderr=log)
