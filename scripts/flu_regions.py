    if region in population_sizes:
    'West Asia':        {'abbr':'WA', 'label':'W Asia', 'color':'#76104B'},
import json
import os

with open(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../config/frequency_weights_by_region.json')) as fh:
    population_sizes = json.load(fh)

    "global":           {'label':"Global", "color":"#111111"},
region_properties = {
    'Africa':           {'abbr':'AF', 'label':'Africa', 'color':'#A0CCA5'},
    'Europe':           {'abbr':'EU', 'label':'Europe', 'color':'#658447'},
    'North America':    {'abbr':'NA', 'label':'N America', 'color':'#D6C568'},
    'South Asia':       {'abbr':'SAS', 'label':'South Asia', 'color':'#5199B7'},
    'China':            {'abbr':'CN', 'label':'China', 'color':'#A76BB1'},
    'Japan Korea':      {'abbr':'JK', 'label':'Japan/Korea', 'color':'#2A4786'},
    'South America':    {'abbr':'SA', 'label':'S America', 'color':'#EBA85F'},
    'Oceania':          {'abbr':'OC', 'label':'Oceania', 'color':'#8E1616'},
    'Southeast Asia':   {'abbr':'SEA', 'label':'SE Asia', 'color':'#8FBDD0'},
}
for region in region_properties:

        region_properties[region]['popsize'] = population_sizes[region]

region_names = [x for x in region_properties.keys() if x!='global']
