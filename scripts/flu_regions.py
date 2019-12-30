region_properties = {
    "global":           {'label':"Global", "color":"#111111"},
    'Africa':           {'abbr':'AF', 'label':'Africa', 'popsize': 1.02, 'color':'#A0CCA5'},
    'Europe':           {'abbr':'EU', 'label':'Europe', 'popsize': 0.74, 'color':'#658447'},
    'North America':    {'abbr':'NA', 'label':'N America', 'popsize': 0.54, 'color':'#D6C568'},
    'China':            {'abbr':'CN', 'label':'China', 'popsize': 1.36, 'color':'#A76BB1'},
    'South Asia':       {'abbr':'SAS', 'label':'South Asia', 'popsize': 1.45, 'color':'#5199B7'},
    'Japan Korea':      {'abbr':'JK', 'label':'Japan/Korea', 'popsize': 0.20, 'color':'#2A4786'},
    'Oceania':          {'abbr':'OC', 'label':'Oceania', 'popsize': 0.04, 'color':'#8E1616'},
    'South America':    {'abbr':'SA', 'label':'S America', 'popsize': 0.41, 'color':'#EBA85F'},
    'Southeast Asia':   {'abbr':'SEA', 'label':'SE Asia', 'popsize': 0.62, 'color':'#8FBDD0'},
    'West Asia':        {'abbr':'WA', 'label':'W Asia', 'popsize': 0.75, 'color':'#76104B'},
}

region_names = [x for x in region_properties.keys() if x != 'global']
