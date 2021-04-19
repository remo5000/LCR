import pydot
from itertools import count
from collections import defaultdict

# To prepare the file for parsing, run these sed substitutions
# To remove dashes in names (messes with parsing)
# sed -E -e ':loop' -e 's/\->(.*)-(.*)\[level/->\1_\2\[level/' -e 't loop' -i .backup robots_net-graph-2014-07-07.dot
# sed -E -e ':loop' -e 's/(.*)-(.*)-> /\1_\2-> /' -e 't loop' -i .backup robots_net-graph-2014-07-07.dot
# To remove spaces in names (messes with labels)
# sed -E -e ':loop' -e 's/-> (.*) (.*) \[level/-> \1\2 \[level/' -e 't loop' -i .backup robots_net-graph-2014-07-07.dot
# sed -E -e ':loop' -e 's/^(.*) (.*) -> /\1\2 -> /' -e 't loop' -i .backup robots_net-graph-2014-07-07.dot
# And test using a python interpreter if it's parsable.
print('Reading .dot file...')
graphs = pydot.graph_from_dot_file('robots_net-graph-2014-07-07.dot')

g = graphs[0]
m = {'Master': 0, 'Apprentice': 1, 'Journeyer': 2, 'Observer': 3}

c = count(1)
id_of_user = defaultdict(lambda: next(c))

print('Writing out to .edge file...')
with open('robots.edge', 'w') as out:
    for edge in g.get_edge_list():
        u, v = edge.obj_dict['points']

        u = id_of_user[u]
        v = id_of_user[v]

        
        if 'level' not in edge.obj_dict['attributes']:
            print(edge.obj_dict)
            continue
            ls = m['Apprentice']
        label_string = edge.obj_dict['attributes']['level'].replace('"', '')
        ls = m[label_string]

        out.write(f'{u} {v} {ls}\n')

