import glob
import json
import objectpath
import pygraphviz
import sqlite3

class ReachGraphAssembler(object):
    def __init__(self):
        self.graph = None
        #data_dir = '/home/ben/Dropbox/postdoc/cwc/src/bioagents/data/'
        #self.drug_db = sqlite3.connect(data_dir + 'drug_targets.db')

    def __del__(self):
        pass
        #if self.drug_db is not None:
        #    self.drug_db.close()
    
    #def find_target_drug(self, target_name):
    #    '''
    #    Find all the drugs that nominally target the target.
    #    '''
    #    target_name_root = target_name.split('/')[0]
    #    target_name_root = target_name_root.split('-')[0]
    #    res = self.drug_db.execute('SELECT name, synonyms, source_id FROM agent '
    #        'WHERE nominal_target LIKE "%%%s%%" ' % target_name_root).fetchall()
    #    drug_ids = [r[2] for r in res if r[2].startswith('HMSL')]
    #    ident_url = 'http://identifiers.org/lincs.molecule/'
    #    drug_urls = [ident_url + d.replace('HMSL','') + '-101' for d in drug_ids]
    #    return drug_urls

    def extend_graph(self, events, entities):
        if self.graph is None:
            self.graph = pygraphviz.AGraph(directed=True, rankdir='TB', splines='ortho')
        self.events_tree = objectpath.Tree(events)
        self.entities_tree = objectpath.Tree(entities)
     
        # Get regulations
        qstr = "$.frames[(@.type is 'regulation')]"
        res = self.events_tree.execute(qstr)
        for r in res:
            controlled_id  = [x for x in objectpath.Tree(r).execute("$.arguments[@.argument_label is \
                'controlled'].arg")][0]
            controlled_frame = [f for f in self.events_tree.execute("$.frames[@.frame_id is '%s']" % \
                controlled_id)][0]
            if controlled_frame['type'] == 'protein-modification':
                self.add_modification(r, controlled_frame)
            else:
                self.add_regulation(r)
        
        # Get activations
        qstr = "$.frames[(@.type is 'activation')]"
        res = self.events_tree.execute(qstr)
        for r in res:
            self.add_activation(r)

        # Get complex assembly
        qstr = "$.frames[(@.type is 'complex-assembly')]"
        res = self.events_tree.execute(qstr)
        for r in res:
            self.add_complex(r)
    
    @staticmethod
    def get_normalized_name(name):
        name = name.strip()
        name = name.upper()
        return name

    def get_grounding(self, node_id):
        qstr = "$.frames[(@.frame_id is \'%s\')].xrefs" % node_id
        xrefs = self.entities_tree.execute(qstr)
        return xrefs

    def add_node(self, node):
        node_id = node['arg']
        node_name = self.get_normalized_name(node['text'])
        xref = [xr for xr in self.get_grounding(node_id)][0][0]
        #drug_urls = self.find_target_drug(node_name)
        if xref['namespace'] == 'uazid':
            node_name_root = node_name.split('/')[0]
            node_name_root = node_name_root.split('-')[0]
            url = 'http://identifiers.org/hgnc.symbol/%s' % node_name_root
        else:
            url = 'http://identifiers.org/%s/%s' % (xref['namespace'], xref['id'])

        #if drug_urls:
        #    url = drug_urls[0]
        color = "#ffeeee"
        if node_name == 'EGF':
            node_label = 'Growth factors'
        else:
            node_label = node_name
        self.graph.add_node(node_name,
                        URL=url,
                        label=node_label,
                        shape='Mrecord',
                        fillcolor=color, style="filled", color="transparent",
                        fontsize="12",
                        margin="0.06,0")

    def add_edge(self, s, t, color='#000000', arrowhead='normal', style='solid'):
        if s is None or t is None:
            return
        else:
            self.add_node(s)
            self.add_node(t)
            s_name = self.get_normalized_name(s['text'])
            t_name = self.get_normalized_name(t['text'])
            self.graph.add_edge(s_name, t_name,
                                arrowhead=arrowhead, color=color, style=style)
        print '"%s" %d "%s"' % (s_name, (1 if style=='solid' else -1), t_name)


    def add_modification(self, controller_frame, controlled_frame):
        args = controller_frame['arguments']
        for a in args:
            if a['argument_label'] == 'controller':
                controller = a
        args = controlled_frame['arguments']
        for a in args:
            if a['argument_label'] == 'theme':
                controlled = a
        arrowhead = 'normal'
        
        if controller_frame['subtype'] == 'negative-regulation':
            color = '#cc0000'
        else:
            color = '#000000'
        arrowhead = 'dot'
      
        self.add_edge(controller, 
                        controlled,
                        arrowhead=arrowhead, color=color)

    def add_regulation(self, frame):
        args = frame['arguments']
        for a in args:
            if a['argument_label'] == 'controller':
                controller = a
            elif a['argument_label'] == 'controlled':
                controlled = a
        
        if frame['subtype'] == 'negative-regulation':
            color = '#cc0000'
        else:
            color = '#000000'
        arrowhead = 'diamond'
        
        self.add_edge(controller,
                        controlled,
                        arrowhead=arrowhead, color=color)

    def add_activation(self, frame):
        args = frame['arguments']
        for a in args:
            if a['argument_label'] == 'controller':
                controller = a
            elif a['argument_label'] == 'controlled':
                controlled = a
        if frame['subtype'] == 'negative-activation':
            color = '#cc0000'
            arrowhead = 'tee'
            style = 'dashed'
        else:
            color = '#000000'
            arrowhead = 'normal'
            style = 'solid'
        self.add_edge(controller,
                        controlled,
                        arrowhead=arrowhead, color=color, style=style)

    def add_complex(self, frame):
        args = frame['arguments']
        color = '#000000'
        arrowhead = 'none'
        style = 'solid'
        self.add_edge(args[0],args[1],
                        arrowhead=arrowhead, color=color, style=style)
    
    def get_string(self):
        graph_string = self.graph.string()
        graph_string = graph_string.replace('\\N', '\\n')
        return graph_string
    
if __name__ == '__main__':
    file_names = glob.glob('txt/*.json')
    ga = GraphAssembler()
    for f in file_names:
        fh = open(f, 'rt')
        txt = fh.read()
        txt = txt.replace('frame-id', 'frame_id')
        txt = txt.replace('argument-label', 'argument_label')

        json_dict = json.loads(txt)
        print f, len(json_dict['events']['frames'])
        fh.close()
        
        ga.extend_graph(json_dict['events'], json_dict['entities'])
    
    with open('graph.dot', 'wt') as fh:
        graph_str = ga.get_string()
        graph_str = graph_str.replace('\\N', '\\n')
        fh.write(graph_str)


