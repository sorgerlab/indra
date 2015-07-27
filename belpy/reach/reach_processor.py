import objectpath

from belpy.statements import *


class ReachProcessor:
    def __init__(self, json_dict):
        self.tree = objectpath.Tree(json_dict)
        self.statements = []

    def get_phosphorylation(self):
        citation = self.tree.execute("$.frames.object-meta.doc-id")
        qstr = "$.frames[(@.type is 'protein-modification') " + \
               "and (@.subtype is 'phosphorylation')]"
        res = self.tree.execute(qstr)
        for r in res:
            frame_id = r['frame-id']
            args = r['arguments']
            site = None
            theme = None
            controller = None

            for a in args:
                if a['argument-label'] == 'theme':
                    theme = a['text']
                elif a['argument-label'] == 'site':
                    site = a['text']
            qstr = "$.frames[(@.type is 'regulation') and " + \
                   "(@.arguments[0].arg is '%s')]" % frame_id
            reg_res = self.tree.execute(qstr)
            for reg in reg_res:
                for a in reg['arguments']:
                    if a['argument-label'] == 'controller':
                        controller = a['text']

            controller_agent = Agent(controller)
            theme_agent = Agent(theme)
            mod = 'Phosphorylation'
            pos = site
            sentence = r['verbose-text']
            evidence = sentence
            # TODO: read $.object-meta.doc-id as citation
            # but dashes don't work with objectpath!
            citation = ''
            annotations = None
            self.statements.append(Phosphorylation(controller_agent,
                                   theme_agent, mod, pos, sentence,
                                   citation, evidence, annotations))
