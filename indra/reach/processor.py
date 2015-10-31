import re
import objectpath

from indra.statements import *

residue_names = {
    'S': 'Serine',
    'T': 'Threonine',
    'Y': 'Tyrosine',
    'SER': 'Serine',
    'THR': 'Threonine',
    'TYR': 'Tyrosine',
    'SERINE': 'Serine',
    'THREONINE': 'Threonine',
    'TYROSINE': 'Tyrosine'
    }


class ReachProcessor:
    def __init__(self, json_dict):
        self.tree = objectpath.Tree(json_dict)
        self.statements = []

    def get_phosphorylation(self):
        citation = self.tree.execute("$.frames.object_meta.doc_id")
        qstr = "$.frames[(@.type is 'protein-modification') " + \
               "and (@.subtype is 'phosphorylation')]"
        res = self.tree.execute(qstr)
        for r in res:
            frame_id = r['frame_id']
            args = r['arguments']
            site = None
            theme = None
            controller = None

            for a in args:
                if a['argument_label'] == 'theme':
                    theme = a['text']
                elif a['argument_label'] == 'site':
                    site = a['text']
            qstr = "$.frames[(@.type is 'regulation') and " + \
                   "(@.arguments[0].arg is '%s')]" % frame_id
            reg_res = self.tree.execute(qstr)
            controller = None
            for reg in reg_res:
                for a in reg['arguments']:
                    if a['argument_label'] == 'controller':
                        controller = a['text']

            if controller is None:
                warnings.warn('Skipping phosphorylation with missing controller.')
                continue

            controller_agent = Agent(controller)
            theme_agent = Agent(theme)
            mod = 'Phosphorylation'
            if site is not None:
                residue, pos = self._parse_site_text(site)
            else:
                residue = ''
                pos = ''
            mod = mod + residue
            sentence = r['verbose-text']
            evidence = sentence
            # TODO: read $.object-meta.doc-id as citation
            # but dashes don't work with objectpath!
            citation = ''
            annotations = None
            self.statements.append(Phosphorylation(controller_agent,
                                   theme_agent, mod, pos, sentence,
                                   citation, evidence, annotations))
    
    def get_complexes(self):
        citation = self.tree.execute("$.frames.object_meta.doc_id")
        qstr = "$.frames[@.type is 'complex-assembly']"
        res = self.tree.execute(qstr)
        for r in res:
            frame_id = r['frame_id']
            args = r['arguments']
            members = []
            for a in args:
                agent = Agent(a['text'])
                members.append(agent)
            self.statements.append(Complex(members))
    
    def _parse_site_text(self, s):
        m = re.match(r'([TYS])[-]?([0-9]+)', s)
        if m is not None:
            residue = residue_names[m.groups()[0]]
            site = m.groups()[1]
            return residue, site

        m = re.match(r'(THR|TYR|SER)[- ]?([0-9]+)', s.upper())    
        if m is not None:
            residue = residue_names[m.groups()[0]]
            site = m.groups()[1]
            return residue, site

        m = re.match(r'(THREONINE|TYROSINE|SERINE)[^0-9]*([0-9]+)', s.upper())
        if m is not None:
            residue = residue_names[m.groups()[0]]
            site = m.groups()[1]
            return residue, site

        m = re.match(r'.*(THREONINE|TYROSINE|SERINE).*', s.upper())
        if m is not None:
            residue = residue_names[m.groups()[0]]
            site = None
            return residue, site
       
        return '', None
