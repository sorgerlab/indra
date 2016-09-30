import os
import pickle
from indra.preassembler.hierarchy_manager import hierarchies
from indra.preassembler import Preassembler, render_stmt_graph,\
                               flatten_evidence
from indra.mechlinker import MechLinker
from indra.assemblers import PysbAssembler, IndexCardAssembler,\
                             EnglishAssembler

def have_file(fname):
    return os.path.exists(fname)

def print_stmts(stmts, file_name):
    with open(file_name, 'wt') as fh:
        for s in stmts:
            agents = s.agent_list()
            db_refs = [('%s(%s)' % (a.name, a.db_refs)) 
                        for a in agents if a is not None]
            db_refs_str = (', '.join(db_refs)).encode('utf-8')
            try:
                fh.write('%s\t%s\t%s\n' %
                        (s, db_refs_str, s.evidence[0].text.encode('utf-8')))
            except UnicodeDecodeError:
                fh.write('%s\t%s\t%s\n' %
                        (s, db_refs_str, s.evidence[0].text))

def is_protein_or_chemical(agent):
    # Default is True if agent is None
    if agent is None:
        return True
    if agent.db_refs.get('UP') is not None or \
        agent.db_refs.get('HGNC') is not None or \
        agent.db_refs.get('CHEBI') is not None or \
        agent.db_refs.get('PFAM-DEF') is not None:
        return True
    return False

background_secs = ['abstract', 'introduction', 'background']

def is_background_knowledge(stmt):
    any_background = False
    for ev in stmt.evidence:
        epi = ev.epistemics
        if epi is not None:
            sec = epi.get('section_type')
            # If there is at least one evidence not from a 
            # background section then we consider this to be
            # a non-background knowledge finding.
            if sec is not None and sec not in background_secs:
                return False
            # If there is at least one evidence that is explicitly
            # from a background section then we keep track of that.
            elif sec in background_secs:
                any_background = True
    # If there is any explicit evidence for this statement being
    # background info (and no evidence otherwise) then we return
    # True, otherwise (for instnace of there is no section info at all)
    # we return False.
    return any_background

def run_assembly(stmts, folder, pmcid):
    indexcard_prefix = folder + '/index_cards/' + pmcid
    otherout_prefix = folder + '/other_outputs/' + pmcid

    # Filter for grounding
    grounded_stmts = []
    for st in stmts:
        if all([is_protein_or_chemical(a) for a in st.agent_list()]):
            grounded_stmts.append(st)

    # Instantiate the Preassembler
    pa = Preassembler(hierarchies)

    pa.add_statements(grounded_stmts)
    print '%d statements collected in total.' % len(pa.stmts)
    unique_stmts = pa.combine_duplicates()
    print '%d statements after combining duplicates.' % len(unique_stmts)
    ml = MechLinker(unique_stmts)
    ml.link_statements()
    pa = Preassembler(hierarchies, ml.statements)
    pa.combine_duplicates()
    related_stmts = pa.combine_related()
    print '%d statements after combining related.' % len(related_stmts)

    with open(otherout_prefix + '.pkl', 'wb') as fh:
        pickle.dump(related_stmts, fh)

    flattened_evidence_stmts = flatten_evidence(related_stmts)

    card_counter = 1
    card_lim = float('inf')
    top_stmts = []
    for st in sorted(flattened_evidence_stmts,
                     key=lambda x: len(x.evidence), reverse=True):
        print len(st.evidence), st

        if is_background_knowledge(st):
            print 'This statement is background knowledge - skipping.'
            continue
        # Assemble IndexCards
        ia = IndexCardAssembler([st])
        ia.make_model()
        if ia.cards:
            ia.save_model(indexcard_prefix + '-%d.json' % card_counter)
            card_counter += 1
            top_stmts.append(st)
            if card_counter > card_lim:
                break

    ea = EnglishAssembler(top_stmts)
    print '======================='
    print ea.make_model()
    print '======================='

    # Print the statement graph
    graph = render_stmt_graph(related_stmts)
    graph.draw(otherout_prefix + '_graph.pdf', prog='dot')
    # Print statement diagnostics
    print_stmts(pa.stmts, otherout_prefix + '_statements.tsv')
    print_stmts(related_stmts, otherout_prefix + '_related_statements.tsv')

    pya = PysbAssembler()
    pya.add_statements(related_stmts)
    model = pya.make_model()

    print 'PySB model has %d monomers and %d rules' %\
        (len(model.monomers), len(model.rules))

