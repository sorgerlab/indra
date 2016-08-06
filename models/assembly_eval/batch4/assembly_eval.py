import os
import pickle
from indra.preassembler.hierarchy_manager import hierarchies
from indra.preassembler import Preassembler, render_stmt_graph,\
                               flatten_evidence
from indra.mechlinker import MechLinker
from indra.assemblers import PysbAssembler, IndexCardAssembler,\
                             EnglishAssembler
from indra.belief import BeliefEngine

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
                continue
            except Exception as e:
                pass

def is_protein_or_chemical(agent):
    '''Return True if the agent is a protein/protein family or chemical.'''
    # Default is True if agent is None
    if agent is None:
        return True
    if agent.db_refs.get('UP') is not None or \
        agent.db_refs.get('HGNC') is not None or \
        agent.db_refs.get('CHEBI') is not None or \
        agent.db_refs.get('PFAM-DEF') is not None:
        return True
    return False

# Sections of the paper that we consider to be background info
background_secs = ['abstract', 'introduction', 'background']

def is_background_knowledge(stmt):
    '''Return True if the Statement is only supported by background knowledge.'''
    any_background = False
    # Iterate over all evidence for the statement
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

def multiple_sources(stmt):
    '''Return True if statement is supported by multiple sources.

    Note: this is currently not used and replaced by BeliefEngine score cutoff
    '''
    sources = list(set([e.source_api for e in stmt.evidence]))
    if len(sources) > 1:
        return True
    return False

def run_assembly(stmts, folder, pmcid):
    '''Run assembly on a list of statements, for a given PMCID.'''
    # Folder for index card output (scored submission)
    indexcard_prefix = folder + '/index_cards/' + pmcid
    # Folder for other outputs (for analysis, debugging)
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

    # Combine duplicates
    unique_stmts = pa.combine_duplicates()
    print '%d statements after combining duplicates.' % len(unique_stmts)

    # Run BeliefEngine on unique statements
    epe = BeliefEngine(pa.unique_stmts)
    epe.set_prior_probs(pa.unique_stmts)

    # Build statement hierarchy
    related_stmts = pa.combine_related()
    # Run BeliefEngine on hierarchy
    epe.set_hierarchy_probs(related_stmts)
    print '%d statements after combining related.' % len(related_stmts)

    # Instantiate the mechanism linker
    ml = MechLinker(related_stmts)
    # Link statements
    linked_stmts = ml.link_statements()
    # Run BeliefEngine on linked statements
    epe.set_linked_probs(linked_stmts)
    # Print linked statements for debugging purposes
    print 'Linked\n====='
    for ls in linked_stmts:
        print ls.inferred_stmt, ls.inferred_stmt.belief
    print '============='

    # Combine all statements including linked ones
    all_statements = ml.statements + [ls.inferred_stmt for ls in linked_stmts]

    # Instantiate a new preassembler
    pa = Preassembler(hierarchies, all_statements)
    # Build hierarchy again
    pa.combine_duplicates()
    # Choose the top-level statements
    related_stmts = pa.combine_related()

    # Dump top-level statements in a pickle
    with open(otherout_prefix + '.pkl', 'wb') as fh:
        pickle.dump(related_stmts, fh)

    # Flatten evidence for statements
    flattened_evidence_stmts = flatten_evidence(related_stmts)

    # Start a card counter
    card_counter = 1
    # We don't limit the number of cards reported in this round
    card_lim = float('inf')
    top_stmts = []
    ###############################################
    # The belief cutoff for statements
    belief_cutoff = 0.3
    ###############################################
    # Sort by amount of evidence
    for st in sorted(flattened_evidence_stmts,
                     key=lambda x: x.belief, reverse=True):
        if belief >= belief_cutoff:
            print belief, st
        if belief < belief_cutoff:
            print 'SKIP', belief, st

        # If it's background knowledge, we skip the statement
        if is_background_knowledge(st):
            print 'This statement is background knowledge - skipping.'
            continue

        # Assemble IndexCards
        ia = IndexCardAssembler([st])
        ia.make_model()
        # If the index card was actually made 
        # (not all statements can be assembled into index cards to
        # this is often not the case)
        if ia.cards:
            # Save the index card json
            ia.save_model(indexcard_prefix + '-%d.json' % card_counter)
            card_counter += 1
            top_stmts.append(st)
            if card_counter > card_lim:
                break

    # Print the English-assembled model for debugging purposes
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
