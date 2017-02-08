from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from os.path import join as pjoin
import os.path
from pysb import Observable
from pysb.export.kappa import KappaExporter
from indra.util import read_unicode_csv
from indra.assemblers import PysbAssembler, IndexCardAssembler
from indra.statements import *
import indra.tools.assemble_corpus as ac

def assemble_pysb(stmts, data_genes, out_file):
    """Return an assembled PySB model."""
    # Filter the INDRA Statements to be put into the model
    stmts = ac.filter_direct(stmts)
    stmts = ac.filter_belief(stmts, 0.95)
    stmts = ac.filter_top_level(stmts)
    stmts = ac.filter_gene_list(stmts, data_genes, 'all')
    stmts = filter_enzyme_kinase(stmts)
    stmts = filter_mod_nokinase(stmts)
    stmts = filter_transcription_factor(stmts)
    # Simplify activity types
    stmts = ac.reduce_activities(stmts)

    # Assemble model
    pa = PysbAssembler()
    pa.add_statements(stmts)
    model = pa.make_model()
    # Add observables
    o = Observable('MAPK1p', model.monomers['MAPK1'](T185='p', Y187='p'))
    model.add_component(o)
    o = Observable('MAPK3p', model.monomers['MAPK3'](T202='p', Y204='p'))
    model.add_component(o)
    o = Observable('GSK3Ap', model.monomers['GSK3A'](S21='p'))
    model.add_component(o)
    o = Observable('GSK3Bp', model.monomers['GSK3B'](S9='p'))
    model.add_component(o)
    o = Observable('RPS6p', model.monomers['RPS6'](S235='p'))
    model.add_component(o)
    o = Observable('EIF4EBP1p', model.monomers['EIF4EBP1'](phospho='p'))
    model.add_component(o)
    o = Observable('JUNp', model.monomers['JUN'](S73='p'))
    model.add_component(o)
    o = Observable('FOXO3p', model.monomers['FOXO3'](S315='p'))
    model.add_component(o)
    o = Observable('AKT1p', model.monomers['AKT1'](S473='p'))
    model.add_component(o)
    o = Observable('AKT2p', model.monomers['AKT2'](S474='p'))
    model.add_component(o)
    o = Observable('AKT3p', model.monomers['AKT3'](S='p'))
    model.add_component(o)
    o = Observable('ELK1p', model.monomers['ELK1'](S383='p'))
    model.add_component(o)
    # Set context
    pa.set_context('SKMEL28_SKIN')
    pa.save_model(out_file)
    base_file, _ = os.path.splitext(out_file)
    pa.export_model('kappa', '%s.ka' % base_file)
    return model


def assemble_index_cards(stmts, out_folder):
    counter = 1
    for st in stmts:
        if isinstance(st, Modification) and st.enz is None:
            continue
        pmids = [ev.pmid for ev in st.evidence if ev.pmid is not None]
        if pmids:
            pmids = ','.join(['PMID%s' % pm for pm in list(set(pmids))])
        else:
            pmids = 'N/A'
        ica = IndexCardAssembler([st], pmc_override=pmids)
        ica.make_model()
        if ica.cards:
            ica.save_model(pjoin(out_folder, 'index_card_%d.json' % counter))
            counter += 1

# These filters are currently in a very basic implementation stage
# where they work well on the statements encountered in the scope of this
# evaluation. They could be generalized and refactored to become
# part of the core of INDRA.
# Note that implicitly these filters also filter out statements in which
# the subject is None.
def filter_enzyme_kinase(stmts_in):
    kinase_table = read_unicode_csv('../../indra/resources/kinases.tsv',
                                    delimiter='\t')
    gene_names = [lin[1] for lin in list(kinase_table)[1:]]
    stmts_out = []
    for st in stmts_in:
        if isinstance(st, Phosphorylation):
            if st.enz is not None:
                if st.enz.name in gene_names:
                    stmts_out.append(st)
        else:
            stmts_out.append(st)
    return stmts_out

def filter_mod_nokinase(stmts_in):
    kinase_table = read_unicode_csv('../../indra/resources/kinases.tsv',
                                    delimiter='\t')
    gene_names = [lin[1] for lin in list(kinase_table)[1:]]
    stmts_out = []
    for st in stmts_in:
        if isinstance(st, Modification) and not \
           isinstance(st, Phosphorylation):
            if st.enz is not None:
                if st.enz.name not in gene_names:
                    stmts_out.append(st)
        else:
            stmts_out.append(st)
    return stmts_out

def filter_transcription_factor(stmts_in):
    tf_table = \
        read_unicode_csv('../../indra/resources/transcription_factors.csv')
    gene_names = [lin[1] for lin in list(tf_table)[1:]]
    stmts_out = []
    for st in stmts_in:
        if isinstance(st, RegulateAmount):
            if st.subj is not None:
                if st.subj.name in gene_names:
                    stmts_out.append(st)
        else:
            stmts_out.append(st)
    return stmts_out
