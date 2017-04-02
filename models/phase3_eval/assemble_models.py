from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import sys
import pickle
from os.path import join as pjoin
from indra.tools import assemble_corpus as ac
from indra.tools.gene_network import GeneNetwork

import process_data, process_r3, process_sparser, process_trips
import read_phosphosite
from assemble_sif import assemble_sif
from assemble_cx import assemble_cx
from assemble_pysb import assemble_pysb

def build_prior(genes, out_file):
    gn = GeneNetwork(genes, 'korkut')
    stmts = gn.get_statements(filter=False)
    ac.dump_statements(stmts, out_file)
    return stmts

def read_extra_sources(out_file):
    sparser_stmts = process_sparser.read_stmts(process_sparser.base_folder)
    sparser_stmts += \
        process_sparser.read_stmts(process_sparser.sentences_folder)
    r3_stmts = process_r3.read_stmts(process_r3.active_forms_files[0])
    trips_stmts = process_trips.read_stmts(process_trips.base_folder)
    phosphosite_stmts, _ = \
        read_phosphosite.read_phosphosite(read_phosphosite.phosphosite_file)
    stmts = trips_stmts + sparser_stmts + r3_stmts + phosphosite_stmts
    ac.dump_statements(stmts, out_file)
    return stmts

def get_prior_genes(fname):
    """Get the list of prior genes."""
    with open(fname, 'rt') as fh:
        genes = fh.read().strip().split('\n')
        return genes

if __name__ == '__main__':
    if len(sys.argv) < 2:
        assemble_models = ['pysb', 'sif', 'cx']
    else:
        model_types = sys.argv[1:]
        if 'all' in model_types:
            assemble_models = ['pysb', 'sif', 'cx']
        else:
            assemble_models = sys.argv[1:]

    print('Assembling the following model types: %s' % \
          ', '.join(assemble_models))
    print('##############')

    outf = 'output/'
    data = process_data.read_data(process_data.data_file)
    data_genes = process_data.get_all_gene_names(data)
    reassemble = False
    if not reassemble:
        stmts = ac.load_statements(pjoin(outf, 'preassembled.pkl'))
    else:
        #prior_stmts = build_prior(data_genes, pjoin(outf, 'prior.pkl'))
        prior_stmts = ac.load_statements(pjoin(outf, 'prior.pkl'))
        prior_stmts = ac.map_grounding(prior_stmts,
                                       save=pjoin(outf, 'gmapped_prior.pkl'))
        reach_stmts = ac.load_statements(pjoin(outf, 'phase3_stmts.pkl'))
        reach_stmts = ac.filter_no_hypothesis(reach_stmts)
        #extra_stmts = ac.load_statements(pjoin(outf, 'extra_stmts.pkl'))
        extra_stmts = read_extra_sources(pjoin(outf, 'extra_stmts.pkl'))
        reading_stmts = reach_stmts + extra_stmts
        reading_stmts = ac.map_grounding(reading_stmts,
                                         save=pjoin(outf, 'gmapped_reading.pkl'))
        stmts = prior_stmts + reading_stmts + extra_stmts

        stmts = ac.filter_grounded_only(stmts)
        stmts = ac.filter_genes_only(stmts, specific_only=False)
        stmts = ac.filter_human_only(stmts)
        stmts = ac.expand_families(stmts)
        stmts = ac.filter_gene_list(stmts, data_genes, 'one')
        stmts = ac.map_sequence(stmts, save=pjoin(outf, 'smapped.pkl'))
        stmts = ac.run_preassembly(stmts, return_toplevel=False,
                                   save=pjoin(outf, 'preassembled.pkl'))

    ### PySB assembly
    if 'pysb' in assemble_models:
        pysb_model = assemble_pysb(stmts, data_genes,
                                   pjoin(outf, 'korkut_model_pysb.py'))
    ### SIF assembly
    if 'sif' in assemble_models:
        sif_str = assemble_sif(stmts, data, pjoin(outf, 'PKN-korkut_all_ab.sif'))
    ### CX assembly
    if 'cx' in assemble_models:
        cxa = assemble_cx(stmts, pjoin(outf, 'korkut_full_high_belief.cx'))
