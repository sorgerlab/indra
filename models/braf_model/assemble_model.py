from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
from pysb import *
import pysb.export
from indra import trips
from indra.assemblers import PysbAssembler

def apply_patch(original, patch):
    orig_lines = [s.strip() for s in original.strip().split('\n')]
    patch_lines = [s.strip() for s in patch.strip().split('\n')]
    remove_lines = []
    add_lines = []
    for p in patch_lines:
        if p.startswith('<'):
            remove_lines.append(p[2:].strip())
        elif p.startswith('>'):
            add_lines.append(p[2:].strip())
    new_lines = []
    for o in orig_lines:
        if o not in remove_lines:
            new_lines.append(o)
    new_lines += add_lines
    new_txt = '\n'.join(new_lines)
    return new_txt

def assemble_model(model_id, reread=False):
    model_name = 'model%d' % model_id
    # If model has already been read, just process the EKB XML
    if os.path.exists(model_name + '.xml') and not reread:
        tp = trips.process_xml(open(model_name + '.xml').read())
    else:
        # Start with the basic model
        model_txt = open('model1.txt').read()
        # Apply patches one by one to get to the current model text
        for j in range(1, model_id):
            patch_txt = open('model%d_from%d.txt' % (j+1, j)).read()
            model_txt = apply_patch(model_txt, patch_txt)
        print('Reading model %d text:' % model_id)
        print(model_txt)
        # Process model text and save result EKB XML
        tp = trips.process_text(model_txt, model_name + '.xml')

    print('Assembling statements:')
    for i, st in enumerate(tp.statements):
        print('%d: %s' % (i, st))
    # Assemble the PySB model
    pa = PysbAssembler()
    pa.add_statements(tp.statements)
    model = pa.make_model(policies='two_step')

    # Set initial conditions
    erk = model.monomers['ERK']
    obs = Observable(b'ERK_p', erk(phospho='p'))
    model.add_component(obs)
    vem = model.monomers['VEMURAFENIB']
    obs = Observable(b'Vem_free', vem(map3k=None))
    model.add_component(obs)
    ras = model.monomers['RAS']
    obs = Observable(b'RAS_active', ras(gtp=ANY))
    model.add_component(obs)
    braf = model.monomers['BRAF']
    obs = Observable(b'BRAF_active', braf(vemurafenib=None))
    model.add_component(obs)
    model.parameters[b'BRAF_0'].value = 0
    egf = model.monomers['EGF']
    obs = Observable(b'EGF_free', egf(erbb=None))
    model.add_component(obs)

    # Add mutated form of BRAF as initial condition
    sites_dict = {}
    for site in braf.sites:
        if site in braf.site_states:
            sites_dict[site] = braf.site_states[site][0]
        else:
            sites_dict[site] = None
    sites_dict['V600'] = 'E'
    model.add_component(Parameter('BRAF_mut_0', 1e5))
    model.initial(braf(**sites_dict), model.parameters['BRAF_mut_0'])

    # Set up model parameters
    model.parameters['kf_ee_bind_1'].value = 1
    model.parameters['kr_ee_bind_1'].value = 0.1
    model.parameters['kf_ee_bind_2'].value = 1
    model.parameters['kr_ee_bind_2'].value = 0.1
    model.parameters['kf_eg_bind_1'].value = 1
    model.parameters['kr_eg_bind_1'].value = 0.1
    model.parameters['kf_gs_bind_1'].value = 1
    model.parameters['kr_gs_bind_1'].value = 0.1
    model.parameters['kf_sr_bind_1'].value = 1
    model.parameters['kr_sr_bind_1'].value = 50
    model.parameters['kf_rg_bind_1'].value = 50
    model.parameters['kr_rg_bind_1'].value = 0.5
    model.parameters['kf_rb_bind_1'].value = 1
    model.parameters['kr_rb_bind_1'].value = 0.5

    model.parameters['kf_vb_bind_1'].value = 10
    model.parameters['kr_vb_bind_1'].value = 1

    model.parameters['kf_bm_bind_1'].value = 1
    model.parameters['kr_bm_bind_1'].value = 0.1
    model.parameters['kc_bm_phosphorylation_1'].value = 3
    model.parameters['kf_pm_bind_1'].value = 1
    model.parameters['kr_pm_bind_1'].value = 0.001
    model.parameters['kc_pm_dephosphorylation_1'].value = 10
    model.parameters['kf_me_bind_1'].value = 1
    model.parameters['kr_me_bind_1'].value = 0.1
    model.parameters['kc_me_phosphorylation_1'].value = 10
    model.parameters['kf_de_bind_1'].value = 1
    model.parameters['kr_de_bind_1'].value = 0.001
    model.parameters['kc_de_dephosphorylation_1'].value = 10


    model.parameters['VEMURAFENIB_0'].value = 0
    model.parameters['EGF_0'].value = 1e3
    model.parameters['EGFR_0'].value = 1e5
    model.parameters['SOS_0'].value = 1e3
    model.parameters['GRB2_0'].value = 1e5
    model.parameters['RAS_0'].value = 2e5
    model.parameters['GTP_0'].value = 1e7
    model.parameters['MEK_0'].value = 1e5
    model.parameters['ERK_0'].value = 1e5
    model.parameters['DUSP6_0'].value = 1e3
    model.parameters['PPP2CA_0'].value = 1e5

    if model_id >= 2:
        model.parameters['PHOSPHATASE_0'].value = 1e2
        model.parameters['kf_es_bind_1'].value = 1e-05
        model.parameters['kr_es_bind_1'].value = 1e-04
        model.parameters['kc_es_phosphorylation_1'].value = 1
        model.parameters['kf_ps_bind_1'].value = 1
        model.parameters['kr_ps_bind_1'].value = 0.1
        model.parameters['kc_ps_dephosphorylation_1'].value = 1e-04

    if model_id >= 3:
        model.parameters['kf_bb_bind_1'].value = 10
        model.parameters['kr_bb_bind_1'].value = 1

    if model_id == 4:
        model.parameters['kf_vb_bind_2'].value = 1e-04

    pa.model = model
    pa.save_model('model%d.py' % model_id)
    return model
