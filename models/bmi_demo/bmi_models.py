import os
import shutil
import matplotlib.pyplot as plt
from indra.sources import eidos
from indra.assemblers.bmi_wrapper import BMIModel
from indra.assemblers import PysbAssembler
from topoflow.framework import emeli

def text_to_stmts(text):
    """Run Eidos reading on a given text and return INDRA Statements."""
    fname = text.replace(' ', '_').replace(',', '_') + '.jsonld'
    if os.path.exists(fname):
        ep = eidos.process_json_ld_file(fname)
    else:
        ep = eidos.process_text(text)
        shutil.move('eidos_output.json', fname)
    return ep.statements


def make_component_repo(bmi_models):
    """Generate files representing each component for EMELI to discover."""
    for m in bmi_models:
        m.export_into_python()
    comps = [m.make_repository_component() for m in bmi_models]
    rep_str = '<repository>%s</repository>' % '\n'.join(comps)
    with open('component_repository.xml', 'w') as fh:
        fh.write(rep_str)
    with open('component_providers.txt', 'w') as fh:
        for m in bmi_models:
            fh.write('%s %s\n' % (m.model.name, m.model.name))


def plot_results(model_dict):
    plt.figure()
    plt.ion()
    for model_name, bmi_model in model_dict.items():
        state_dict = {}
        times = [tc[0] for tc in bmi_model.time_course]
        for var_name, var_id in bmi_model.species_name_map.items():
            state_dict[var_name] = [tc[1][var_id] for tc in
                                    bmi_model.time_course]
            plt.plot(times, state_dict[var_name], label=var_name)
    #plt.xlim([1, 1001])
    plt.legend()
    plt.show()


# Variables of interest in topoflow/components/met_base.py
# _output_var_names
# 'atmosphere_water__rainfall_volume_flux'
# 'land_surface__temperature'
# 'land_surface_net-shortwave-radiation__energy_flux'

if __name__ == '__main__':
    model_txts = ['rainfall causes floods',
        'floods cause displacement, and displacement reduces access to food']
    stmts = [text_to_stmts(t) for t in model_txts]
    bmi_models = []
    root_vars = [['rainfall'], []]
    for idx, model_stmts in enumerate(stmts):
        pa = PysbAssembler()
        pa.add_statements(model_stmts)
        model = pa.make_model()
        model.name = 'model%d' % idx
        bm = BMIModel(model, root_vars=root_vars[idx], stop_time=5000)
        bmi_models.append(bm)
    make_component_repo(bmi_models)

    # Make the EMELI framework and run it
    f = emeli.framework()
    f.run_model(cfg_prefix='component', cfg_directory='.',
                driver_comp_name=bmi_models[0].model.name)
    plot_results(f.comp_set)
