import os
import sys
import shutil
import matplotlib.pyplot as plt
from indra.sources import eidos
from indra.statements import DecreaseAmount, Concept
from indra.assemblers.bmi_wrapper import BMIModel
from indra.assemblers import PysbAssembler
from topoflow.framework import emeli
import eval_model

# This is a specific configuration of the Topoflow model we will use
# topoflow_config = '/Users/ben/src/topoflow/topoflow/examples/Treynor_Iowa_30m'
# cfg_prefix = 'June_20_67'
topoflow_config = '/Users/ben/src/topoflow/topoflow/examples/Gel-Aliab'
cfg_prefix = 'Test2'

def text_to_stmts(text):
    """Run Eidos reading on a given text and return INDRA Statements."""
    # We use some caching here so that sentences we have already read
    # are not re-read.
    fname = text.replace(' ', '_').replace(',', '_') + '.jsonld'
    if os.path.exists(fname):
        ep = eidos.process_json_ld_file(fname)
    else:
        ep = eidos.process_text(text)
        shutil.move('eidos_output.json', fname)
    return ep.statements


def make_component_repo(bmi_models, topo):
    """Generate files representing each component for EMELI to discover."""
    # First we make the BMI models self-export into a python file from
    # which they can be loaded
    for m in bmi_models:
        m.export_into_python()
    # We next generate the component repository XML that EMELI requires
    # to discover components in the simulation workflow
    comps = [m.make_repository_component() for m in bmi_models]
    comp_xmls = '\n'.join(comps)
    # If Topoflow components are included in the simulation then an additional
    # XML block describing these components needs to be appended
    if topo:
        with open('met_comp.xml', 'r') as fh:
            met_comp = fh.read()
        comp_xmls += met_comp
    rep_str = '<repository>%s</repository>' % comp_xmls
    # Finally, the repository XML is saved
    with open('component_repository.xml', 'w') as fh:
        fh.write(rep_str)
    # EMELI also requires a file which lists the names of models, in general
    # this file can be used to choose between alternative providers of
    # models with the same variables
    with open('component_providers.txt', 'w') as fh:
        for m in bmi_models:
            fh.write('%s %s\n' % (m.model.name, m.model.name))


def plot_results(model_dict, vars_to_plot=None):
    """Plot the results of a simulation."""
    plt.figure()
    plt.ion()
    for model_name, bmi_model in model_dict.items():
        if model_name not in ['indra_model0', 'indra_model1',
                              'indra_eval_model']:
            continue
        state_dict = {}
        times = [tc[0] for tc in bmi_model.time_course]
        for var_name, var_id in bmi_model.species_name_map.items():
            if not vars_to_plot or var_name in vars_to_plot:
                state_dict[var_name] = [tc[1][var_id] for tc in
                                        bmi_model.time_course]
                short_var_name = var_name[:20]
                plt.plot(times, state_dict[var_name], label=short_var_name)
    plt.legend()
    plt.show()


# Variables of interest in topoflow/components/met_base.py
# _output_var_names
# 'atmosphere_water__rainfall_volume_flux'
# 'land_surface__temperature'
# 'land_surface_net-shortwave-radiation__energy_flux'

if __name__ == '__main__':
    # Choose which demo we are running. A new Python session needs to be
    # started in between demos to clear the modules loaded by EMELI.
    if len(sys.argv) < 2:
        demo_idx = 1
    else:
        demo_idx = int(sys.argv[1])


    if demo_idx in (1, 2):
        # The content of the two independent models described in natural language
        model_txts = ['rainfall causes floods',
            'floods cause displacement, and displacement reduces access to food']
        # Read with Eidos and extract INDRA Statements for each NL model
        stmts = [text_to_stmts(t) for t in model_txts]
        # It makes sense to assume that floods go away naturally and that
        # displacement is also decreased naturally over time. These are not
        # really needed for the demo but are conceptually interesting to think
        # about in terms of information missing from the original descriptions.
        stmts[0].append(DecreaseAmount(None, Concept('flood')))
        stmts[1].append(DecreaseAmount(None, Concept('displacement')))
        # We now create the variable mappings and the assumed "root" variables
        # for each demo case. In Demo 1, "rainfall" is assumed to be a root
        # variable that is assumed to be fixed. In demo 2, "rainfall" is mapped
        # to a corresponding Topoflow variable.
    elif demo_idx == 3:
        stmts = [eval_model.stmts]
    bmi_models = []
    if demo_idx == 1:
        out_name_maps = [{}, {}]
        input_vars = [[], ['flood']]
    elif demo_idx == 2:
        out_name_maps = [{'atmosphere_water__rainfall_volume_flux': 'rainfall'},
                         {}]
        input_vars = [['rainfall'], ['flood']]
    else:
        out_name_maps = [{'atmosphere_water__rainfall_volume_flux': 
                          'Precipitation'}]
        input_vars = [['Precipitation']]
    # We now assemble PySB models from the INDRA Statements and then
    # instantiate these models as BMI-wrapped models along with a simulator
    for idx, model_stmts in enumerate(stmts):
        pa = PysbAssembler()
        pa.add_statements(model_stmts)
        model = pa.make_model()
        if demo_idx in (1, 2):
            model.name = 'indra_model%d' % idx
        else:
            model.name = 'indra_eval_model'
        bm = BMIModel(model, inputs=input_vars[idx], stop_time=50000,
                      outside_name_map=out_name_maps[idx])
        bmi_models.append(bm)

    # Example 1: two NL models co-simulated
    if demo_idx == 1:
        # We make the model component repository without Topoflow
        make_component_repo(bmi_models, False)
        # We instantiate the EMELI framework and then run the simulations
        f = emeli.framework()
        f.run_model(cfg_prefix='component', cfg_directory='.',
                    driver_comp_name=bmi_models[0].model.name)
        # Finally plot the results
        plot_results(f.comp_set)

    # Example 2: two NL models + Topoflow co-simulated
    elif demo_idx == 2:
        # We make the model component repository with Topoflow
        make_component_repo(bmi_models, True)
        f = emeli.framework()
        # We instantiate the EMELI framework and then run the simulations
        f.run_model(cfg_prefix=cfg_prefix, cfg_directory=topoflow_config,
                    driver_comp_name=bmi_models[0].model.name)
        # Finally plot the results
        plot_results(f.comp_set)

    # Example 3: evaluation model + Topoflow co-simulated
    elif demo_idx == 3:
        # We make the model component repository with Topoflow
        make_component_repo(bmi_models, True)
        f = emeli.framework()
        # We instantiate the EMELI framework and then run the simulations
        f.run_model(cfg_prefix=cfg_prefix, cfg_directory=topoflow_config,
                    driver_comp_name=bmi_models[0].model.name)
        # Finally plot the results
        plot_results(f.comp_set, ['Precipitation',
                                  'Crop_technology'])
