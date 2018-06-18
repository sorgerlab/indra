import os
import shutil
from indra.sources import eidos
from indra.assemblers.bmi_wrapper import BMIModel
from indra.assemblers import PysbAssembler


def text_to_stmts(text):
    fname = text.replace(' ', '_') + '.jsonld'
    if os.path.exists(fname):
        ep = eidos.process_json_ld_file(fname)
    else:
        ep = eidos.process_text(text)
        shutil.move('eidos_output.json', fname)
    return ep.statements


def make_component_repo(bmi_models):
    for m in bmi_models:
        m.export_into_python()
    comps = [m.make_repository_component() for m in bmi_models]
    rep_str = '<repository>%s</repository>' % '\n'.join(comps)
    with open('component_repository.xml', 'w') as fh:
        fh.write(rep_str)
    with open('component_providers.txt', 'w') as fh:
        fh.write('model0 model0\nmodel1 model1')



if __name__ == '__main__':
    model_txts = ['rainfall causes floods', 'floods cause displacement']
    stmts = [text_to_stmts(t) for t in model_txts]
    bmi_models = []
    for idx, model_stmts in enumerate(stmts):
        pa = PysbAssembler()
        pa.add_statements(model_stmts)
        model = pa.make_model()
        model.name = 'model%d' % idx
        bm = BMIModel(model)
        bmi_models.append(bm)
    make_component_repo(bmi_models)
