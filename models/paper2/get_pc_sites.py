import pickle
from indra.statements import Agent, ModCondition
from indra.sources.biopax import processor as bpc
from indra.sources.biopax import pathway_commons_client as pcc

owl_pattern = '/home/bmg16/data/pathwaycommons/PathwayCommons9.%s.BIOPAX.owl'
dbs = ['psp', 'pid', 'reactome', 'kegg', 'panther', 'hprd']

def get_proteins(mod_feature):
    # Assume it's a physical entity feature
    pe = obj.getFeatureOf().toArray()
    return pe

for db in dbs:
    owl_file = owl_pattern % db
    print('Reading %s...' % owl_file)
    model = pcc.owl_to_model(owl_file)
    mf_class = bpc._bpimpl('ModificationFeature')

    objs = model.getObjects().toArray()

    agents = []
    for obj in objs:
        if not isinstance(obj, mf_class):
            continue
        try:
            mc = bpc.BiopaxProcessor._extract_mod_from_feature(obj)
        except Exception as e:
            print('ERROR: ' + str(e))
            continue
        if not mc or not mc.residue or not mc.position:
            continue

        proteins = get_proteins(obj)
        if not proteins:
            continue
        for protein in proteins:
            name = bpc.BiopaxProcessor._get_element_name(protein)
            db_refs = bpc.BiopaxProcessor._get_db_refs(protein)
            agent = Agent(name, mods=[mc], db_refs=db_refs)
            reactions = protein.getParticipantOf().toArray()
            if not reactions:
                upstream = protein.getMemberPhysicalEntityOf().toArray()
                for u in upstream:
                    reactions += u.getParticipantOf().toArray()
            for reaction in reactions:
                controls = reaction.getControlledOf().toArray()
                if not controls:
                    agents.append(agent)
                for contr in controls:
                    agents.append(agent)

    with open('pc_%s_modified_agents.pkl' % db, 'wb') as fh:
        pickle.dump(agents, fh)
