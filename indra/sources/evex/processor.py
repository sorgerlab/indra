import tqdm
from indra.statements import *
from indra.ontology.standardize import get_standard_agent

type_mappings = {
    'Binding': Complex,
    'Catalysis of DNA methylation': Methylation,
    'Catalysis of acetylation': Acetylation,
    'Catalysis of glycosylation': Glycosylation,
    'Catalysis of hydroxylation': Hydroxylation,
    'Catalysis of methylation': Methylation,
    'Catalysis of phosphorylation': Phosphorylation,
    'Catalysis of ubiquitination': Ubiquitination,
    'Indirect_catalysis of acetylation': Acetylation,
    'Indirect_catalysis of methylation': Methylation,
    'Indirect_catalysis of ubiquitination': Ubiquitination,
    'Indirect_regulation': None,
    'Indirect_regulation of binding': None,
    'Indirect_regulation of catabolism': None,
    'Indirect_regulation of expression': IncreaseAmount,
    'Indirect_regulation of localization': None,
    'Indirect_regulation of phosphorylation': Phosphorylation,
    'Indirect_regulation of transcription': IncreaseAmount,
    'Regulation': None,
    'Regulation of binding': None,
    'Regulation of catabolism': None,
    'Regulation of expression': IncreaseAmount,
    'Regulation of localization': None,
    'Regulation of phosphorylation': Phosphorylation,
    'Regulation of transcription': IncreaseAmount,
}


class EvexProcessor:
    def __init__(self, relations_table, standoff_index):
        self.relations_table = relations_table
        self.statements = []

    def process_statements(self):
        for row in tqdm.tqdm(self.relations_table.itertuples(),
                             total=len(self.relations_table)):
            stmt_type = type_mappings.get(row.refined_type)
            if not stmt_type:
                continue
            subj_agent = get_standard_agent(
                'EGID:%s' % str(row.source_entrezgene_id),
                db_refs={'EGID': str(row.source_entrezgene_id)}
            )
            obj_agent = get_standard_agent(
                'EGID:%s' % str(row.target_entrezgene_id),
                db_refs={'EGID': str(row.target_entrezgene_id)}
            )

            if stmt_type == Complex:
                stmt = Complex([subj_agent, obj_agent])
            else:
                stmt = stmt_type(subj_agent, obj_agent)
            self.statements.append(stmt)