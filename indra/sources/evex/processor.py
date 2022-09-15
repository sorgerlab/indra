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
    def __init__(self, relations_table, articles_table, standoff_index):
        self.relations_table = relations_table
        self.articles_table = articles_table
        self.article_lookup = self.build_article_lookup()
        self.statements = []

    def process_statements(self):
        for row in tqdm.tqdm(self.relations_table.itertuples(),
                             total=len(self.relations_table),
                             desc='Processing Evex relations'):
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

            article_prefix, article_id = \
                self.article_lookup.get(row.general_event_id)

            text_refs = {article_prefix: article_id}
            pmid = article_id if article_prefix == 'PMID' else None

            ev = Evidence(source_api='evex',
                          pmid=pmid,
                          text_refs=text_refs)

            if stmt_type == Complex:
                stmt = Complex([subj_agent, obj_agent], evidence=[ev])
            else:
                stmt = stmt_type(subj_agent, obj_agent, evidence=[ev])
            self.statements.append(stmt)

    def build_article_lookup(self):
        article_lookup = {}
        for row in self.articles_table.itertuples():
            prefix, article_id = row.article_id.split(': ')
            if prefix == 'PMCID':
                if not article_id.startswith('PMC'):
                    article_id = 'PMC' + article_id
                article_lookup[row.general_event_id] = ('PMCID', article_id)
            elif prefix == 'PMID':
                article_lookup[row.general_event_id] = ('PMID', article_id)
            else:
                ValueError('Unexpected article type: %s' % prefix)
        return article_lookup