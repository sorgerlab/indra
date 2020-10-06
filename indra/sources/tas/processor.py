__all__ = ['TasProcessor']

import logging
from indra.statements import Inhibition, Agent, Evidence
from indra.statements.validate import assert_valid_db_refs
from indra.ontology.standardize import standardize_name_db_refs
from indra.databases import hgnc_client, chembl_client, lincs_client


logger = logging.getLogger(__name__)


CLASS_MAP = {'1': 'Kd < 100nM', '2': '100nM < Kd < 1uM',
             '3': '1uM < Kd < 10uM', '10': 'Kd > 10uM'}


lincs_client_obj = lincs_client.LincsClient()


class TasProcessor(object):
    """A processor for the Target Affinity Spectrum data table."""
    def __init__(self, data, affinity_class_limit=2, named_only=False,
                 standardized_only=False):
        self._data = data
        self.affinity_class_limit = affinity_class_limit
        self.named_only = named_only
        self.standardized_only = standardized_only

        self.statements = []
        for row in data:
            # Skip rows that are above the affinity class limit
            if int(row['tas']) > affinity_class_limit:
                continue
            self._process_row(row)
        return

    def _process_row(self, row):
        drugs = self._extract_drugs(row['compound_ids'], row['lspci_id'])
        prot = self._extract_protein(row['entrez_gene_symbol'],
                                     row['entrez_gene_id'])
        evidences = self._make_evidences(row['tas'], row['references'])
        # NOTE: there are several entries in this data set that refer to
        # non-human Entrez genes, e.g.
        # https://www.ncbi.nlm.nih.gov/gene/3283880
        # We skip these for now because resources for Entrez-based
        # mappings for non-human genes are not integrated, and would cause
        # pre-assembly issues.
        if 'HGNC' not in prot.db_refs:
            return
        for drug in drugs:
            self.statements.append(Inhibition(drug, prot, evidence=evidences))

    def _extract_drugs(self, compound_ids, lspci_id):
        drugs = []
        for id_ in compound_ids.split('|'):
            db_refs = {'LSPCI': lspci_id}
            if id_.startswith('CHEMBL'):
                db_refs['CHEMBL'] = id_
            elif id_.startswith('HMSL'):
                db_refs['HMS-LINCS'] = id_.split('HMSL')[1]
            else:
                logger.warning('Unhandled ID type: %s' % id_)
            # Name standardization finds correct names but because
            # ChEMBL is incomplete as a local resource, we don't
            # universally standardize its names, instead, we look
            # it up explicitly when necessary.
            name, db_refs = standardize_name_db_refs(db_refs)
            if name is None:
                # This is one way to detect that the drug could not be
                # standardized beyond just its name so in the
                # standardized_only condition, we skip this drug
                if self.standardized_only:
                    continue
                elif 'HMS-LINCS' in db_refs:
                    name = \
                        lincs_client_obj.get_small_molecule_name(
                            db_refs['HMS-LINCS'])
                elif 'CHEMBL' in db_refs:
                    name = chembl_client.get_chembl_name(db_refs['CHEMBL'])
            # If name is still None, we just use the ID as the name
            if name is None:
                # With the named_only restriction, we skip drugs without
                # a proper name.
                if self.named_only:
                    continue
                name = id_
            assert_valid_db_refs(db_refs)
            drugs.append(Agent(name, db_refs=db_refs))
        drugs = list({agent.matches_key():
                      agent for agent in drugs}.values())
        return drugs

    def _extract_protein(self, name, gene_id):
        refs = {'EGID': gene_id}
        hgnc_id = hgnc_client.get_hgnc_from_entrez(gene_id)
        if hgnc_id is not None:
            refs['HGNC'] = hgnc_id
        standard_name, db_refs = standardize_name_db_refs(refs)
        if standard_name:
            name = standard_name
        assert_valid_db_refs(db_refs)
        return Agent(name, db_refs=db_refs)

    def _make_evidences(self, class_min, references):
        evidences = []
        for reference in references.split('|'):
            pmid, source_id, text_refs = None, None, None
            annotations = {'class_min': CLASS_MAP[class_min]}
            ref, id_ = reference.split(':')
            if ref == 'pubmed':
                pmid = id_
                text_refs = {'PMID': pmid}
            elif ref == 'doi':
                text_refs = {'DOI': id_}
            else:
                source_id = reference
            ev = Evidence(source_api='tas', source_id=source_id, pmid=pmid,
                          annotations=annotations, epistemics={'direct': True},
                          text_refs=text_refs)
            evidences.append(ev)
        return evidences
