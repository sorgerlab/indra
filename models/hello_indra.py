# exported from PySB model 'None'

from pysb import Model, Monomer, Parameter, Expression, Compartment, Rule, Observable, Initial, MatchOnce, Annotation, ANY, WILD

Model()

Monomer(u'MAP2K1')
Monomer(u'MAPK1', [u'T185', u'Y187'], {u'Y187': [u'u', u'p'], u'T185': [u'u', u'p']})

Parameter(u'kf_mm_phosphorylation_1', 1e-06)
Parameter(u'kf_mm_phosphorylation_2', 1e-06)
Parameter(u'MAP2K1_0', 1000.0)
Parameter(u'MAPK1_0', 1000.0)

Rule(u'MAP2K1_phosphorylation_MAPK1_T185', MAP2K1() + MAPK1(T185=u'u') >> MAP2K1() + MAPK1(T185=u'p'), kf_mm_phosphorylation_1)
Rule(u'MAP2K1_phosphorylation_MAPK1_Y187', MAP2K1() + MAPK1(Y187=u'u') >> MAP2K1() + MAPK1(Y187=u'p'), kf_mm_phosphorylation_2)

Initial(MAP2K1(), MAP2K1_0)
Initial(MAPK1(T185=u'u', Y187=u'u'), MAPK1_0)

Annotation(MAP2K1, u'http://identifiers.org/hgnc/HGNC:6840', u'is')
Annotation(MAP2K1, u'http://identifiers.org/ncit/C52823', u'is')
Annotation(MAP2K1, u'http://identifiers.org/uniprot/Q02750', u'is')
Annotation(MAPK1, u'http://identifiers.org/hgnc/HGNC:6871', u'is')
Annotation(MAPK1, u'http://identifiers.org/ncit/C52872', u'is')
Annotation(MAPK1, u'http://identifiers.org/uniprot/P28482', u'is')
Annotation(MAP2K1_phosphorylation_MAPK1_T185, u'MAP2K1', u'rule_has_subject')
Annotation(MAP2K1_phosphorylation_MAPK1_T185, u'MAPK1', u'rule_has_object')
Annotation(MAP2K1_phosphorylation_MAPK1_T185, u'd837945f-70cb-4acb-84fc-c583fd16e86d', u'from_indra_statement')
Annotation(MAP2K1_phosphorylation_MAPK1_Y187, u'MAP2K1', u'rule_has_subject')
Annotation(MAP2K1_phosphorylation_MAPK1_Y187, u'MAPK1', u'rule_has_object')
Annotation(MAP2K1_phosphorylation_MAPK1_Y187, u'f3c8a0b4-db9f-4218-8f96-d00ccb369c76', u'from_indra_statement')

