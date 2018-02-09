# exported from PySB model 'None'

from pysb import Model, Monomer, Parameter, Expression, Compartment, Rule, Observable, Initial, MatchOnce, Annotation, ANY, WILD

Model()

Monomer(u'MDM2')
Monomer(u'TP53')

Parameter(u'kf_tm_synth_1', 4.0)
Parameter(u'Ka_tm_synth_1', 10000.0)
Parameter(u'n_tm_synth_1', 1.0)
Parameter(u'MDM2_0', 10000.0)
Parameter(u'TP53_0', 10000.0)

Observable(u'TP53_synthesizes_MDM2_subj_obs', TP53())

Expression(u'TP53_synthesizes_MDM2_rate', kf_tm_synth_1*TP53_synthesizes_MDM2_subj_obs**(-1 + n_tm_synth_1)*(TP53_synthesizes_MDM2_subj_obs**n_tm_synth_1 + Ka_tm_synth_1**n_tm_synth_1)**(-1))

Rule(u'TP53_synthesizes_MDM2', TP53() >> TP53() + MDM2(), TP53_synthesizes_MDM2_rate)

Initial(MDM2(), MDM2_0)
Initial(TP53(), TP53_0)

Annotation(TP53_synthesizes_MDM2, u'd264fe48-8ee8-49ba-95b1-f34485116b3c', u'from_indra_statement')
Annotation(TP53_synthesizes_MDM2, u'MDM2', u'rule_has_object')
Annotation(TP53_synthesizes_MDM2, u'TP53', u'rule_has_subject')

