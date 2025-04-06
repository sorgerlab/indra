"""
Statements represent mechanistic relationships between biological agents.

Statement classes follow an inheritance hierarchy, with all Statement types
inheriting from the parent class :py:class:`Statement`. At
the next level in the hierarchy are the following classes:

Open Domain

- :py:class:`Event`
- :py:class:`Influence`
- :py:class:`Association`

Biological Domain

- :py:class:`Complex`
- :py:class:`Modification`
- :py:class:`SelfModification`
- :py:class:`RegulateActivity`
- :py:class:`RegulateAmount`
- :py:class:`ActiveForm`
- :py:class:`Translocation`
- :py:class:`Gef`
- :py:class:`Gap`
- :py:class:`Conversion`

There are several types of Statements representing post-translational
modifications that further inherit from
:py:class:`Modification`:

- :py:class:`Phosphorylation`
- :py:class:`Dephosphorylation`
- :py:class:`Ubiquitination`
- :py:class:`Deubiquitination`
- :py:class:`Sumoylation`
- :py:class:`Desumoylation`
- :py:class:`Hydroxylation`
- :py:class:`Dehydroxylation`
- :py:class:`Acetylation`
- :py:class:`Deacetylation`
- :py:class:`Glycosylation`
- :py:class:`Deglycosylation`
- :py:class:`Farnesylation`
- :py:class:`Defarnesylation`
- :py:class:`Geranylgeranylation`
- :py:class:`Degeranylgeranylation`
- :py:class:`Palmitoylation`
- :py:class:`Depalmitoylation`
- :py:class:`Myristoylation`
- :py:class:`Demyristoylation`
- :py:class:`Ribosylation`
- :py:class:`Deribosylation`
- :py:class:`Methylation`
- :py:class:`Demethylation`
- :py:class:`Carboxylation`
- :py:class:`Decarboxylation`
- :py:class:`Amidation`
- :py:class:`Deamidation`
- :py:class:`Sulfation`
- :py:class:`Desulfation`
- :py:class:`Neddylation`
- :py:class:`Deneddylation`
- :py:class:`Alkylation`
- :py:class:`Dealkylation`
- :py:class:`Prenylation`
- :py:class:`Deprenylation`

There are additional subtypes of :py:class:`SelfModification`:

- :py:class:`Autophosphorylation`
- :py:class:`Transphosphorylation`

Interactions between proteins are often described simply in terms of their
effect on a protein's "activity", e.g., "Active MEK activates ERK", or "DUSP6
inactives ERK".  These types of relationships are indicated by the
:py:class:`RegulateActivity` abstract base class which has subtypes

- :py:class:`Activation`
- :py:class:`Inhibition`

while the :py:class:`RegulateAmount` abstract base class has subtypes

- :py:class:`IncreaseAmount`
- :py:class:`DecreaseAmount`

Statements involve one or more *Concepts*, which, depending on the
semantics of the Statement, are typically biological *Agents*,
such as proteins, represented by the class :py:class:`Agent`. (However,
:py:class`Influence` statements involve two or more :py:class`Event` objects,
each of which takes a :py:class`Concept` as an argument.)

Agents can have several types of context specified on them including

- a specific post-translational modification state (indicated by one or
  more instances of :py:class:`ModCondition`),
- other bound Agents (:py:class:`BoundCondition`),
- mutations (:py:class:`MutCondition`),
- an activity state (:py:class:`ActivityCondition`), and
- cellular location

The *active* form of an agent (in terms of its post-translational modifications
or bound state) is indicated by an instance of the class
:py:class:`ActiveForm`.

Grounding and DB references
---------------------------
Agents also carry grounding information which links them to database entries.
These database references are represented as a dictionary in the `db_refs`
attribute of each Agent. The dictionary can have multiple entries. For
instance, INDRA's input Processors produce genes and proteins that carry both
UniProt and HGNC IDs in db_refs, whenever possible. FamPlex provides a name
space for protein families that are typically used in the literature.  More
information about FamPlex can be found here:
https://github.com/sorgerlab/famplex

In general, the capitalized version of any identifiers.org name space (see
https://registry.identifiers.org/ for full list) can be used in db_refs with
a few cases where INDRA's internal db_refs name space is different from the
identifiers.org name space (e.g., UP vs uniprot). These special cases can
be programmatically mapped between INDRA and identifiers.org using the
`identifiers_mappings` and `identifiers_reverse` dictionaries in the
`indra.databases.identifiers` module.

Examples of the most commonly encountered db_refs name spaces and IDs are
listed below.

+------------------------+------------------+------------------------------+
| Type                   | Database         | Example                      |
+========================+==================+==============================+
| Gene/Protein           | HGNC             | {'HGNC': '11998'}            |
+------------------------+------------------+------------------------------+
| Gene/Protein           | UniProt          | {'UP': 'P04637'}             |
+------------------------+------------------+------------------------------+
| Protein chain          | UniProt          | {'UPPRO': 'PRO_0000435839'}  |
+------------------------+------------------+------------------------------+
| Gene/Protein           | Entrez           | {'EGID': '5583'}             |
+------------------------+------------------+------------------------------+
| Gene/Protein family    | FamPlex          | {'FPLX': 'ERK'}              |
+------------------------+------------------+------------------------------+
| Gene/Protein family    | InterPro         | {'IP': 'IPR000308'}          |
+------------------------+------------------+------------------------------+
| Gene/Protein family    | Pfam             | {'PF': 'PF00071'}            |
+------------------------+------------------+------------------------------+
| Gene/Protein family    | NextProt family  | {'NXPFA': '03114'}           |
+------------------------+------------------+------------------------------+
| Chemical               | ChEBI            | {'CHEBI': 'CHEBI:63637'}     |
+------------------------+------------------+------------------------------+
| Chemical               | PubChem          | {'PUBCHEM': '42611257'}      |
+------------------------+------------------+------------------------------+
| Chemical               | LINCS            | {'LINCS': '42611257'}        |
+------------------------+------------------+------------------------------+
| Metabolite             | HMDB             | {'HMDB': 'HMDB00122'}        |
+------------------------+------------------+------------------------------+
| Process, location, etc.| GO               | {'GO': 'GO:0006915'}         |
+------------------------+------------------+------------------------------+
| Process, disease, etc. | MeSH             | {'MESH': 'D008113'}          |
+------------------------+------------------+------------------------------+
| Disease                | Disease Ontology | {'DOID': 'DOID:8659'}        |
+------------------------+------------------+------------------------------+
| Phenotypic abnormality | Human Pheno. Ont.| {'HP': 'HP:0031296'}         |
+------------------------+------------------+------------------------------+
| Experimental factors   | Exp. Factor Ont. | {'EFO': '0007820'}           |
+------------------------+------------------+------------------------------+
| General terms          | NCIT             | {'NCIT': 'C28597'}           |
+------------------------+------------------+------------------------------+
| Raw text               | TEXT             | {'TEXT': 'Nf-kappaB'}        |
+------------------------+------------------+------------------------------+


The evidence for a given Statement, which could include relevant citations,
database identifiers, and passages of text from the scientific literature, is
contained in one or more :py:class:`Evidence` objects associated with the
Statement.


JSON serialization of INDRA Statements
--------------------------------------

Statements can be serialized into JSON and deserialized from JSON to allow
their exchange in a platform-independent way. We also provide a JSON
schema (see http://json-schema.org to learn about schemas) in
https://raw.githubusercontent.com/sorgerlab/indra/master/indra/resources/statements_schema.json
which can be used to validate INDRA Statements JSONs.

Some validation tools include:

- jsonschema
    a Python package to validate JSON content with respect to
    a schema
- ajv-cli
    Available at https://www.npmjs.com/package/ajv-cli
    Install with "npm install -g ajv-cli" and then validate with:
    ajv -s statements_schema.json -d file_to_validate.json. This tool
    provides more sophisticated and better interpretable output than
    jsonschema.
- Web based tools
    There are a variety of web-based tools for validation with JSON schemas,
    including https://www.jsonschemavalidator.net
"""

from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from future.utils import python_2_unicode_compatible

__all__ = [
    # Condition classes
    'BoundCondition', 'MutCondition', 'ModCondition', 'ActivityCondition',

    # Statement classes
    'Statement', 'Modification', 'AddModification', 'RemoveModification',
    'SelfModification', 'Phosphorylation', 'Autophosphorylation',
    'Transphosphorylation', 'Dephosphorylation', 'Hydroxylation',
    'Dehydroxylation', 'Sumoylation', 'Desumoylation', 'Acetylation',
    'Deacetylation', 'Glycosylation', 'Deglycosylation', 'Ribosylation',
    'Deribosylation', 'Ubiquitination', 'Deubiquitination', 'Farnesylation',
    'Defarnesylation', 'Geranylgeranylation', 'Degeranylgeranylation',
    'Palmitoylation', 'Depalmitoylation', 'Myristoylation', 'Demyristoylation',
    'Methylation', 'Demethylation', 'RegulateActivity', 'Inhibition',
    'Activation', 'GtpActivation', 'ActiveForm', 'HasActivity', 'Gef', 'Gap',
    'Complex', 'Translocation', 'RegulateAmount', 'DecreaseAmount',
    'IncreaseAmount', 'Influence', 'Conversion', 'Unresolved',
    'Association', 'Event', 'Migration',
    'Carboxylation', 'Decarboxylation', 'Amidation', 'Deamidation', 'Sulfation',
    'Desulfation', 'Neddylation', 'Deneddylation', 'Alkylation', 'Dealkylation',
    'Prenylation', 'Deprenylation',

    # Error classes
    'InputError', 'UnresolvedUuidError', 'InvalidLocationError',
    'InvalidResidueError', 'NotAStatementName',

    # Other classes
    'Concept', 'Agent', 'Evidence', 'QualitativeDelta', 'QuantitativeState',

    # Context classes
    'BioContext', 'WorldContext', 'TimeContext', 'RefContext', 'Context',
    'MovementContext',

    # Functions and values
    'stmts_from_json', 'get_unresolved_support_uuids', 'stmts_to_json',
    'stmts_from_json_file', 'stmts_to_json_file', 'stmt_from_json',
    'stmt_from_json_str', 'get_valid_residue',
    'draw_stmt_graph', 'get_all_descendants','make_statement_camel',
    'amino_acids', 'amino_acids_reverse', 'activity_types',
    'modtype_to_modclass',
    'modclass_to_modtype', 'modtype_conditions', 'modtype_to_inverse',
    'modclass_to_inverse', 'get_statement_by_name', 'make_hash', 'stmt_type',
    'default_ns_order', 'mk_str', 'pretty_print_stmts', 'print_stmt_summary',
    'set_pretty_print_max_width'
    ]

import abc
import sys
import uuid
import logging
import networkx
import itertools
from copy import deepcopy
from collections import Counter, OrderedDict as _o
from .util import *
from .concept import *
from .context import *
from .evidence import *
from .resources import *
from .delta import *

logger = logging.getLogger(__name__)


try:  # Python 2
    basestring
except NameError:  # Python 3
    basestring = str


class Statement(object):
    """The parent class of all statements.

    Parameters
    ----------
    evidence : None or :py:class:`Evidence` or list of :py:class:`Evidence`
        If a list of Evidence objects is passed to the constructor, the
        value is set to this list. If a bare Evidence object is passed,
        it is enclosed in a list. If no evidence is passed (the default),
        the value is set to an empty list.
    supports : list of :py:class:`Statement`
        Statements that this Statement supports.
    supported_by : list of :py:class:`Statement`
        Statements supported by this statement.
    """

    _agent_order = NotImplemented

    def __init__(self, evidence=None, supports=None, supported_by=None):
        if evidence is None:
            self.evidence = []
        elif isinstance(evidence, Evidence):
            self.evidence = [evidence]
        elif isinstance(evidence, list):
            self.evidence = evidence
        else:
            raise ValueError('evidence must be an Evidence object, a list '
                             '(of Evidence objects), or None.')

        # Initialize supports/supported_by fields, which should be lists
        self.supports = supports if supports else []
        self.supported_by = supported_by if supported_by else []
        self.belief = 1
        self.uuid = '%s' % uuid.uuid4()
        self._full_hash = None
        self._shallow_hash = None
        return

    def matches_key(self):
        raise NotImplementedError("Method must be implemented in child class.")

    def matches(self, other):
        return self.matches_key() == other.matches_key()

    def get_hash(self, shallow=True, refresh=False, matches_fun=None):
        """Get a hash for this Statement.

        There are two types of hash, "shallow" and "full". A shallow hash is
        as unique as the information carried by the statement, i.e. it is a hash
        of the `matches_key`. This means that differences in source, evidence,
        and so on are not included. As such, it is a shorter hash (14 nibbles).
        The odds of a collision among all the statements we expect to encounter
        (well under 10^8) is ~10^-9 (1 in a billion). Checks for collisions can
        be done by using the matches keys.

        A full hash includes, in addition to the matches key, information from
        the evidence of the statement. These hashes will be equal if the two
        Statements came from the same sentences, extracted by the same reader,
        from the same source. These hashes are correspondingly longer (16
        nibbles). The odds of a collision for an expected less than 10^10
        extractions is ~10^-9 (1 in a billion).

        Note that a hash of the Python object will also include the `uuid`, so
        it will always be unique for every object.

        Parameters
        ----------
        shallow : bool
            Choose between the shallow and full hashes described above. Default
            is true (e.g. a shallow hash).
        refresh : bool
            Used to get a new copy of the hash. Default is false, so the hash,
            if it has been already created, will be read from the attribute.
            This is primarily used for speed testing.
        matches_fun : Optional[function]
            A function which takes a Statement as argument and returns a string
            matches key which is then hashed. If not provided the Statement's
            built-in matches_key method is used.

        Returns
        -------
        hash : int
            A long integer hash.
        """
        if shallow:
            # If the hash attribute is missing or we want to refresh it, or
            # we are passing in a custom matches function (since we don't know
            # if the cached value is consistent with it), then we
            # recalculate the hash. Otherwise we just return the cached value.
            if not hasattr(self, '_shallow_hash') or self._shallow_hash is None \
                    or refresh or matches_fun:
                matches_key = matches_fun(self) if matches_fun else \
                    self.matches_key()
                self._shallow_hash = make_hash(matches_key, 14)
            ret = self._shallow_hash
        else:
            # If the hash attribute is missing or we want to refresh it, or
            # we are passing in a custom matches function (since we don't know
            # if the cached value is consistent with it), then we
            # recalculate the hash. Otherwise we just return the cached value.
            if not hasattr(self, '_full_hash') or self._full_hash is None \
                    or refresh or matches_fun:
                ev_mk_list = sorted([ev.matches_key() for ev in self.evidence])
                matches_key = matches_fun(self) if matches_fun else \
                    self.matches_key()
                self._full_hash = \
                    make_hash(matches_key + str(ev_mk_list), 16)
            ret = self._full_hash
        return ret

    def _tag_evidence(self):
        """Set all the Evidence stmt_tag to my deep matches-key hash."""
        h = self.get_hash(shallow=False)
        for ev in self.evidence:
            ev.stmt_tag = h
        return

    def agent_list_with_bound_condition_agents(self):
        # Returns the list of agents both directly participating in the
        # statement and referenced through bound conditions.
        l = self.agent_list()
        for a in self.agent_list():
            if a is not None:
                bc_agents = [bc.agent for bc in a.bound_conditions]
                l.extend(bc_agents)
        return l

    def agent_list(self, deep_sorted=False):
        """Get the canonicalized agent list."""
        ag_list = []
        for ag_name in self._agent_order:
            ag_attr = getattr(self, ag_name)
            if isinstance(ag_attr, Concept) or ag_attr is None:
                ag_list.append(ag_attr)
            elif isinstance(ag_attr, list):
                if not all([isinstance(ag, Concept) for ag in ag_attr]):
                    raise TypeError("Expected all elements of list to be Agent "
                                    "and/or Concept, but got: %s"
                                    % {type(ag) for ag in ag_attr})
                if deep_sorted:
                    ag_attr = sorted_agents(ag_attr)
                ag_list.extend(ag_attr)
            else:
                raise TypeError("Expected type Agent, Concept, or list, got "
                                "type %s." % type(ag_attr))
        return ag_list

    def real_agent_list(self):
        """Return all agents in the statement that are not None."""
        return [a for a in self.agent_list() if a is not None]

    def entities_match(self, other):
        self_key = self.entities_match_key()
        other_key = other.entities_match_key()
        if len(self_key) != len(other_key):
            return False
        for self_agent, other_agent in zip(self_key, other_key):
            if self_agent is None or other_agent is None:
                continue
            if self_agent != other_agent:
                return False
        return True

    def entities_match_key(self):
        key = tuple(a.entity_matches_key() if a is not None
                    else None for a in self.agent_list())
        return key

    def print_supports(self):
        print('%s supported_by:' % str(self))
        if self.supported_by:
            print('-->')
            for s in self.supported_by:
                s.print_supports()

    def __repr__(self):
        if sys.version_info[0] >= 3:
            return str(self)
        else:
            return str(self).encode('utf-8')

    def equals(self, other):
        return self.types_equals(other) \
            and self.agents_equal(other) \
            and self.evidence_equals(other)

    def types_equals(self, other):
        if stmt_type(self) != stmt_type(other):
            return False
        return True

    def agents_equal(self, other):
        if len(self.agent_list()) == len(other.agent_list()):
            for s, o in zip(self.agent_list(), other.agent_list()):
                if (s is None and o is not None) or \
                        (s is not None and o is None):
                    return False
                if s is not None and o is not None and not s.equals(o):
                    return False
        else:
            return False
        return True

    def evidence_equals(self, other):
        if len(self.evidence) == len(other.evidence):
            for s, o in zip(self.evidence, other.evidence):
                if not s.equals(o):
                    return False
        else:
            return False
        return True

    def contradicts(self, other, ontology):
        # Placeholder for implementation in subclasses
        return False

    def to_json(self, use_sbo=False, matches_fun=None):
        """Return serialized Statement as a JSON dict.

        Parameters
        ----------
        use_sbo : Optional[bool]
            If True, SBO annotations are added to each applicable element of
            the JSON. Default: False
        matches_fun : Optional[function]
            A custom function which, if provided, is used to construct the
            matches key which is then hashed and put into the return value.
            Default: None

        Returns
        -------
        json_dict : dict
            The JSON-serialized INDRA Statement.
        """
        stmt_type = type(self).__name__
        # Original comment: For backwards compatibility, could be removed later
        all_stmts = [self] + self.supports + self.supported_by
        for st in all_stmts:
            if not hasattr(st, 'uuid'):
                st.uuid = '%s' % uuid.uuid4()
        ##################
        json_dict = _o(type=stmt_type)
        json_dict['belief'] = self.belief
        if self.evidence:
            evidence = [ev.to_json() for ev in self.evidence]
            json_dict['evidence'] = evidence
        json_dict['id'] = '%s' % self.uuid
        json_dict['matches_hash'] = \
            '%s' % self.get_hash(shallow=True, refresh=True,
                                 matches_fun=matches_fun)
        if self.supports:
            json_dict['supports'] = \
                ['%s' % st.uuid for st in self.supports]
        if self.supported_by:
            json_dict['supported_by'] = \
                ['%s' % st.uuid for st in self.supported_by]

        def get_sbo_term(cls):
            sbo_term = stmt_sbo_map.get(cls.__name__.lower())
            while not sbo_term:
                cls = cls.__bases__[0]
                sbo_term = stmt_sbo_map.get(cls.__name__.lower())
            return sbo_term

        if use_sbo:
            sbo_term = get_sbo_term(self.__class__)
            json_dict['sbo'] = \
                'https://identifiers.org/SBO:%s' % sbo_term
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        stmt_type = json_dict.get('type')
        stmt_cls = getattr(sys.modules[__name__], stmt_type)
        stmt = stmt_cls._from_json(json_dict)
        evidence = json_dict.get('evidence', [])
        stmt.evidence = [Evidence._from_json(ev) for ev in evidence]
        stmt.supports = json_dict.get('supports', [])[:]
        stmt.supported_by = json_dict.get('supported_by', [])[:]
        stmt.belief = json_dict.get('belief', 1.0)
        stmt_id = json_dict.get('id')
        if not stmt_id:
            stmt_id = '%s' % uuid.uuid4()
        stmt.uuid = stmt_id
        json_matches_hash = json_dict.get('matches_hash')
        if json_matches_hash:
            # This code is so central that we need to handle the corner case
            # that a non-int string hash is there without erroring
            try:
                stmt._shallow_hash = int(json_matches_hash)
            except ValueError:
                pass
        return stmt

    def to_graph(self):
        """Return Statement as a networkx graph."""
        def json_node(graph, element, prefix):
            if not element:
                return None
            node_id = '|'.join(prefix)
            if isinstance(element, list):
                graph.add_node(node_id, label='')
                # Enumerate children and add nodes and connect to anchor node
                for i, sub_element in enumerate(element):
                    sub_id = json_node(graph, sub_element, prefix + ['%s' % i])
                    if sub_id:
                        graph.add_edge(node_id, sub_id, label='')
            elif isinstance(element, dict):
                graph.add_node(node_id, label='')
                # Add node recursively for each element
                # Connect to this node with edge label according to key
                for k, v in element.items():
                    if k == 'id':
                        continue
                    elif k == 'name':
                        graph.nodes[node_id]['label'] = v
                        continue
                    elif k == 'type':
                        graph.nodes[node_id]['label'] = v
                        continue

                    sub_id = json_node(graph, v, prefix + ['%s' % k])
                    if sub_id:
                        graph.add_edge(node_id, sub_id, label=('%s' % k))
            else:
                if isinstance(element, basestring) and \
                   element.startswith('http'):
                    element = element.split('/')[-1]
                graph.add_node(node_id, label=('%s' % str(element)))
            return node_id
        jd = self.to_json()
        graph = networkx.DiGraph()
        json_node(graph, jd, ['%s' % self.uuid])
        return graph

    def make_generic_copy(self, deeply=False):
        """Make a new matching Statement with no provenance.

        All agents and other attributes besides evidence, uuid, supports, and
        supported_by will be copied over, and a new uuid will be assigned.
        Thus, the new Statement will satisfy `new_stmt.matches(old_stmt)`.

        If `deeply` is set to True, all the attributes will be deep-copied,
        which is comparatively slow. Otherwise, attributes of this statement
        may be altered by changes to the new matching statement.
        """
        if deeply:
            kwargs = deepcopy(self.__dict__)
        else:
            kwargs = self.__dict__.copy()
        for attr in ['evidence', 'uuid', 'supports', 'supported_by',
                     'is_activation']:
            kwargs.pop(attr, None)
        my_belief = kwargs.pop('belief', 1)
        my_hash = kwargs.pop('_full_hash', None)
        my_shallow_hash = kwargs.pop('_shallow_hash', None)
        for attr in self._agent_order:
            attr_value = kwargs.get(attr)
            if isinstance(attr_value, list):
                kwargs[attr] = sorted_agents(attr_value)
        new_instance = self.__class__(**kwargs)
        new_instance._full_hash = my_hash
        new_instance._shallow_hash = my_shallow_hash
        new_instance.belief = my_belief
        return new_instance

    def flip_polarity(self, agent_idx=None):
        """If applicable, flip the polarity of the statement"""
        pass


@python_2_unicode_compatible
class Modification(Statement):
    """Generic statement representing the modification of a protein.

    Parameters
    ----------
    enz : :py:class:`indra.statement.Agent`
        The enzyme involved in the modification.
    sub : :py:class:`indra.statement.Agent`
        The substrate of the modification.
    residue : str or None
        The amino acid residue being modified, or None if it is unknown or
        unspecified.
    position : str or None
        The position of the modified amino acid, or None if it is unknown or
        unspecified.
    evidence : None or :py:class:`Evidence` or list of :py:class:`Evidence`
        Evidence objects in support of the modification.
    """

    _agent_order = ['enz', 'sub']

    def __init__(self, enz, sub, residue=None, position=None, evidence=None):
        super(Modification, self).__init__(evidence)
        self.enz = enz
        self.sub = sub
        self.residue = get_valid_residue(residue)
        if isinstance(position, int):
            self.position = str(position)
        else:
            self.position = position

    def matches_key(self):
        if self.enz is None:
            enz_key = None
        else:
            enz_key = self.enz.matches_key()
        key = (stmt_type(self, True), enz_key, self.sub.matches_key(),
               str(self.residue), str(self.position))
        return mk_str(key)

    def set_agent_list(self, agent_list):
        if len(agent_list) != 2:
            raise ValueError("Modification has two agents in agent_list.")
        self.enz = agent_list[0]
        self.sub = agent_list[1]

    def refinement_of(self, other, ontology, entities_refined=False):
        # Make sure the statement types match
        if stmt_type(self) != stmt_type(other):
            return False

        # Check agent arguments
        if self.enz is None and other.enz is None:
            enz_refinement = True
        elif self.enz is None and other.enz is not None:
            enz_refinement = False
        elif self.enz is not None and other.enz is None:
            enz_refinement = True
        else:
            enz_refinement = self.enz.refinement_of(other.enz, ontology,
                entities_refined=entities_refined)
        sub_refinement = self.sub.refinement_of(other.sub, ontology,
            entities_refined=entities_refined)
        if not (enz_refinement and sub_refinement):
            return False
        # For this to be a refinement of the other, the modifications either
        # have to match or have this one be a subtype of the other; in
        # addition, the sites have to match, or this one has to have site
        # information and the other one not.
        residue_matches = (other.residue is None or
                           (self.residue == other.residue))
        position_matches = (other.position is None or
                            (self.position == other.position))
        return residue_matches and position_matches

    def equals(self, other):
        matches = super(Modification, self).equals(other)
        matches = (matches and (self.residue == other.residue)
                   and (self.position == other.position))
        return matches

    def contradicts(self, other, ontology):
        # If the modifications are not the opposite polarity of the
        # same subtype
        if not modclass_to_inverse[self.__class__] == other.__class__:
            return False
        # Skip all instances of not fully specified modifications
        agents = (self.enz, self.sub, other.enz, other.sub)
        if not all(a is not None for a in agents):
            return False
        # If the entities don't match, they can't be contradicting
        # Here we check pairs of agents at each "position" and
        # make sure they are the same or they are refinements of each other
        for self_agent, other_agent in zip(self.agent_list(),
                                           other.agent_list()):
            if not (self_agent.entity_matches(other_agent) or \
                    self_agent.refinement_of(other_agent, ontology) or \
                    other_agent.refinement_of(self_agent, ontology)):
                return False
        # At this point the entities definitely match so we need to
        # check the specific site that is being modified
        if self.residue == other.residue and self.position == other.position:
            return True
        else:
            return False

    def _get_mod_condition(self):
        """Return a ModCondition corresponding to this Modification."""
        mod_type = modclass_to_modtype[self.__class__]
        if isinstance(self, RemoveModification):
            mod_type = modtype_to_inverse[mod_type]
        mc = ModCondition(mod_type, self.residue, self.position, True)
        return mc

    def to_json(self, use_sbo=False, matches_fun=None):
        generic = super(Modification, self).to_json(use_sbo, matches_fun)
        json_dict = _o(type=generic['type'])
        if self.enz is not None:
            json_dict['enz'] = self.enz.to_json()
            if use_sbo:
                # enzymatic catalyst
                json_dict['enz']['sbo'] = \
                    'https://identifiers.org/SBO:0000460'
        if self.sub is not None:
            json_dict['sub'] = self.sub.to_json()
            if use_sbo:
                # substrate
                json_dict['sub']['sbo'] = \
                    'https://identifiers.org/SBO:0000015'
        if self.residue is not None:
            json_dict['residue'] = self.residue
        if self.position is not None:
            json_dict['position'] = self.position
        json_dict.update(generic)
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        enz = json_dict.get('enz')
        sub = json_dict.get('sub')
        residue = json_dict.get('residue')
        position = json_dict.get('position')
        evidence = json_dict.get('evidence', [])
        if enz:
            enz = Agent._from_json(enz)
        if sub:
            sub = Agent._from_json(sub)
        stmt = cls(enz, sub, residue, position)
        return stmt

    def __str__(self):
        res_str = (', %s' % self.residue) if self.residue is not None else ''
        pos_str = (', %s' % self.position) if self.position is not None else ''
        s = ("%s(%s, %s%s%s)" %
                 (type(self).__name__, self.enz, self.sub, res_str, pos_str))
        return s


class AddModification(Modification):
    pass


class RemoveModification(Modification):
    pass


@python_2_unicode_compatible
class SelfModification(Statement):
    """Generic statement representing the self-modification of a protein.

    Parameters
    ----------
    enz : :py:class:`indra.statement.Agent`
        The enzyme involved in the modification, which is also the substrate.
    residue : str or None
        The amino acid residue being modified, or None if it is unknown or
        unspecified.
    position : str or None
        The position of the modified amino acid, or None if it is unknown or
        unspecified.
    evidence : None or :py:class:`Evidence` or list of :py:class:`Evidence`
        Evidence objects in support of the modification.
    """

    _agent_order = ['enz']

    def __init__(self, enz, residue=None, position=None, evidence=None):
        super(SelfModification, self).__init__(evidence)
        self.enz = enz
        self.residue = get_valid_residue(residue)
        if isinstance(position, int):
            self.position = str(position)
        else:
            self.position = position

    def __str__(self):
        res_str = (', %s' % self.residue) if self.residue is not None else ''
        pos_str = (', %s' % self.position) if self.position is not None else ''
        s = ("%s(%s%s%s)" %
             (type(self).__name__, self.enz, res_str, pos_str))
        return s

    def matches_key(self):
        key = (stmt_type(self, True), self.enz.matches_key(),
               str(self.residue), str(self.position))
        return mk_str(key)

    def set_agent_list(self, agent_list):
        if len(agent_list) != 1:
            raise ValueError("SelfModification has one agent.")
        self.enz = agent_list[0]

    def refinement_of(self, other, ontology, entities_refined=False):
        # Make sure the statement types match
        if stmt_type(self) != stmt_type(other):
            return False

        # Check agent arguments
        if not self.enz.refinement_of(other.enz, ontology,
                                      entities_refined=entities_refined):
            return False
        # For this to be a refinement of the other, the modifications either
        # have to match or have this one be a subtype of the other; in
        # addition, the sites have to match, or this one has to have site
        # information and the other one not.
        residue_matches = (other.residue is None or
                           (self.residue == other.residue))
        position_matches = (other.position is None or
                            (self.position == other.position))
        return residue_matches and position_matches

    def equals(self, other):
        matches = super(SelfModification, self).equals(other)
        matches = (matches and self.residue == other.residue
                   and self.position == other.position)
        return matches

    def _get_mod_condition(self):
        """Return a ModCondition corresponding to this Modification."""
        mod_type = modclass_to_modtype[self.__class__]
        mc = ModCondition(mod_type, self.residue, self.position, True)
        return mc

    def to_json(self, use_sbo=False, matches_fun=None):
        generic = super(SelfModification, self).to_json(use_sbo, matches_fun)
        json_dict = _o(type=generic['type'])
        if self.enz is not None:
            json_dict['enz'] = self.enz.to_json()
            if use_sbo:
                # enzymatic catalyst
                json_dict['enz']['sbo'] = \
                    'https://identifiers.org/SBO:0000460'
        if self.residue is not None:
            json_dict['residue'] = self.residue
        if self.position is not None:
            json_dict['position'] = self.position
        json_dict.update(generic)
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        enz = json_dict.get('enz')
        residue = json_dict.get('residue')
        position = json_dict.get('position')
        if enz:
            enz = Agent._from_json(enz)
        stmt = cls(enz, residue, position)
        return stmt


class Phosphorylation(AddModification):
    """Phosphorylation modification.

    Examples
    --------
    MEK (MAP2K1) phosphorylates ERK (MAPK1) at threonine 185:

    >>> mek = Agent('MAP2K1')
    >>> erk = Agent('MAPK1')
    >>> phos = Phosphorylation(mek, erk, 'T', '185')
    """


class Autophosphorylation(SelfModification):
    """Intramolecular autophosphorylation, i.e., in *cis*.

    Examples
    --------
    p38 bound to TAB1 cis-autophosphorylates itself (see :pmid:`19155529`).

    >>> tab1 = Agent('TAB1')
    >>> p38_tab1 = Agent('P38', bound_conditions=[BoundCondition(tab1)])
    >>> autophos = Autophosphorylation(p38_tab1)
    """


class Transphosphorylation(SelfModification):
    """Autophosphorylation in *trans.*

    Transphosphorylation assumes that a kinase is already bound to a substrate
    (usually of the same molecular species), and phosphorylates it in an
    intra-molecular fashion. The enz property of the statement must have
    exactly one bound_conditions entry, and we assume that enz phosphorylates
    this molecule. The bound_neg property is ignored here.
    """


class Dephosphorylation(RemoveModification):
    """Dephosphorylation modification.

    Examples
    --------
    DUSP6 dephosphorylates ERK (MAPK1) at T185:

    >>> dusp6 = Agent('DUSP6')
    >>> erk = Agent('MAPK1')
    >>> dephos = Dephosphorylation(dusp6, erk, 'T', '185')
    """


class Hydroxylation(AddModification):
    """Hydroxylation modification."""


class Dehydroxylation(RemoveModification):
    """Dehydroxylation modification."""


class Sumoylation(AddModification):
    """Sumoylation modification."""


class Desumoylation(RemoveModification):
    """Desumoylation modification."""


class Acetylation(AddModification):
    """Acetylation modification."""


class Deacetylation(RemoveModification):
    """Deacetylation modification."""


class Glycosylation(AddModification):
    """Glycosylation modification."""


class Deglycosylation(RemoveModification):
    """Deglycosylation modification."""


class Ribosylation(AddModification):
    """Ribosylation modification."""


class Deribosylation(RemoveModification):
    """Deribosylation modification."""


class Ubiquitination(AddModification):
    """Ubiquitination modification."""


class Deubiquitination(RemoveModification):
    """Deubiquitination modification."""


class Farnesylation(AddModification):
    """Farnesylation modification."""


class Defarnesylation(RemoveModification):
    """Defarnesylation modification."""


class Geranylgeranylation(AddModification):
    """Geranylgeranylation modification."""


class Degeranylgeranylation(RemoveModification):
    """Degeranylgeranylation modification."""


class Palmitoylation(AddModification):
    """Palmitoylation modification."""


class Depalmitoylation(RemoveModification):
    """Depalmitoylation modification."""


class Myristoylation(AddModification):
    """Myristoylation modification."""


class Demyristoylation(RemoveModification):
    """Demyristoylation modification."""


class Methylation(AddModification):
    """Methylation modification."""


class Demethylation(RemoveModification):
    """Demethylation modification."""


class Carboxylation(AddModification):
    """Carboxylation modification."""


class Decarboxylation(RemoveModification):
    """Decarboxylation modification."""


class Amidation(AddModification):
    """Amidation modification."""


class Deamidation(RemoveModification):
    """Deamidation modification."""


class Sulfation(AddModification):
    """Sulfation modification."""


class Desulfation(RemoveModification):
    """Desulfation modification."""


class Neddylation(AddModification):
    """Neddylation modification."""


class Deneddylation(RemoveModification):
    """Deneddylation modification."""


class Alkylation(AddModification):
    """Alkylation modification."""


class Dealkylation(RemoveModification):
    """Dealkylation modification."""


class Prenylation(AddModification):
    """Prenylation modification."""


class Deprenylation(RemoveModification):
    """Deprenylation modification."""


@python_2_unicode_compatible
class RegulateActivity(Statement):
    """Regulation of activity.

    This class implements shared functionality of Activation and Inhibition
    statements and it should not be instantiated directly.
    """

    # The constructor here is an abstractmethod so that this class cannot
    # be directly instantiated.
    __metaclass__ = abc.ABCMeta

    _agent_order = ['subj', 'obj']

    @abc.abstractmethod
    def __init__(self):
        pass

    def __setstate__(self, state):
        if 'subj_activity' in state:
            logger.warning('Pickle file is out of date!')
        state.pop('subj_activity', None)
        self.__dict__.update(state)

    def matches_key(self):
        key = (stmt_type(self, True), self.subj.matches_key(),
               self.obj.matches_key(), str(self.obj_activity),
               str(self.is_activation))
        return mk_str(key)

    def set_agent_list(self, agent_list):
        if len(agent_list) != 2:
            raise ValueError("%s has two agents." % self.__class__.__name__)
        self.subj = agent_list[0]
        self.obj = agent_list[1]

    def refinement_of(self, other, ontology, entities_refined=False):
        # Make sure the statement types match
        if stmt_type(self) != stmt_type(other):
            return False
        if self.is_activation != other.is_activation:
            return False
        if not (self.subj.refinement_of(other.subj, ontology,
                                        entities_refined=entities_refined) and
                self.obj.refinement_of(other.obj, ontology,
                                       entities_refined=entities_refined)):
            return False

        return ((self.obj_activity == other.obj_activity) or
                 ontology.isa('INDRA_ACTIVITIES', self.obj_activity,
                              'INDRA_ACTIVITIES', other.obj_activity))

    def contradicts(self, other, ontology):
        # If they aren't opposite classes, it's not a contradiction
        if {self.__class__, other.__class__} != {Activation, Inhibition}:
            return False

        # If they aren't opposite classes, it's not a contradiction
        if self.is_activation == other.is_activation:
            return False
        # Skip all instances of not fully specified statements
        agents = (self.subj, self.obj, other.subj, other.obj)
        if not all(a is not None for a in agents):
            return False
        # If the entities don't match, they can't be contradicting
        # Here we check pairs of agents at each "position" and
        # make sure they are the same or they are refinements of each other
        for self_agent, other_agent in zip(self.agent_list(),
                                           other.agent_list()):
            if not (self_agent.entity_matches(other_agent) or \
                    self_agent.refinement_of(other_agent, ontology) or \
                    other_agent.refinement_of(self_agent, ontology)):
                return False
        # Otherwise they are contradicting
        return True

    def to_json(self, use_sbo=False, matches_fun=None):
        generic = super(RegulateActivity, self).to_json(use_sbo, matches_fun)
        json_dict = _o(type=generic['type'])
        if self.subj is not None:
            json_dict['subj'] = self.subj.to_json()
            if use_sbo:
                if self.is_activation:
                    json_dict['subj']['sbo'] = \
                        'https://identifiers.org/SBO:0000459'  # stimulator
                else:
                    json_dict['subj']['sbo'] = \
                        'https://identifiers.org/SBO:0000020'  # inhibitor
        if self.obj is not None:
            json_dict['obj'] = self.obj.to_json()
            if use_sbo:
                if self.is_activation:
                    json_dict['obj']['sbo'] = \
                        'https://identifiers.org/SBO:0000643'  # stimulated
                else:
                    json_dict['obj']['sbo'] = \
                        'https://identifiers.org/SBO:0000642'  # inhibited
        if self.obj_activity is not None:
            json_dict['obj_activity'] = self.obj_activity
        json_dict.update(generic)
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        subj = json_dict.get('subj')
        obj = json_dict.get('obj')
        obj_activity = json_dict.get('obj_activity')
        if subj:
            subj = Agent._from_json(subj)
        if obj:
            obj = Agent._from_json(obj)
        stmt = cls(subj, obj, obj_activity)
        return stmt

    def __str__(self):
        obj_act_str = ', %s' % self.obj_activity if \
            self.obj_activity != 'activity' else ''
        s = ("%s(%s, %s%s)" %
             (type(self).__name__, self.subj,
              self.obj, obj_act_str))
        return s

    def __repr__(self):
        return self.__str__()

    def equals(self, other):
        matches = super(RegulateActivity, self).equals(other)
        matches = (matches and self.obj_activity == other.obj_activity
                   and self.is_activation == other.is_activation)
        return matches

    def _get_activity_condition(self):
        """Return ActivityCondition corresponding to this RegulateActivity."""
        return ActivityCondition(self.obj_activity, self.is_activation)


class Inhibition(RegulateActivity):
    """Indicates that a protein inhibits or deactivates another protein.

    This statement is intended to be used for physical interactions where the
    mechanism of inhibition is not explicitly specified, which is often the
    case for descriptions of mechanisms extracted from the literature.

    Parameters
    ----------
    subj : :py:class:`Agent`
        The agent responsible for the change in activity, i.e., the "upstream"
        node.
    obj : :py:class:`Agent`
        The agent whose activity is influenced by the subject, i.e., the
        "downstream" node.
    obj_activity : Optional[str]
        The activity of the obj Agent that is affected, e.g., its "kinase"
        activity.
    evidence : None or :py:class:`Evidence` or list of :py:class:`Evidence`
        Evidence objects in support of the modification.
    """
    def __init__(self, subj, obj, obj_activity='activity', evidence=None):
        super(RegulateActivity, self).__init__(evidence)
        self.subj = subj
        self.obj = obj
        if obj_activity not in activity_types:
            logger.warning('Invalid activity type: %s' % obj_activity)
        self.obj_activity = obj_activity
        self.is_activation = False


class Activation(RegulateActivity):
    """Indicates that a protein activates another protein.

    This statement is intended to be used for physical interactions where the
    mechanism of activation is not explicitly specified, which is often the
    case for descriptions of mechanisms extracted from the literature.

    Parameters
    ----------
    subj : :py:class:`Agent`
        The agent responsible for the change in activity, i.e., the "upstream"
        node.
    obj : :py:class:`Agent`
        The agent whose activity is influenced by the subject, i.e., the
        "downstream" node.
    obj_activity : Optional[str]
        The activity of the obj Agent that is affected, e.g., its "kinase"
        activity.
    evidence : None or :py:class:`Evidence` or list of :py:class:`Evidence`
        Evidence objects in support of the modification.

    Examples
    --------

    MEK (MAP2K1) activates the kinase activity of ERK (MAPK1):

    >>> mek = Agent('MAP2K1')
    >>> erk = Agent('MAPK1')
    >>> act = Activation(mek, erk, 'kinase')
    """

    def __init__(self, subj, obj, obj_activity='activity', evidence=None):
        super(RegulateActivity, self).__init__(evidence)
        self.subj = subj
        self.obj = obj
        if obj_activity not in activity_types:
            logger.warning('Invalid activity type: %s' % obj_activity)
        self.obj_activity = obj_activity
        self.is_activation = True


class GtpActivation(Activation):
    pass


@python_2_unicode_compatible
class ActiveForm(Statement):
    """Specifies conditions causing an Agent to be active or inactive.

    Types of conditions influencing a specific type of biochemical activity can
    include modifications, bound Agents, and mutations.

    Parameters
    ----------
    agent : :py:class:`Agent`
        The Agent in a particular active or inactive state. The sets
        of ModConditions, BoundConditions, and MutConditions on the given
        Agent instance indicate the relevant conditions.
    activity : str
        The type of activity influenced by the given set of conditions,
        e.g., "kinase".
    is_active : bool
        Whether the conditions are activating (True) or inactivating (False).
    """

    _agent_order = ['agent']

    def __init__(self, agent, activity, is_active, evidence=None):
        super(ActiveForm, self).__init__(evidence)
        self.agent = agent
        if agent.activity is not None:
            logger.warning('Agent in ActiveForm should not have ' +
                           'ActivityConditions.')
            agent.activity = None
        if activity not in activity_types:
            logger.warning('Invalid activity type: %s' % activity)
        self.activity = activity
        self.is_active = is_active

    def matches_key(self):
        key = (stmt_type(self, True), self.agent.matches_key(),
               str(self.activity), str(self.is_active))
        return mk_str(key)

    def set_agent_list(self, agent_list):
        if len(agent_list) != 1:
            raise ValueError("ActiveForm has one agent.")
        self.agent = agent_list[0]

    def refinement_of(self, other, ontology, entities_refined=False):
        # Make sure the statement types match
        if stmt_type(self) != stmt_type(other):
            return False

        # Check agent arguments
        if not self.agent.refinement_of(other.agent, ontology,
                                        entities_refined=entities_refined):
            return False

        # Make sure that the relationships and activities match
        if (self.is_active == other.is_active) and \
            (self.activity == other.activity
             or ontology.isa('INDRA_ACTIVITIES', self.activity,
                             'INDRA_ACTIVITIES', other.activity)):
                return True
        else:
            return False

    def contradicts(self, other, ontology):
        # Make sure the statement types match
        if stmt_type(self) != stmt_type(other):
            return False
        # Check that the polarity is constradicting up front
        # TODO: we could also check for cases where the polarities are
        # the same but some of the state conditions have opposite
        # polarities, for instance, if the presence/lack of the
        # same modification activates a given agent, that could be
        # considered a contradiction.
        if self.is_active == other.is_active:
            return False
        # Check that the activity types are the same
        # TODO: we could check for refinements here
        if self.activity != other.activity:
            return False
        # If the agents are exactly the same, this is a contradiction
        if self.agent.matches_key() == other.agent.matches_key():
            return True
        # Otherwise, if the two agents are related at the level of entities
        # and their state is exactly the same, then they contradict
        if self.agent.state_matches_key() == other.agent.state_matches_key():
            if self.agent.isa(other.agent, ontology) or \
                other.agent.isa(self.agent, ontology):
                return True
        return False

    def to_json(self, use_sbo=False, matches_fun=None):
        generic = super(ActiveForm, self).to_json(use_sbo, matches_fun)
        json_dict = _o(type=generic['type'])
        json_dict.update({'agent': self.agent.to_json(),
                          'activity': self.activity,
                          'is_active': self.is_active})
        if use_sbo:
            json_dict['agent']['sbo'] = \
                'https://identifiers.org/SBO:0000644'  # modified
        json_dict.update(generic)
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        agent = json_dict.get('agent')
        if agent:
            agent = Agent._from_json(agent)
        else:
            logger.error('ActiveForm statement missing agent')
            return None
        activity = json_dict.get('activity')
        is_active = json_dict.get('is_active')
        if activity is None:
            logger.warning('ActiveForm activity missing, defaulting ' +
                           'to `activity`')
            activity = 'activity'
        if is_active is None:
            logger.warning('ActiveForm is_active missing, defaulting ' +
                           'to True')
            is_active = True
        stmt = cls(agent, activity, is_active)
        return stmt

    def __str__(self):
        s = ("ActiveForm(%s, %s, %s)" %
             (self.agent, self.activity, self.is_active))
        return s

    def equals(self, other):
        matches = super(ActiveForm, self).equals(other)
        matches = (matches and self.activity == other.activity
                   and self.is_active == other.is_active)
        return matches


@python_2_unicode_compatible
class HasActivity(Statement):
    """States that an Agent has or doesn't have a given activity type.

    With this Statement, one cane express that a given protein is a kinase, or,
    for instance, that it is a transcription factor. It is also possible to
    construct negative statements with which one epxresses, for instance,
    that a given protein is not a kinase.

    Parameters
    ----------
    agent : :py:class:`Agent`
        The Agent that that statement is about. Note that the detailed state
        of the Agent is not relevant for this type of statement.
    activity : str
        The type of activity, e.g., "kinase".
    has_activity : bool
        Whether the given Agent has the given activity (True) or
        not (False).
    """

    _agent_order = ['agent']

    def __init__(self, agent, activity, has_activity, evidence=None):
        super(HasActivity, self).__init__(evidence)
        if agent.activity is not None:
            logger.warning('Agent in HasActivity should not have ' +
                           'ActivityConditions.')
            agent.activity = None
        self.agent = agent
        if activity not in activity_types:
            logger.warning('Invalid activity type: %s' % activity)
        self.activity = activity
        self.has_activity = has_activity

    def matches_key(self):
        key = (stmt_type(self, True), self.agent.matches_key(),
               str(self.activity), str(self.has_activity))
        return mk_str(key)

    def set_agent_list(self, agent_list):
        if len(agent_list) != 1:
            raise ValueError("HasActivity has one agent.")
        self.agent = agent_list[0]

    def refinement_of(self, other, ontology, entities_refined=False):
        # Make sure the statement types match
        if stmt_type(self) != stmt_type(other):
            return False

        # Check agent arguments
        if not self.agent.refinement_of(other.agent, ontology,
                                        entities_refined=entities_refined):
            return False

        # Make sure that the relationships and activities match
        if (self.has_activity == other.has_activity) and \
            (self.activity == other.activity
             or ontology.isa('INDRA_ACTIVITIES', self.activity,
                             'INDRA_ACTIVITIES', other.activity)):
                return True
        else:
            return False

    def __str__(self):
        s = ("HasActivity(%s, %s, %s)" %
             (self.agent, self.activity, self.has_activity))
        return s

    def equals(self, other):
        matches = super(HasActivity, self).equals(other)
        matches = (matches and self.activity == other.activity
                   and self.has_activity == other.has_activity)
        return matches


@python_2_unicode_compatible
class Gef(Statement):
    """Exchange of GTP for GDP on a small GTPase protein mediated by a GEF.

    Represents the generic process by which a guanosine exchange factor (GEF)
    catalyzes nucleotide exchange on a GTPase protein.

    Parameters
    ----------
    gef : :py:class:`Agent`
        The guanosine exchange factor.
    ras : :py:class:`Agent`
        The GTPase protein.

    Examples
    --------
    SOS1 catalyzes nucleotide exchange on KRAS:

    >>> sos = Agent('SOS1')
    >>> kras = Agent('KRAS')
    >>> gef = Gef(sos, kras)
    """

    _agent_order = ['gef', 'ras']

    def __init__(self, gef, ras, evidence=None):
        super(Gef, self).__init__(evidence)
        self.gef = gef
        self.ras = ras

    def matches_key(self):
        key = (stmt_type(self, True), self.gef.matches_key(),
               self.ras.matches_key())
        return mk_str(key)

    def set_agent_list(self, agent_list):
        if len(agent_list) != 2:
            raise ValueError("Gef has two agents.")
        self.gef = agent_list[0]
        self.ras = agent_list[1]

    def __str__(self):
        s = "Gef(%s, %s)" % (self.gef.name, self.ras.name)
        return s

    def refinement_of(self, other, ontology, entities_refined=False):
        # Make sure the statement types match
        if stmt_type(self) != stmt_type(other):
            return False
        # Check the GEF
        if self.gef.refinement_of(other.gef, ontology,
                                  entities_refined=entities_refined) and \
           self.ras.refinement_of(other.ras, ontology,
                                  entities_refined=entities_refined):
            return True
        else:
            return False

    def equals(self, other):
        matches = super(Gef, self).equals(other)
        return matches

    def to_json(self, use_sbo=False, matches_fun=None):
        generic = super(Gef, self).to_json(use_sbo, matches_fun)
        json_dict = _o(type=generic['type'])
        if self.gef is not None:
            json_dict['gef'] = self.gef.to_json()
            if use_sbo:
                json_dict['gef']['sbo'] = \
                    'https://identifiers.org/SBO:0000013'  # catalyst
        if self.ras is not None:
            json_dict['ras'] = self.ras.to_json()
            if use_sbo:
                json_dict['ras']['sbo'] = \
                    'https://identifiers.org/SBO:0000015'  # substrate
        json_dict.update(generic)
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        gef = json_dict.get('gef')
        ras = json_dict.get('ras')
        if gef:
            gef = Agent._from_json(gef)
        if ras:
            ras = Agent._from_json(ras)
        stmt = cls(gef, ras)
        return stmt


@python_2_unicode_compatible
class Gap(Statement):
    """Acceleration of a GTPase protein's GTP hydrolysis rate by a GAP.

    Represents the generic process by which a GTPase activating protein (GAP)
    catalyzes GTP hydrolysis by a particular small GTPase protein.

    Parameters
    ----------
    gap : :py:class:`Agent`
        The GTPase activating protein.
    ras : :py:class:`Agent`
        The GTPase protein.

    Examples
    --------
    RASA1 catalyzes GTP hydrolysis on KRAS:

    >>> rasa1 = Agent('RASA1')
    >>> kras = Agent('KRAS')
    >>> gap = Gap(rasa1, kras)
    """

    _agent_order = ['gap', 'ras']

    def __init__(self, gap, ras, evidence=None):
        super(Gap, self).__init__(evidence)
        self.gap = gap
        self.ras = ras

    def matches_key(self):
        key = (stmt_type(self, True), self.gap.matches_key(),
               self.ras.matches_key())
        return mk_str(key)

    def set_agent_list(self, agent_list):
        if len(agent_list) != 2:
            raise ValueError("Gap has two agents.")
        self.gap = agent_list[0]
        self.ras = agent_list[1]

    def refinement_of(self, other, ontology, entities_refined=False):
        # Make sure the statement types match
        if stmt_type(self) != stmt_type(other):
            return False

        # Check the GAP
        if self.gap.refinement_of(other.gap, ontology,
                                  entities_refined=entities_refined) and \
           self.ras.refinement_of(other.ras, ontology,
                                  entities_refined=entities_refined):
            return True
        else:
            return False

    def __str__(self):
        s = "Gap(%s, %s)" % (self.gap.name, self.ras.name)
        return s

    def equals(self, other):
        matches = super(Gap, self).equals(other)
        return matches

    def to_json(self, use_sbo=False, matches_fun=None):
        generic = super(Gap, self).to_json(use_sbo, matches_fun)
        json_dict = _o(type=generic['type'])
        if self.gap is not None:
            json_dict['gap'] = self.gap.to_json()
            if use_sbo:
                json_dict['gap']['sbo'] = \
                    'https://identifiers.org/SBO:0000013'  # catalyst
        if self.ras is not None:
            json_dict['ras'] = self.ras.to_json()
            if use_sbo:
                json_dict['ras']['sbo'] = \
                    'https://identifiers.org/SBO:0000015'  # substrate
        json_dict.update(generic)
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        gap = json_dict.get('gap')
        ras = json_dict.get('ras')
        if gap:
            gap = Agent._from_json(gap)
        if ras:
            ras = Agent._from_json(ras)
        stmt = cls(gap, ras)
        return stmt


@python_2_unicode_compatible
class Complex(Statement):
    """A set of proteins observed to be in a complex.

    Parameters
    ----------
    members : list of :py:class:`Agent`
        The set of proteins in the complex.

    Examples
    --------
    BRAF is observed to be in a complex with RAF1:

    >>> braf = Agent('BRAF')
    >>> raf1 = Agent('RAF1')
    >>> cplx = Complex([braf, raf1])
    """

    _agent_order = ['members']

    def __init__(self, members, evidence=None):
        super(Complex, self).__init__(evidence)
        self.members = members

    def matches_key(self):
        key = (stmt_type(self, True), tuple(m.matches_key()
                                            for m in self.sorted_members()))
        return mk_str(key)

    def sorted_members(self):
        return sorted(self.members, key=lambda x: x.matches_key())

    def entities_match_key(self):
        key = tuple(a.entity_matches_key() if a is not None
                    else None for a in sorted(self.members,
                                              key=lambda x: x.matches_key()))
        return key

    def set_agent_list(self, agent_list):
        self.members = agent_list

    def __str__(self):
        s = '%s(%s)' % (type(self).__name__,
                        (', '.join([('%s' % m) for m in self.members])))
        return s

    def refinement_of(self, other, ontology, entities_refined=False):
        # Make sure the statement types match
        if stmt_type(self) != stmt_type(other):
            return False
        # Make sure the length of the members list is the same. Note that this
        # treats Complex([A, B, C]) as distinct from Complex([A, B]), rather
        # than as a refinement.
        if len(self.members) != len(other.members):
            return False

        def match_members(self_members, other_members):
            # First build a bipartite graph of refinement links
            G = networkx.Graph()
            for (self_idx, self_member), (other_idx, other_member) in \
                itertools.product(enumerate(self_members),
                                  enumerate(other_members)):
                if self_member.refinement_of(other_member, ontology):
                    G.add_edge('S%d' % self_idx, 'O%d' % other_idx)
            # Then find a maximal matching in the bipartite graph
            match = networkx.algorithms.max_weight_matching(G)
            # If every member has a pair, it is a valid refinement
            return len(match) == len(self_members)

        return match_members(self.members, other.members)

    def to_json(self, use_sbo=False, matches_fun=None):
        generic = super(Complex, self).to_json(use_sbo, matches_fun)
        json_dict = _o(type=generic['type'])
        members = [m.to_json() for m in self.members]
        json_dict['members'] = members
        json_dict.update(generic)
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        members = json_dict.get('members')
        members = [Agent._from_json(m) for m in members]
        stmt = cls(members)
        return stmt


@python_2_unicode_compatible
class Translocation(Statement):
    """The translocation of a molecular agent from one location to another.

    Parameters
    ----------
    agent : :py:class:`Agent`
        The agent which translocates.
    from_location : Optional[str]
        The location from which the agent translocates. This must
        be a valid GO cellular component name (e.g. "cytoplasm")
        or ID (e.g. "GO:0005737").
    to_location : Optional[str]
        The location to which the agent translocates. This must
        be a valid GO cellular component name or ID.
    """
    _agent_order = ['agent']

    def __init__(self, agent, from_location=None, to_location=None,
                 evidence=None):
        super(Translocation, self).__init__(evidence)
        self.agent = agent
        self.from_location = from_location
        self.to_location = to_location

    def set_agent_list(self, agent_list):
        if len(agent_list) != 1:
            raise ValueError("Translocation has 1 agent")
        self.agent = agent_list[0]

    def __str__(self):
        s = ("Translocation(%s, %s, %s)" %
             (self.agent, self.from_location, self.to_location))
        return s

    def refinement_of(self, other, ontology, entities_refined=False):
        from indra.databases import go_client
        # Make sure the statement types match
        if stmt_type(self) != stmt_type(other):
            return False
        # Check several conditions for refinement
        ref1 = self.agent.refinement_of(other.agent, ontology,
                                        entities_refined=entities_refined)
        ofl = go_client.get_go_id_from_label(other.from_location)
        sfl = go_client.get_go_id_from_label(self.from_location)
        otl = go_client.get_go_id_from_label(other.to_location)
        stl = go_client.get_go_id_from_label(self.to_location)
        ref2 = (other.from_location is None or
                self.from_location == other.from_location or
                ontology.isa_or_partof('GO', sfl, 'GO', ofl))
        ref3 = (other.to_location is None or
                self.to_location == other.to_location or
                ontology.isa_or_partof('GO', stl, 'GO', otl))
        return ref1 and ref2 and ref3

    def equals(self, other):
        matches = super(Translocation, self).equals(other)
        matches = matches and (self.from_location == other.from_location)
        matches = matches and (self.to_location == other.to_location)
        return matches

    def matches_key(self):
        key = (stmt_type(self, True), self.agent.matches_key(),
               str(self.from_location), str(self.to_location))
        return mk_str(key)

    def to_json(self, use_sbo=False, matches_fun=None):
        generic = super(Translocation, self).to_json(use_sbo, matches_fun)
        json_dict = _o(type=generic['type'])
        json_dict['agent'] = self.agent.to_json()
        if self.from_location is not None:
            json_dict['from_location'] = self.from_location
        if self.to_location is not None:
            json_dict['to_location'] = self.to_location
        json_dict.update(generic)
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        agent = json_dict.get('agent')
        if agent:
            agent = Agent._from_json(agent)
        else:
            logger.error('Translocation statement missing agent')
            return None
        from_location = json_dict.get('from_location')
        to_location = json_dict.get('to_location')
        stmt = cls(agent, from_location, to_location)
        return stmt


@python_2_unicode_compatible
class RegulateAmount(Statement):
    """Superclass handling operations on directed, two-element interactions."""
    _agent_order = ['subj', 'obj']

    def __init__(self, subj, obj, evidence=None):
        super(RegulateAmount, self).__init__(evidence)
        self.subj = subj
        if obj is None:
            raise ValueError('Object of %s cannot be None.' %
                             type(self).__name__)
        self.obj = obj

    def matches_key(self):
        if self.subj is None:
            subj_key = None
        else:
            subj_key = self.subj.matches_key()
        key = (stmt_type(self, True), subj_key, self.obj.matches_key())
        return mk_str(key)

    def set_agent_list(self, agent_list):
        if len(agent_list) != 2:
            raise ValueError("%s has two agents in agent_list." %
                             type(self).__name__)
        self.subj = agent_list[0]
        self.obj = agent_list[1]

    def to_json(self, use_sbo=False, matches_fun=None):
        generic = super(RegulateAmount, self).to_json(use_sbo, matches_fun)
        json_dict = _o(type=generic['type'])
        if self.subj is not None:
            json_dict['subj'] = self.subj.to_json()
            if use_sbo:
                if isinstance(self, IncreaseAmount):
                    json_dict['subj']['sbo'] = \
                        'https://identifiers.org/SBO:0000459'  # stimulator
                else:
                    json_dict['subj']['sbo'] = \
                        'https://identifiers.org/SBO:0000020'  # inhibitor
        if self.obj is not None:
            json_dict['obj'] = self.obj.to_json()
            if use_sbo:
                if isinstance(self, IncreaseAmount):
                    json_dict['obj']['sbo'] = \
                        'https://identifiers.org/SBO:0000011'  # product
                else:
                    json_dict['obj']['sbo'] = \
                        'https://identifiers.org/SBO:0000010'  # reactant
        json_dict.update(generic)
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        subj = json_dict.get('subj')
        obj = json_dict.get('obj')
        if subj:
            subj = Agent._from_json(subj)
        if obj:
            obj = Agent._from_json(obj)
        stmt = cls(subj, obj)
        return stmt

    def refinement_of(self, other, ontology, entities_refined=False):
        # Make sure the statement types match
        if stmt_type(self) != stmt_type(other):
            return False

        # Check agent arguments
        if self.subj is None and other.subj is None:
            subj_refinement = True
        elif self.subj is None and other.subj is not None:
            subj_refinement = False
        elif self.subj is not None and other.subj is None:
            subj_refinement = True
        else:
            subj_refinement = \
                self.subj.refinement_of(other.subj, ontology,
                                        entities_refined=entities_refined)
        obj_refinement = \
            self.obj.refinement_of(other.obj, ontology,
                                   entities_refined=entities_refined)
        return subj_refinement and obj_refinement

    def equals(self, other):
        matches = super(RegulateAmount, self).equals(other)
        return matches

    def contradicts(self, other, ontology):
        # If they aren't opposite classes, it's not a contradiction
        if {self.__class__, other.__class__} != \
                {IncreaseAmount, DecreaseAmount}:
            return False
        # Skip all instances of not fully specified statements
        agents = (self.subj, self.obj, other.subj, other.obj)
        if not all(a is not None for a in agents):
            return False
        # If the entities don't match, they can't be contradicting
        # Here we check pairs of agents at each "position" and
        # make sure they are the same or they are refinements of each other
        for self_agent, other_agent in zip(self.agent_list(),
                                           other.agent_list()):
            if not (self_agent.entity_matches(other_agent) or \
                    self_agent.refinement_of(other_agent, ontology) or \
                    other_agent.refinement_of(self_agent, ontology)):
                return False
        # Otherwise they are contradicting
        return True

    def __str__(self):
        s = ("%s(%s, %s)" % (type(self).__name__, self.subj, self.obj))
        return s


class DecreaseAmount(RegulateAmount):
    """Degradation of a protein, possibly mediated by another protein.

    Note that this statement can also be used to represent inhibitors of
    synthesis (e.g., cycloheximide).

    Parameters
    ----------
    subj : :py:class:`indra.statement.Agent`
        The protein mediating the degradation.
    obj : :py:class:`indra.statement.Agent`
        The protein that is degraded.
    evidence : None or :py:class:`Evidence` or list of :py:class:`Evidence`
        Evidence objects in support of the degradation statement.
    """


class IncreaseAmount(RegulateAmount):
    """Synthesis of a protein, possibly mediated by another protein.

    Parameters
    ----------
    subj : :py:class:`indra.statement.Agent`
        The protein mediating the synthesis.
    obj : :py:class:`indra.statement.Agent`
        The protein that is synthesized.
    evidence : None or :py:class:`Evidence` or list of :py:class:`Evidence`
        Evidence objects in support of the synthesis statement.
    """


class Influence(Statement):
    """An influence on the quantity of a concept of interest.

    Parameters
    ----------
    subj : :py:class:`indra.statement.Event`
        The event which acts as the influencer.
    obj : :py:class:`indra.statement.Event`
        The event which acts as the influencee.
    evidence : None or :py:class:`Evidence` or list of :py:class:`Evidence`
        Evidence objects in support of the statement.
    """
    _agent_order = ['subj', 'obj']

    def agent_list(self, deep_sorted=False):
        return [self.subj.concept, self.obj.concept]

    def __init__(self, subj, obj, evidence=None):
        super(Influence, self).__init__(evidence)
        self.subj = subj
        self.obj = obj

    def refinement_of(self, other, ontology, entities_refined=False):
        # Make sure the statement types match
        if stmt_type(self) != stmt_type(other):
            return False
        delta_refinement = self.delta_refinement_of(other)
        if not delta_refinement:
            return False
        # Check agent arguments
        subj_refinement = self.subj.concept.refinement_of(
            other.subj.concept, ontology, entities_refined=entities_refined)
        obj_refinement = self.obj.concept.refinement_of(
            other.obj.concept, ontology, entities_refined=entities_refined)
        return (subj_refinement and obj_refinement)

    def delta_refinement_of(self, other):
        op = other.overall_polarity()
        sp = self.overall_polarity()
        # If we have "less" polarity here than in other then it's
        # not a refinement.
        if self.polarity_count() < other.polarity_count():
            delta_refinement = False
        # If we have some polarity and the other doesn't then it's
        # always a refinement
        elif sp is not None and op is None:
            delta_refinement = True
        # Otherwise we need to check if the overall polarity matches, if it
        # does then this is a refinement. Otherwise it isn't.
        else:
            delta_refinement = (op == sp)
        return delta_refinement

    def equals(self, other):
        equals = super(Influence, self).equals(other) and \
            self.subj.equals(other.subj) and self.obj.equals(other.obj)
        return equals

    def matches_key(self):
        # With polarities, here, the goal is to match overall polarity
        # if both polarities are given, i.e. +/+ matches -/-. Also, if only
        # one polarity is given, we match the overall polarity e.g.
        # None/+, +/None will match.
        key = (stmt_type(self, True),
               self.subj.matches_key(),
               self.obj.matches_key(),
               self.polarity_count(),
               self.overall_polarity()
               )
        return mk_str(key)

    def contradicts(self, other, ontology):
        # Make sure the statement types match
        if stmt_type(self) != stmt_type(other):
            return False

        # Determine some refinements and opposites up front
        subj_ref = self.subj.concept.refinement_of(other.subj.concept,
                                                   ontology) or \
            other.subj.concept.refinement_of(self.subj.concept, ontology)
        obj_ref = self.obj.concept.refinement_of(other.obj.concept,
                                                 ontology) or \
            other.obj.concept.refinement_of(self.obj.concept, ontology)
        subj_opp = self.subj.concept.is_opposite(other.subj.concept,
                                                 ontology)
        obj_opp = self.obj.concept.is_opposite(other.obj.concept, ontology)
        sp = self.overall_polarity()
        op = other.overall_polarity()

        # If all entities are "compatible" or mutually opposites
        # but the polarities are explicitly different then this is
        # a contradiction
        if (subj_ref and obj_ref) or (subj_opp and obj_opp):
            if sp is not None and op is not None and sp != op:
                return True
        # If one entity is the opposite and the other compatible and the
        # polarities are the same then this is a contradiction
        if (subj_ref and obj_opp) or (subj_opp and obj_ref):
            if sp is not None and op is not None and sp == op:
                return True

        return False

    def overall_polarity(self):
        # Set p1 and p2 to None / 1 / -1 depending on polarity
        p1 = self.subj.delta.polarity
        p2 = self.obj.delta.polarity
        if p1 is None and p2 is None:
            pol = None
        elif p2 is None:
            pol = p1
        elif p1 is None:
            pol = p2
        else:
            pol = p1 * p2
        return pol

    def flip_polarity(self, agent_idx):
        if agent_idx == 0:
            self.subj.flip_polarity()
        elif agent_idx == 1:
            self.obj.flip_polarity()

    def polarity_count(self):
        return ((1 if self.subj.delta.polarity is not None else 0) +
                (1 if self.obj.delta.polarity is not None else 0))

    def to_json(self, use_sbo=False, matches_fun=None):
        generic = super(Influence, self).to_json(use_sbo, matches_fun)
        json_dict = _o(type=generic['type'],
                       subj=self.subj.to_json(with_evidence=False,
                                              matches_fun=matches_fun),
                       obj=self.obj.to_json(with_evidence=False,
                                            matches_fun=matches_fun))
        json_dict.update(generic)
        json_dict['matches_hash'] = int(json_dict['matches_hash'])
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        subj = json_dict.get('subj')
        obj = json_dict.get('obj')
        if subj:
            subj = Statement._from_json(subj)
        if obj:
            obj = Statement._from_json(obj)
        stmt = cls(subj, obj)
        return stmt

    def __repr__(self):
        if sys.version_info[0] >= 3:
            return self.__str__()
        else:
            return self.__str__().encode('utf-8')

    def __str__(self):
        s = "%s(%s, %s)" % (type(self).__name__, str(self.subj),
                            str(self.obj))
        return s


class Association(Complex):
    """A set of events associated with each other without causal relationship.

    Parameters
    ----------
    members : list of :py:class:Event
        A list of events associated with each other.
    evidence : None or :py:class:`Evidence` or list of :py:class:`Evidence`
        Evidence objects in support of the modification.
    """
    _agent_order = ['members']

    def __init__(self, members, evidence=None):
        if len(members) != 2:
            raise ValueError('Association Statement can only have 2 members, '
                             '%d were given.' % len(members))
        super().__init__(members, evidence)

    def matches_key(self):
        key = (stmt_type(self, True),
               tuple(m.matches_key() for m in self.sorted_members()),
               self.polarity_count(),
               self.overall_polarity()
               )
        return mk_str(key)

    def overall_polarity(self):
        p1 = self.members[0].delta.polarity
        p2 = self.members[1].delta.polarity
        if p1 is None and p2 is None:
            pol = None
        elif p2 is None:
            pol = p1
        elif p1 is None:
            pol = p2
        else:
            pol = p1 * p2
        return pol

    def polarity_count(self):
        return sum(
            1 if m.delta.polarity is not None else 0 for m in self.members)

    def agent_list(self, deep_sorted=False):
        members = self.members if not deep_sorted else \
            sorted_agents(self.members)
        return [m.concept for m in members]

    def refinement_of(self, other, ontology, entities_refined=False):
        members_refinement = \
            super().refinement_of(other, ontology,
                                  entities_refined=entities_refined)
        op = other.overall_polarity()
        sp = self.overall_polarity()
        if self.polarity_count() < other.polarity_count():
            delta_refinement = False
        elif sp is not None and op is None:
            delta_refinement = True
        else:
            delta_refinement = (op == sp)
        return (members_refinement and delta_refinement)

    def equals(self, other):

        def members_equal(a1, a2):
            if len(a1.members) == len(a2.members):
                for m1, m2 in zip(a1.members, a2.members):
                    if (m1 is None and m2 is not None) or \
                            (m1 is not None and m2 is None):
                        return False
                    if m1 is not None and m2 is not None and not m1.equals(m2):
                        return False
            else:
                return False
            return True

        equals = super().equals(other) and members_equal(self, other)
        return equals

    def contradicts(self, other, ontology):
        if stmt_type(self) != stmt_type(other):
            return False

        def match_members(self_members, other_members):
            rel_types = {'refinement_of': 0, 'is_opposite': 0}
            G = networkx.Graph()
            for (self_idx, self_member), (other_idx, other_member) in \
                itertools.product(enumerate(self_members),
                                  enumerate(other_members)):
                if ('S%d' % self_idx) in G or ('O%d' % other_idx) in G:
                    continue
                if self_member.concept.refinement_of(other_member.concept,
                                                     ontology):
                    G.add_edge('S%d' % self_idx, 'O%d' % other_idx)
                    rel_types['refinement_of'] += 1
                elif self_member.concept.is_opposite(other_member.concept,
                                                     ontology):
                    G.add_edge('S%d' % self_idx, 'O%d' % other_idx)
                    rel_types['is_opposite'] += 1
            return rel_types

        rel_types = match_members(self.members, other.members)
        sp = self.overall_polarity()
        op = other.overall_polarity()

        # If all entities are "compatible" or mutually opposites
        # but the polarities are explicitly different then this is
        # a contradiction
        if set(rel_types.values()) == set([0, 2]):
            if sp is not None and op is not None and sp != op:
                return True
        # If one entity is the opposite and the other compatible and the
        # polarities are the same then this is a contradiction
        if set(rel_types.values()) == set([1, 1]):
            if sp is not None and op is not None and sp == op:
                return True

        return False

    def flip_polarity(self, agent_idx):
        self.members[agent_idx].flip_polarity()

    def to_json(self, use_sbo=False, matches_fun=None):
        # Get generic from two inheritance levels above - from Statement class
        generic = super(Complex, self).to_json(use_sbo, matches_fun)
        json_dict = _o(type=generic['type'])
        members = [m.to_json(with_evidence=False, matches_fun=matches_fun)
                   for m in self.members]
        json_dict['members'] = members
        json_dict.update(generic)
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        members = json_dict.get('members')
        members = [Statement._from_json(m) for m in members]
        stmt = cls(members)
        return stmt


class Conversion(Statement):
    """Conversion of molecular species mediated by a controller protein.

    Parameters
    ----------
    subj : :py:class:`indra.statement.Agent`
        The protein mediating the conversion.
    obj_from : list of :py:class:`indra.statement.Agent`
        The list of molecular species being consumed by the conversion.
    obj_to : list of :py:class:`indra.statement.Agent`
        The list of molecular species being created by the conversion.
    evidence : None or :py:class:`Evidence` or list of :py:class:`Evidence`
        Evidence objects in support of the synthesis statement.
    """

    _agent_order = ['subj', 'obj_from', 'obj_to']

    def __init__(self, subj, obj_from=None, obj_to=None, evidence=None):
        super(Conversion, self).__init__(evidence=evidence)
        self.subj = subj
        self.obj_from = obj_from if obj_from is not None else []
        if isinstance(obj_from, Agent):
            self.obj_from = [obj_from]
        self.obj_to = obj_to if obj_to is not None else []
        if isinstance(obj_to, Agent):
            self.obj_to = [obj_to]

    def matches_key(self):
        keys = [stmt_type(self, True)]
        keys += [self.subj.matches_key() if self.subj else None]
        keys += [agent.matches_key() for agent in sorted_agents(self.obj_to)]
        keys += [agent.matches_key() for agent in sorted_agents(self.obj_from)]
        return mk_str(keys)

    def set_agent_list(self, agent_list):
        num_obj_from = len(self.obj_from)
        num_obj_to = len(self.obj_to)
        if len(agent_list) != 1 + num_obj_from + num_obj_to:
            raise Exception('Conversion agent number must be preserved '
                            'when setting agent list.')
        self.subj = agent_list[0]
        self.obj_from = agent_list[1:num_obj_from+1]
        self.obj_to = agent_list[num_obj_from+1:]

    def to_json(self, use_sbo=False, matches_fun=None):
        generic = super(Conversion, self).to_json(use_sbo, matches_fun)
        json_dict = _o(type=generic['type'])
        if self.subj is not None:
            json_dict['subj'] = self.subj.to_json()
            if use_sbo:
                json_dict['subj']['sbo'] = \
                    'https://identifiers.org/SBO:0000013'  # catalyst
        json_dict['obj_from'] = [o.to_json() for o in self.obj_from]
        if use_sbo:
            for of in json_dict['obj_from']:
                of['sbo'] = \
                    'https://identifiers.org/SBO:0000010'  # reactant
        json_dict['obj_to'] = [o.to_json() for o in self.obj_to]
        if use_sbo:
            for ot in json_dict['obj_to']:
                ot['sbo'] = \
                    'https://identifiers.org/SBO:0000011'  # product
        json_dict.update(generic)
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        subj = json_dict.get('subj')
        obj_from = json_dict.get('obj_from')
        obj_to = json_dict.get('obj_to')
        if subj:
            subj = Agent._from_json(subj)
        if obj_from:
            obj_from = [Agent._from_json(o) for o in obj_from]
        if obj_to:
            obj_to = [Agent._from_json(o) for o in obj_to]
        stmt = cls(subj, obj_from, obj_to)
        return stmt

    def refinement_of(self, other, ontology, entities_refined=False):
        # Make sure the statement types match
        if stmt_type(self) != stmt_type(other):
            return False

        if self.subj is None and other.subj is None:
            subj_refinement = True
        elif self.subj is None and other.subj is not None:
            subj_refinement = False
        elif self.subj is not None and other.subj is None:
            subj_refinement = True
        else:
            subj_refinement = self.subj.refinement_of(other.subj, ontology)

        def refinement_agents(lst1, lst2):
            if len(lst1) != len(lst2):
                return False
            # Check that every agent in other is refined in self, but
            # only once!
            self_match_indices = set([])
            for other_agent in lst2:
                for self_agent_ix, self_agent in enumerate(lst1):
                    if self_agent_ix in self_match_indices:
                        continue
                    if self_agent.refinement_of(other_agent, ontology):
                        self_match_indices.add(self_agent_ix)
                        break
            if len(self_match_indices) != len(lst2):
                return False
            return True

        obj_from_refinement = refinement_agents(self.obj_from, other.obj_from)
        obj_to_refinement = refinement_agents(self.obj_to, other.obj_to)

        return (subj_refinement and obj_from_refinement and obj_to_refinement)

    def equals(self, other):
        matches = super(Conversion, self).equals(other)
        return matches

    def __str__(self):
        s = ("%s(%s, %s, %s)" % (type(self).__name__, self.subj, self.obj_from,
                                 self.obj_to))
        return s


class Event(Statement):
    """An event representing the change of a Concept.

    Attributes
    ----------
    concept : indra.statements.concept.Concept
        The concept over which the event is defined.
    delta : indra.statements.delta.Delta
        Represents a change in the concept, with a polarity
        and an adjectives entry.
    context : indra.statements.context.Context
        The context associated with the event.
    """
    _agent_order = ['concept']

    def __init__(self, concept, delta=None, context=None, evidence=None,
                 supports=None, supported_by=None):
        super().__init__(evidence, supports, supported_by)
        self.concept = concept
        self.delta = delta if delta else QualitativeDelta(
                polarity=None, adjectives=None)
        self.context = context

    def matches_key(self):
        mk = (self.concept.matches_key(),)
        return str(mk)

    def refinement_of(self, other, ontology, entities_refined=False,
                      ignore_polarity=False):
        concept_ref = \
            self.concept.refinement_of(other.concept, ontology,
                                       entities_refined=entities_refined)
        if ignore_polarity:
            pol_ref = True
        else:
            pol_ref = (self.delta.polarity and not other.delta.polarity) or \
                self.delta.polarity == other.delta.polarity
        return concept_ref and pol_ref

    def equals(self, other):
        return self.concept.equals(other.concept) and \
            self.delta.equals(other.delta)

    def to_json(self, with_evidence=True, use_sbo=False, matches_fun=None):
        generic = super(Event, self).to_json(False, matches_fun)
        json_dict = _o(type=generic['type'],
                       concept=self.concept.to_json(),
                       delta=self.delta.to_json())
        if self.context:
            json_dict['context'] = self.context.to_json()
        if with_evidence and 'evidence' in generic:
            json_dict['evidence'] = generic['evidence']
        json_dict.update(generic)
        json_dict['matches_hash'] = int(json_dict['matches_hash'])
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        concept = Concept._from_json(json_dict['concept'])
        delta = QualitativeDelta.from_json(json_dict['delta'])
        context_entry = json_dict.get('context')
        if context_entry:
            context = Context.from_json(json_dict['context'])
        else:
            context = None
        stmt = cls(concept, delta=delta,
                   context=context)
        return stmt

    def flip_polarity(self, agent_idx=None):
        # If we have an explicit polarity, flip it, otherwise do nothing
        if self.delta.polarity == 1:
            self.delta.polarity = -1
        elif self.delta.polarity == -1:
            self.delta.polarity = 1

    def __str__(self):
        return '%s(%s)' % (type(self).__name__, self.concept.name)


class Migration(Event):
    """A special class of Event representing Migration."""
    _agent_order = ['concept']

    def __init__(self, concept, delta=None, context=None, evidence=None,
                 supports=None, supported_by=None):
        self.delta = delta if delta else QuantitativeState(entity='person')
        super().__init__(concept, self.delta, context, evidence, supports,
                         supported_by)

    @classmethod
    def _from_json(cls, json_dict):
        concept = Concept._from_json(json_dict['concept'])
        delta = QuantitativeState.from_json(json_dict['delta'])
        context_entry = json_dict.get('context')
        if context_entry:
            context = Context.from_json(json_dict['context'])
        else:
            context = None
        stmt = cls(concept, delta=delta, context=context)
        return stmt


class Unresolved(Statement):
    """A special statement type used in support when a uuid can't be resolved.

    When using the `stmts_from_json` method, it is sometimes not possible to
    resolve the uuid found in `support` and `supported_by` in the json
    representation of an indra statement. When this happens, this class is used
    as a place-holder, carrying only the uuid of the statement.
    """

    def __init__(self, uuid_str=None, shallow_hash=None, full_hash=None):
        super(Unresolved, self).__init__()
        self.uuid = uuid_str
        self._shallow_hash = shallow_hash
        self._full_hash = full_hash
        assert self.uuid or self._shallow_hash or self._full_hash,\
            "Some identifying information must be given."

    def __str__(self):
        if self.uuid:
            return "%s(uuid=%s)" % (type(self).__name__, self.uuid)
        elif self._shallow_hash:
            return "%s(shallow_hash=%s)" % (type(self).__name__,
                                            self._shallow_hash)
        else:
            return "%s(full_hash=%s)" % (type(self).__name__,
                                         self._full_hash)


# Mapping between modification type strings and subclasses of Modification
modtype_to_modclass = {str(cls.__name__.lower()): cls for cls in
                       AddModification.__subclasses__() +
                       RemoveModification.__subclasses__()}
# Add modification as a generic type
modtype_to_modclass['modification'] = Modification

modclass_to_modtype = {cls: str(cls.__name__.lower()) for cls in
                       AddModification.__subclasses__() +
                       RemoveModification.__subclasses__()}
# Add modification as a generic type
modclass_to_modtype[Modification] = 'modification'
modclass_to_modtype[Autophosphorylation] = 'phosphorylation'
modclass_to_modtype[Transphosphorylation] = 'phosphorylation'

# These are the modification types that are valid in ModConditions
modtype_conditions = [modclass_to_modtype[mt] for mt in
                      AddModification.__subclasses__()]
modtype_conditions.append('modification')


def _get_mod_inverse_maps():
    modtype_to_inverse = {}
    modclass_to_inverse = {}
    for cls in AddModification.__subclasses__():
        modtype = modclass_to_modtype[cls]
        modtype_inv = 'de' + modtype
        cls_inv = modtype_to_modclass[modtype_inv]
        modtype_to_inverse[modtype] = modtype_inv
        modtype_to_inverse[modtype_inv] = modtype
        modclass_to_inverse[cls] = cls_inv
        modclass_to_inverse[cls_inv] = cls
    return modtype_to_inverse, modclass_to_inverse


modtype_to_inverse, modclass_to_inverse = _get_mod_inverse_maps()


from .io import *
from .agent import *


stmt_sbo_map = {
    'acetylation': '0000215',
    'glycosylation': '0000217',
    'hydroxylation': '0000233',
    'methylation': '0000214',
    'myristoylation': '0000219',
    'palmitoylation': '0000218',
    'phosphorylation': '0000216',
    'farnesylation': '0000222',
    'geranylgeranylation': '0000223',
    'ubiquitination': '0000224',
    'dephosphorylation': '0000330',
    'addmodification': '0000210',  # addition of a chemical group
    'removemodification': '0000211',  # removal of a chemical group
    'modification': '0000182',  # conversion
    'conversion': '0000182',  # conversion
    'autophosphorylation': '0000216',  # phosphorylation
    'transphosphorylation': '0000216',  # phosphorylation
    'decreaseamount': '0000179',  # degradation
    'increaseamount': '0000183',  # transcription
    'complex': '0000526',  # protein complex formation
    'translocation': '0000185',  # transport reaction
    'regulateactivity': '0000182',  # conversion
    'activeform': '0000412',  # biological activity
    'rasgef': '0000172',  # catalysis
    'rasgap': '0000172',  # catalysis
    'statement': '0000231'  # occuring entity representation
    }


def sorted_agents(agent_list):
    return sorted(agent_list, key=lambda ag: ag.matches_key())


def get_all_descendants(parent):
    """Get all the descendants of a parent class, recursively."""
    children = parent.__subclasses__()
    descendants = children[:]
    for child in children:
        descendants += get_all_descendants(child)
    return descendants


# In the future, when hierarchy is no longer determined by sub-classing, this
# function should be altered to account for the change.
def get_type_hierarchy(s):
    """Get the sequence of parents from `s` to Statement.

    Parameters
    ----------
    s : a class or instance of a child of Statement
        For example the statement `Phosphorylation(MEK(), ERK())` or just the
        class `Phosphorylation`.

    Returns
    -------
    parent_list : list[types]
        A list of the types leading up to Statement.

    Examples
    --------
        >> s = Phosphorylation(MAPK1(), Elk1())
        >> get_type_hierarchy(s)
        [Phosphorylation, AddModification, Modification, Statement]
        >> get_type_hierarchy(AddModification)
        [AddModification, Modification, Statement]
    """
    tp = type(s) if not isinstance(s, type) else s
    p_list = [tp]
    for p in tp.__bases__:
        if p is not Statement:
            p_list.extend(get_type_hierarchy(p))
        else:
            p_list.append(p)
    return p_list


class NotAStatementName(Exception):
    pass


def get_statement_by_name(stmt_name):
    """Get a statement class given the name of the statement class."""
    stmt_classes = get_all_descendants(Statement)
    for stmt_class in stmt_classes:
        if stmt_class.__name__.lower() == stmt_name.lower():
            return stmt_class
    raise NotAStatementName('\"%s\" is not recognized as a statement type!'
                            % stmt_name)


def make_statement_camel(stmt_name):
    """Makes a statement name match the case of the corresponding statement."""
    return get_statement_by_name(stmt_name).__name__


def get_unresolved_support_uuids(stmts):
    """Get uuids unresolved in support from stmts from stmts_from_json."""
    return {s.uuid for stmt in stmts for s in stmt.supports + stmt.supported_by
            if isinstance(s, Unresolved)}


def stmt_type(obj, mk=True):
    """Return standardized, backwards compatible object type String.

    This is a temporary solution to make sure type comparisons and
    matches keys of Statements and related classes are backwards
    compatible.
    """
    if isinstance(obj, Statement) and mk:
        return type(obj)
    else:
        return type(obj).__name__


def mk_str(mk):
    """Replace class path for backwards compatibility of matches keys."""
    return str(mk).replace('indra.statements.statements', 'indra.statements')
