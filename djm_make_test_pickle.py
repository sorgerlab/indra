from indra.tools.extract_grounding_map import TmpDebugTestClass
import pickle
from indra.statements import Agent

a1 = Agent('a1', db_refs={'TEXT': 'a1', 'number': 1})
a2 = Agent('a2', db_refs={'TEXT': 'a2', 'number': 2})
a3 = Agent('a3', db_refs={'TEXT': 'a3', 'number': 3})
a4 = Agent('a4', db_refs={'TEXT': 'a4', 'number': 4})

s1 = TmpDebugTestClass([a1, a2])
s2 = TmpDebugTestClass(a3)
s3 = TmpDebugTestClass([a4])

d = {}
d['foo'] = [s1, s2]
d['bar'] = [s3]

pickle.dump(d, open('test.p', 'wb'))
