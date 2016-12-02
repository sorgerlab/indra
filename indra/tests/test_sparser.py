from indra import sparser

xml_str = '''
 <interpretation>
 <sentence-text>MEK1 phosphorylates ERK1</sentence-text>
 <sem>
 <ref category="phosphorylate">
 <var name="agent">
 <ref category="protein" name="MP2K1_HUMAN" uid="UP:MP2K1_HUMAN"/>
 </var>
 <var name="substrate">
 <ref category="protein" name="MK03_HUMAN" uid="UP:MK03_HUMAN"/>
 </var>
 <var name="present"><ref category="present"></ref></var>
 </ref>
 </sem>
</interpretation>
'''

def test_invalid_xml():
    sp = sparser.process_xml('xyz')
    assert(sp is None)

def test_phosphorylation():
    sp = sparser.process_xml(xml_str)
    print(sp.tree)
    assert(len(sp.statements) == 1)
