import shephard
from shephard.tools import attribute_tools
from shephard import Proteome


def test_cast_basic_protein():
    
    p = Proteome([])
    p.add_protein('PRTEINSEQWENCE', 'prot a', 'test_1', attributes={'value_1':'1', 'value_2':'3'})

    # can cast all attributes 
    assert p.protein('test_1').attribute('value_1') == '1'
    attribute_tools.cast_attributes(p.protein('test_1'))
    assert p.protein('test_1').attribute('value_1') == 1.0

    attribute_tools.cast_attributes(p.protein('test_1'),cast_type=int)
    assert p.protein('test_1').attribute('value_1') == 1

    attribute_tools.cast_attributes(p.protein('test_1'),cast_type=str)
    assert p.protein('test_1').attribute('value_1') == '1'


def test_cast_basic_protein_with_exclude():
    
    p = Proteome([])
    p.add_protein('PRTEINSEQWENCE', 'prot a', 'test_1', attributes={'value_1':'1', 'value_2':'3', 'string_1':'hello'})

    # can cast all attributes 
    assert p.protein('test_1').attribute('value_1') == '1'
    attribute_tools.cast_attributes(p.protein('test_1'), exclude=['string_1'])
    assert p.protein('test_1').attribute('value_1') == 1.0
    assert p.protein('test_1').attribute('string_1') == 'hello'
    

    attribute_tools.cast_attributes(p.protein('test_1'), cast_type=int, exclude=['string_1'])
    assert p.protein('test_1').attribute('value_1') == 1
    assert p.protein('test_1').attribute('string_1') == 'hello'

    attribute_tools.cast_attributes(p.protein('test_1'),cast_type=str, exclude=['string_1'])
    assert p.protein('test_1').attribute('value_1') == '1'
    assert p.protein('test_1').attribute('string_1') == 'hello'

    # don't need eclude if everything is already a string
    attribute_tools.cast_attributes(p.protein('test_1'),cast_type=str)
    assert p.protein('test_1').attribute('value_1') == '1'
    assert p.protein('test_1').attribute('string_1') == 'hello'


def test_cast_basic_protein_with_include():
    
    p = Proteome([])
    p.add_protein('PRTEINSEQWENCE', 'prot a', 'test_1', attributes={'value_1':'1', 'value_2':'3', 'string_1':'hello'})

    # can cast all attributes 
    assert p.protein('test_1').attribute('value_1') == '1'
    attribute_tools.cast_attributes(p.protein('test_1'), include=['value_1'])
    assert p.protein('test_1').attribute('value_1') == 1.0
    assert p.protein('test_1').attribute('value_2') == '3'
    assert p.protein('test_1').attribute('string_1') == 'hello'

    attribute_tools.cast_attributes(p.protein('test_1'), cast_type=int, include=['value_1'])
    assert p.protein('test_1').attribute('value_1') == 1
    assert p.protein('test_1').attribute('value_2') == '3'
    assert p.protein('test_1').attribute('string_1') == 'hello'

    attribute_tools.cast_attributes(p.protein('test_1'), cast_type=int, include=['value_1', 'value_2'])
    assert p.protein('test_1').attribute('value_1') == 1
    assert p.protein('test_1').attribute('value_2') == 3
    assert p.protein('test_1').attribute('string_1') == 'hello'
    
    attribute_tools.cast_attributes(p.protein('test_1'),cast_type=str, include=['value_1'])
    assert p.protein('test_1').attribute('value_1') == '1'
    assert p.protein('test_1').attribute('value_2') == 3
    assert p.protein('test_1').attribute('string_1') == 'hello'


## to do - same tests for domains, sites, proteomes...
