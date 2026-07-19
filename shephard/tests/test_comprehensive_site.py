"""
SHEPHARD comprehensive test suite : Site

Detailed tests covering the user-facing Site API:

  * properties (position/residue/protein/site_type/symbol/value)
  * value/symbol casting & None handling
  * update_site_value / update_site_symbol
  * attributes (add/get/remove/safe)
  * get_local_sequence_context
  * get_domains (with offset)
  * get_track_values / get_track_value / get_track_symbols / get_track_symbol
  * dunders

Holehouse Lab - Washington University in St. Louis
"""

import pytest

from shephard.proteome import Proteome
from shephard.protein import Protein
from shephard.exceptions import SiteException

SEQ = 'ACDEFGHIKLMNPQRSTVWY'


@pytest.fixture
def prot():
    P = Proteome([{'sequence': SEQ, 'name': 'demo',
                   'unique_ID': 'P1', 'attributes': {}}])
    return P.protein('P1')


@pytest.fixture
def site(prot):
    prot.add_site(10, 'phosphosite', symbol='P', value=0.75,
                  attributes={'enzyme': 'CK2'})
    return prot.site(10)[0]


# ---------------------------------------------------------------------------
# properties / casting
# ---------------------------------------------------------------------------
class TestSiteProperties:

    def test_basic_properties(self, site):
        assert site.position == 10
        assert site.site_type == 'phosphosite'
        assert site.symbol == 'P'
        assert site.value == 0.75
        assert isinstance(site.protein, Protein)

    def test_residue(self, site):
        assert site.residue == SEQ[9]    # position 10 -> index 9

    def test_value_cast_to_float(self, prot):
        prot.add_site(2, 't', value='3')
        s = prot.site(2)[0]
        assert s.value == 3.0 and isinstance(s.value, float)

    def test_symbol_cast_to_str(self, prot):
        prot.add_site(2, 't', symbol=7)
        assert prot.site(2)[0].symbol == '7'

    def test_none_value_and_symbol(self, prot):
        prot.add_site(2, 't')
        s = prot.site(2)[0]
        assert s.value is None
        assert s.symbol is None

    def test_bad_attributes_raises(self, prot):
        with pytest.raises(SiteException):
            prot.add_site(2, 't', attributes='not-a-dict')


# ---------------------------------------------------------------------------
# update value / symbol
# ---------------------------------------------------------------------------
class TestSiteUpdate:

    def test_update_value(self, site):
        site.update_site_value(9.0)
        assert site.value == 9.0
        site.update_site_value(None)
        assert site.value is None

    def test_update_value_casts(self, site):
        site.update_site_value('2.5')
        assert site.value == 2.5

    def test_update_symbol(self, site):
        site.update_site_symbol('Q')
        assert site.symbol == 'Q'
        site.update_site_symbol(None)
        assert site.symbol is None


# ---------------------------------------------------------------------------
# attributes
# ---------------------------------------------------------------------------
class TestSiteAttributes:

    def test_constructor_attribute(self, site):
        assert site.attribute('enzyme') == 'CK2'
        assert 'enzyme' in site.attributes

    def test_add_duplicate_safe(self, site):
        site.add_attribute('k', 1)
        with pytest.raises(SiteException):
            site.add_attribute('k', 2)
        site.add_attribute('k', 2, safe=False)
        assert site.attribute('k') == 2

    def test_missing_safe(self, site):
        with pytest.raises(SiteException):
            site.attribute('absent')
        assert site.attribute('absent', safe=False) is None

    def test_remove(self, site):
        site.add_attribute('k', 1)
        site.remove_attribute('k')
        assert 'k' not in site.attributes
        with pytest.raises(Exception):
            site.remove_attribute('k')
        site.remove_attribute('k', safe=False)


# ---------------------------------------------------------------------------
# sequence context
# ---------------------------------------------------------------------------
class TestSiteSequenceContext:

    def test_local_context_default_offset(self, site):
        # position 10, offset 5 -> SEQ[4:15]
        assert site.get_local_sequence_context(offset=5) == SEQ[4:15]

    def test_local_context_truncates(self, prot):
        prot.add_site(1, 't', value=1)
        s = prot.site(1)[0]
        assert s.get_local_sequence_context(offset=5) == SEQ[0:6]


# ---------------------------------------------------------------------------
# site -> domains
# ---------------------------------------------------------------------------
class TestSiteDomains:

    def test_get_domains(self, prot):
        prot.add_domain(5, 15, 'IDR')
        prot.add_site(10, 't', value=1)
        s = prot.site(10)[0]
        doms = s.get_domains()
        assert len(doms) == 1
        assert doms[0].domain_type == 'IDR'

    def test_get_domains_offset_extends_reach(self, prot):
        prot.add_domain(5, 9, 'IDR')
        prot.add_site(11, 't', value=1)
        s = prot.site(11)[0]
        assert len(s.get_domains(offset=0)) == 0
        assert len(s.get_domains(offset=3)) == 1


# ---------------------------------------------------------------------------
# site -> tracks
# ---------------------------------------------------------------------------
class TestSiteTracks:

    def test_get_track_value(self, prot):
        prot.add_track('v', values=list(range(20)))
        prot.add_site(10, 't', value=1)
        s = prot.site(10)[0]
        assert s.get_track_value('v') == 9.0     # position 10 -> values[9]

    def test_get_track_values_offset(self, prot):
        prot.add_track('v', values=list(range(20)))
        prot.add_site(10, 't', value=1)
        s = prot.site(10)[0]
        # offset 2 around position 10 -> residues 8..12 -> values 7..11
        assert s.get_track_values('v', offset=2) == [7.0, 8.0, 9.0, 10.0, 11.0]

    def test_get_track_symbol(self, prot):
        prot.add_track('s', symbols=list('ABCDEFGHIJKLMNOPQRST'))
        prot.add_site(3, 't', value=1)
        s = prot.site(3)[0]
        assert s.get_track_symbol('s') == 'C'

    def test_get_track_symbols_offset(self, prot):
        prot.add_track('s', symbols=list('ABCDEFGHIJKLMNOPQRST'))
        prot.add_site(3, 't', value=1)
        s = prot.site(3)[0]
        assert s.get_track_symbols('s', offset=1) == ['B', 'C', 'D']

    def test_missing_track_safe(self, prot):
        prot.add_site(3, 't', value=1)
        s = prot.site(3)[0]
        with pytest.raises(Exception):
            s.get_track_values('absent')
        assert s.get_track_values('absent', safe=False) is None
        assert s.get_track_symbols('absent', safe=False) is None


# ---------------------------------------------------------------------------
# dunders
# ---------------------------------------------------------------------------
class TestSiteDunders:

    def test_repr(self, site):
        r = repr(site)
        assert 'Site' in r and 'phosphosite' in r
