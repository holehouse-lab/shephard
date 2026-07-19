"""
SHEPHARD comprehensive test suite : Domain

Detailed tests covering the user-facing Domain API:

  * construction & validation (start>end, out-of-range, attribute typing)
  * properties (start/end/protein/sequence/domain_type/domain_name)
  * update_domain_name
  * attributes (add/get/remove/safe)
  * inside_domain / domain_overlap
  * site accessors (sites, site_positions, site(), get_sites_by_type)
  * track accessors (get_track_values / get_track_symbols + error paths)
  * dunders (__len__, __repr__)

Holehouse Lab - Washington University in St. Louis
"""

import pytest

from shephard.proteome import Proteome
from shephard.protein import Protein
from shephard.domain import Domain
from shephard.exceptions import DomainException

SEQ = 'ACDEFGHIKLMNPQRSTVWY'   # 20 residues


@pytest.fixture
def prot():
    P = Proteome([{'sequence': SEQ, 'name': 'demo',
                   'unique_ID': 'P1', 'attributes': {}}])
    return P.protein('P1')


@pytest.fixture
def domain(prot):
    prot.add_domain(5, 12, 'IDR', attributes={'score': 0.8})
    return prot.get_domains_by_type('IDR')[0]


# ---------------------------------------------------------------------------
# construction / validation
# ---------------------------------------------------------------------------
class TestDomainConstruction:

    def test_start_greater_than_end_raises(self, prot):
        with pytest.raises(DomainException):
            prot.add_domain(10, 5, 'bad')

    def test_out_of_range_raises(self, prot):
        with pytest.raises(Exception):
            prot.add_domain(1, 99, 'bad')
        with pytest.raises(Exception):
            prot.add_domain(0, 5, 'bad')

    def test_bad_attributes_raises(self, prot):
        with pytest.raises(DomainException):
            prot.add_domain(1, 5, 'd', attributes='not-a-dict')

    def test_single_residue_domain(self, prot):
        prot.add_domain(7, 7, 'point')
        d = prot.get_domains_by_type('point')[0]
        assert d.start == d.end == 7
        assert len(d) == 1


# ---------------------------------------------------------------------------
# properties
# ---------------------------------------------------------------------------
class TestDomainProperties:

    def test_basic_properties(self, domain):
        assert domain.start == 5
        assert domain.end == 12
        assert domain.domain_type == 'IDR'
        assert isinstance(domain.protein, Protein)
        assert domain.domain_name == 'IDR_5_12'

    def test_sequence_property(self, domain):
        # inclusive 1-based slice -> SEQ[4:12]
        assert domain.sequence == SEQ[4:12]
        assert len(domain.sequence) == 8

    def test_update_domain_name(self, domain):
        domain.update_domain_name('renamed')
        assert domain.domain_name == 'renamed'


# ---------------------------------------------------------------------------
# attributes
# ---------------------------------------------------------------------------
class TestDomainAttributes:

    def test_constructor_attribute(self, domain):
        assert domain.attribute('score') == 0.8
        assert 'score' in domain.attributes

    def test_add_get_duplicate_safe(self, domain):
        domain.add_attribute('k', 1)
        with pytest.raises(DomainException):
            domain.add_attribute('k', 2)
        domain.add_attribute('k', 2, safe=False)
        assert domain.attribute('k') == 2

    def test_missing_attribute_safe(self, domain):
        with pytest.raises(DomainException):
            domain.attribute('absent')
        assert domain.attribute('absent', safe=False) is None

    def test_remove_attribute(self, domain):
        domain.add_attribute('k', 1)
        domain.remove_attribute('k')
        assert 'k' not in domain.attributes
        with pytest.raises(Exception):
            domain.remove_attribute('k')
        domain.remove_attribute('k', safe=False)


# ---------------------------------------------------------------------------
# geometry helpers
# ---------------------------------------------------------------------------
class TestDomainGeometry:

    def test_inside_domain(self, domain):
        assert domain.inside_domain(5) is True
        assert domain.inside_domain(12) is True
        assert domain.inside_domain(4) is False
        assert domain.inside_domain(13) is False

    def test_domain_overlap_true(self, prot):
        prot.add_domain(1, 10, 'A')
        prot.add_domain(5, 15, 'B')
        a = prot.get_domains_by_type('A')[0]
        b = prot.get_domains_by_type('B')[0]
        assert a.domain_overlap(b) is True

    def test_domain_overlap_false(self, prot):
        prot.add_domain(1, 5, 'A')
        prot.add_domain(10, 15, 'B')
        a = prot.get_domains_by_type('A')[0]
        b = prot.get_domains_by_type('B')[0]
        assert a.domain_overlap(b) is False


# ---------------------------------------------------------------------------
# site accessors on a domain
# ---------------------------------------------------------------------------
class TestDomainSites:

    def test_sites_within_domain(self, prot):
        prot.add_domain(5, 12, 'IDR')
        d = prot.get_domains_by_type('IDR')[0]
        prot.add_site(6, 'phos', value=1)
        prot.add_site(10, 'acet', value=1)
        prot.add_site(1, 'phos', value=1)        # outside domain
        assert len(d.sites) == 2
        assert sorted(d.site_positions) == [6, 10]

    def test_site_lookup(self, prot):
        prot.add_domain(5, 12, 'IDR')
        d = prot.get_domains_by_type('IDR')[0]
        prot.add_site(6, 'phos', value=1)
        assert len(d.site(6)) == 1
        # in-range but no site -> DomainException (regression: not KeyError)
        with pytest.raises(DomainException):
            d.site(7)
        # outside the domain boundary -> DomainException
        with pytest.raises(DomainException):
            d.site(1)

    def test_get_sites_by_type(self, prot):
        prot.add_domain(5, 12, 'IDR')
        d = prot.get_domains_by_type('IDR')[0]
        prot.add_site(6, 'phos', value=1)
        prot.add_site(8, 'phos', value=1)
        prot.add_site(9, 'acet', value=1)
        res = d.get_sites_by_type('phos')
        flat = [s for v in res.values() for s in v]
        assert len(flat) == 2


# ---------------------------------------------------------------------------
# track accessors on a domain
# ---------------------------------------------------------------------------
class TestDomainTracks:

    def test_get_track_values(self, prot):
        prot.add_domain(5, 12, 'IDR')
        d = prot.get_domains_by_type('IDR')[0]
        prot.add_track('v', values=list(range(20)))
        vals = d.get_track_values('v')
        # values_region is 1-based inclusive over residues 5..12
        assert vals == [float(i) for i in range(4, 12)]

    def test_get_track_symbols(self, prot):
        prot.add_domain(5, 12, 'IDR')
        d = prot.get_domains_by_type('IDR')[0]
        prot.add_track('s', symbols=list('ABCDEFGHIJKLMNOPQRST'))
        syms = d.get_track_symbols('s')
        assert syms == list('ABCDEFGHIJKLMNOPQRST')[4:12]

    def test_get_track_values_missing_safe(self, prot):
        prot.add_domain(5, 12, 'IDR')
        d = prot.get_domains_by_type('IDR')[0]
        with pytest.raises(Exception):
            d.get_track_values('absent')
        assert d.get_track_values('absent', safe=False) is None

    def test_get_track_values_on_symbols_track_raises(self, prot):
        prot.add_domain(5, 12, 'IDR')
        d = prot.get_domains_by_type('IDR')[0]
        prot.add_track('s', symbols=list('X' * 20))
        with pytest.raises(DomainException):
            d.get_track_values('s')


# ---------------------------------------------------------------------------
# dunders
# ---------------------------------------------------------------------------
class TestDomainDunders:

    def test_len(self, domain):
        assert len(domain) == 8           # residues 5..12 inclusive

    def test_repr(self, domain):
        r = repr(domain)
        assert 'Domain' in r and 'IDR' in r
