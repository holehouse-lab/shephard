"""
SHEPHARD comprehensive test suite : Protein

Detailed tests covering the user-facing Protein API:

  * identity / sequence properties and 1-based indexing semantics
  * residue / get_sequence_region / get_sequence_context (+truncation, indices)
  * check_sequence_is_valid / convert_to_valid (copy & length-change safety)
  * protein attributes
  * track creation/retrieval (add_track, build_track*, get_track_*),
  * domain creation/retrieval (add_domain(s), autoname, get_domains_by_*)
  * site creation/retrieval (add_site, get_sites_by_*), return_list, wiggle
  * remove_track / remove_domain / remove_site
  * dunder methods

Holehouse Lab - Washington University in St. Louis
"""

import pytest

import shephard
from shephard.proteome import Proteome
from shephard.protein import Protein
from shephard.domain import Domain
from shephard.site import Site
from shephard.track import Track
from shephard.exceptions import ProteinException

SEQ = 'ACDEFGHIKLMNPQRSTVWY'   # 20 residues, all standard


def _proteome(seq=SEQ, uid='P1'):
    return Proteome([{'sequence': seq, 'name': 'demo',
                      'unique_ID': uid, 'attributes': {}}])


@pytest.fixture
def prot():
    return _proteome().protein('P1')


# ---------------------------------------------------------------------------
# identity / sequence
# ---------------------------------------------------------------------------
class TestProteinIdentityAndSequence:

    def test_identity(self, prot):
        assert prot.unique_ID == 'P1'
        assert prot.name == 'demo'
        assert isinstance(prot.proteome, Proteome)

    def test_sequence_and_len(self, prot):
        assert prot.sequence == SEQ
        assert len(prot) == 20

    def test_one_based_indexing(self, prot):
        # residue 1 is the first residue, not the leading sentinel
        assert prot.residue(1) == 'A'
        assert prot.residue(20) == 'Y'

    def test_residue_out_of_range_raises(self, prot):
        with pytest.raises(ProteinException):
            prot.residue(0)
        with pytest.raises(ProteinException):
            prot.residue(21)

    def test_get_sequence_region_inclusive(self, prot):
        assert prot.get_sequence_region(1, 3) == 'ACD'
        assert prot.get_sequence_region(18, 20) == 'VWY'
        assert prot.get_sequence_region(5, 5) == 'F'

    def test_get_sequence_region_invalid(self, prot):
        with pytest.raises(ProteinException):
            prot.get_sequence_region(0, 3)
        with pytest.raises(ProteinException):
            prot.get_sequence_region(1, 99)

    def test_get_sequence_context_centered(self, prot):
        assert prot.get_sequence_context(10, offset=2) == SEQ[7:12]

    def test_get_sequence_context_truncates_at_termini(self, prot):
        assert prot.get_sequence_context(1, offset=5) == SEQ[0:6]
        assert prot.get_sequence_context(20, offset=5) == SEQ[14:20]

    def test_get_sequence_context_return_indices(self, prot):
        s, p1, p2 = prot.get_sequence_context(10, offset=2, return_indices=True)
        assert s == SEQ[7:12]
        assert (p1, p2) == (8, 12)

    def test_check_sequence_is_valid(self):
        assert _proteome(SEQ).protein('P1').check_sequence_is_valid() is True
        assert _proteome('ACDEBXZ').protein('P1').check_sequence_is_valid() is False

    def test_convert_to_valid_inplace(self):
        p = _proteome('ACDEBUXZ').protein('P1')   # same length conversions
        p.convert_to_valid()
        assert p.sequence == 'ACDENCGQ'
        assert p.check_sequence_is_valid() is True

    def test_convert_to_valid_copy_does_not_mutate(self):
        p = _proteome('ACDEBUXZ').protein('P1')
        out = p.convert_to_valid(copy=True)
        assert out == 'ACDENCGQ'
        assert p.sequence == 'ACDEBUXZ'   # original untouched

    def test_convert_to_valid_length_change_safe(self):
        p = _proteome('ACDE*-').protein('P1')
        with pytest.raises(ProteinException):
            p.convert_to_valid()                 # would shorten -> raise
        # safe=False permits the (length-changing) conversion
        p.convert_to_valid(safe=False)
        assert p.sequence == 'ACDE'


# ---------------------------------------------------------------------------
# protein attributes
# ---------------------------------------------------------------------------
class TestProteinAttributes:

    def test_add_get(self, prot):
        prot.add_attribute('gene', 'WASL')
        assert 'gene' in prot.attributes
        assert prot.attribute('gene') == 'WASL'

    def test_duplicate_safe(self, prot):
        prot.add_attribute('k', 1)
        with pytest.raises(ProteinException):
            prot.add_attribute('k', 2)
        prot.add_attribute('k', 2, safe=False)
        assert prot.attribute('k') == 2

    def test_attribute_missing_safe(self, prot):
        with pytest.raises(ProteinException):
            prot.attribute('absent')
        assert prot.attribute('absent', safe=False) is None

    def test_remove(self, prot):
        prot.add_attribute('k', 1)
        prot.remove_attribute('k')
        assert 'k' not in prot.attributes
        with pytest.raises(ProteinException):
            prot.remove_attribute('k')
        prot.remove_attribute('k', safe=False)


# ---------------------------------------------------------------------------
# tracks
# ---------------------------------------------------------------------------
class TestProteinTracks:

    def test_add_values_track(self, prot):
        prot.add_track('v', values=list(range(20)))
        t = prot.track('v')
        assert isinstance(t, Track)
        assert t.track_type == 'values'
        assert prot.get_track_values('v') == [float(i) for i in range(20)]
        assert prot.track_names == ['v']
        assert len(prot.tracks) == 1

    def test_add_symbols_track(self, prot):
        prot.add_track('s', symbols=list('X' * 20))
        assert prot.track('s').track_type == 'symbols'
        assert prot.get_track_symbols('s') == ['X'] * 20

    def test_add_track_length_mismatch_raises(self, prot):
        with pytest.raises(Exception):
            prot.add_track('bad', values=[1, 2, 3])

    def test_add_track_duplicate_safe(self, prot):
        prot.add_track('v', values=[0.0] * 20)
        with pytest.raises(ProteinException):
            prot.add_track('v', values=[1.0] * 20)
        prot.add_track('v', values=[2.0] * 20, safe=False)
        assert prot.get_track_values('v') == [2.0] * 20

    def test_get_track_values_subregion(self, prot):
        prot.add_track('v', values=list(range(20)))
        assert prot.get_track_values('v', start=1, end=3) == [0.0, 1.0, 2.0]

    def test_track_missing_safe(self, prot):
        with pytest.raises(ProteinException):
            prot.track('absent')
        assert prot.track('absent', safe=False) is None
        assert prot.get_track_values('absent', safe=False) is None

    def test_build_track_values_from_sequence(self, prot):
        def charge(seq):
            return [1 if r in 'KR' else 0 for r in seq]
        prot.build_track_values_from_sequence('charge', charge)
        vals = prot.get_track_values('charge')
        assert vals[SEQ.index('K')] == 1.0
        assert vals[SEQ.index('A')] == 0.0

    def test_build_track_values_from_sequence_with_dict(self, prot):
        def fx(seq, d):
            return [1 if r in d['pos'] else 0 for r in seq]
        prot.build_track_values_from_sequence('c', fx, input_dictionary={'pos': 'KR'})
        assert prot.get_track_values('c')[SEQ.index('R')] == 1.0

    def test_build_track_symbols_from_sequence(self, prot):
        def classify(seq):
            return ''.join('+' if r in 'KR' else '0' for r in seq)
        prot.build_track_symbols_from_sequence('cls', classify)
        syms = prot.get_track_symbols('cls')
        assert syms[SEQ.index('K')] == '+'

    def test_build_track_generic(self, prot):
        def deffn(data):
            return {'values': [float(x) for x in data], 'symbols': None}
        prot.build_track('g', list(range(20)), deffn)
        assert prot.get_track_values('g')[5] == 5.0

    def test_remove_track(self, prot):
        prot.add_track('v', values=[0.0] * 20)
        t = prot.track('v')
        prot.remove_track(t)
        assert 'v' not in prot.track_names
        assert 'v' not in prot.proteome.unique_track_names

    def test_remove_track_none_safe(self, prot):
        prot.remove_track(None, safe=False)   # tolerated
        with pytest.raises(ProteinException):
            prot.remove_track(None, safe=True)


# ---------------------------------------------------------------------------
# domains
# ---------------------------------------------------------------------------
class TestProteinDomains:

    def test_add_domain_and_accessors(self, prot):
        prot.add_domain(1, 5, 'IDR')
        prot.add_domain(10, 15, 'FD')
        assert len(prot.domains) == 2
        # domains sorted by start
        assert [d.start for d in prot.domains] == [1, 10]
        assert sorted(prot.domain_types) == ['FD', 'IDR']
        assert len(prot.domain_names) == 2

    def test_add_domain_duplicate_safe(self, prot):
        prot.add_domain(1, 5, 'IDR')
        with pytest.raises(ProteinException):
            prot.add_domain(1, 5, 'IDR')

    def test_add_domain_autoname_allows_overlap(self, prot):
        prot.add_domain(1, 5, 'IDR')
        prot.add_domain(1, 5, 'IDR', autoname=True)
        assert len(prot.domains) == 2

    def test_add_domains_bulk(self, prot):
        prot.add_domains([{'start': 1, 'end': 4, 'domain_type': 'A'},
                          {'start': 6, 'end': 9, 'domain_type': 'B',
                           'attributes': {'x': 1}}])
        assert len(prot.domains) == 2
        b = prot.get_domains_by_type('B')[0]
        assert b.attribute('x') == 1

    def test_build_domain(self, prot):
        def deffn(_):
            return [{'start': 2, 'end': 8, 'domain_type': 'built'}]
        prot.build_domain(None, deffn)
        assert prot.get_domains_by_type('built')[0].start == 2

    def test_domain_lookup_by_name(self, prot):
        prot.add_domain(1, 5, 'IDR')
        name = prot.domain_names[0]
        assert isinstance(prot.domain(name), Domain)
        with pytest.raises(ProteinException):
            prot.domain('no-such')
        assert prot.domain('no-such', safe=False) is None

    def test_get_domains_by_position(self, prot):
        prot.add_domain(5, 10, 'D')
        assert len(prot.get_domains_by_position(7)) == 1
        assert len(prot.get_domains_by_position(1)) == 0
        assert len(prot.get_domains_by_position(1, wiggle=10)) == 1

    def test_get_domains_by_position_and_type(self, prot):
        prot.add_domain(5, 10, 'D')
        prot.add_domain(5, 10, 'E', autoname=True)
        assert len(prot.get_domains_by_position_and_type(7, 'D')) == 1

    def test_get_domains_by_range_modes(self, prot):
        prot.add_domain(50 if False else 5, 10, 'D')  # domain 5-10
        # internal: query fully inside domain
        assert prot.get_domains_by_range(6, 8, mode='internal')
        # internal: query straddling -> not found
        assert not prot.get_domains_by_range(1, 20, mode='internal')
        # overlap-strict: query fully covers domain
        assert prot.get_domains_by_range(1, 20, mode='overlap-strict')
        # overlap: partial straddle counts
        assert prot.get_domains_by_range(1, 7, mode='overlap')
        assert not prot.get_domains_by_range(1, 7, mode='overlap-strict')

    def test_get_domains_by_range_bad_mode(self, prot):
        prot.add_domain(5, 10, 'D')
        with pytest.raises(Exception):
            prot.get_domains_by_range(1, 5, mode='nonsense')

    def test_get_domains_by_range_negative_wiggle(self, prot):
        prot.add_domain(5, 10, 'D')
        with pytest.raises(ProteinException):
            prot.get_domains_by_range(5, 10, wiggle=-1)

    def test_get_domains_by_type_perfect_and_partial(self, prot):
        prot.add_domain(1, 4, 'IDR_long')
        prot.add_domain(6, 9, 'IDR_short')
        assert len(prot.get_domains_by_type('IDR_long')) == 1
        assert len(prot.get_domains_by_type('IDR', perfect_match=False)) == 2

    def test_remove_domain(self, prot):
        prot.add_domain(1, 5, 'IDR')
        d = prot.get_domains_by_type('IDR')[0]
        prot.remove_domain(d)
        assert len(prot.domains) == 0

    def test_remove_domain_none_safe(self, prot):
        prot.remove_domain(None, safe=False)
        with pytest.raises(ProteinException):
            prot.remove_domain(None, safe=True)


# ---------------------------------------------------------------------------
# sites
# ---------------------------------------------------------------------------
class TestProteinSites:

    def test_add_site_and_accessors(self, prot):
        prot.add_site(3, 'phos', symbol='P', value=0.9)
        prot.add_site(3, 'acet', symbol='A', value=1.0)
        prot.add_site(10, 'phos', value=0.5)
        assert len(prot.sites) == 3
        assert prot.site_positions == [3, 10]
        assert sorted(prot.site_types) == ['acet', 'phos']
        assert len(prot.site(3)) == 2

    def test_add_site_out_of_range_raises(self, prot):
        with pytest.raises(ProteinException):
            prot.add_site(0, 'phos')
        with pytest.raises(ProteinException):
            prot.add_site(99, 'phos')

    def test_site_value_and_symbol_casting(self, prot):
        prot.add_site(2, 'phos', symbol=5, value='1.5')
        s = prot.site(2)[0]
        assert s.value == 1.5 and isinstance(s.value, float)
        assert s.symbol == '5'

    def test_site_missing_safe(self, prot):
        with pytest.raises(ProteinException):
            prot.site(5)
        assert prot.site(5, safe=False) == []

    def test_get_sites_by_position(self, prot):
        prot.add_site(5, 'phos', value=1)
        d = prot.get_sites_by_position(5)
        assert isinstance(d, dict) and 5 in d
        lst = prot.get_sites_by_position(5, return_list=True)
        assert isinstance(lst, list) and len(lst) == 1

    def test_get_sites_by_position_wiggle(self, prot):
        prot.add_site(5, 'phos', value=1)
        assert prot.get_sites_by_position(3, wiggle=2, return_list=True)
        assert not prot.get_sites_by_position(1, wiggle=1, return_list=True)

    def test_get_sites_by_range(self, prot):
        prot.add_site(3, 'phos', value=1)
        prot.add_site(8, 'acet', value=1)
        assert len(prot.get_sites_by_range(1, 5, return_list=True)) == 1
        assert len(prot.get_sites_by_range(1, 20, return_list=True)) == 2

    def test_get_sites_by_type(self, prot):
        prot.add_site(3, 'phos', value=1)
        prot.add_site(8, 'phos', value=1)
        prot.add_site(9, 'acet', value=1)
        assert len(prot.get_sites_by_type('phos', return_list=True)) == 2
        assert len(prot.get_sites_by_type(['phos', 'acet'], return_list=True)) == 3

    def test_get_sites_by_type_and_range(self, prot):
        prot.add_site(3, 'phos', value=1)
        prot.add_site(18, 'phos', value=1)
        got = prot.get_sites_by_type_and_range('phos', 1, 10, return_list=True)
        assert len(got) == 1

    def test_remove_site(self, prot):
        prot.add_site(3, 'phos', value=1)
        s = prot.site(3)[0]
        prot.remove_site(s)
        assert len(prot.sites) == 0
        assert 3 not in prot.site_positions

    def test_remove_site_none_safe(self, prot):
        prot.remove_site(None, safe=False)
        with pytest.raises(ProteinException):
            prot.remove_site(None, safe=True)


# ---------------------------------------------------------------------------
# dunders
# ---------------------------------------------------------------------------
class TestProteinDunders:

    def test_len(self, prot):
        assert len(prot) == 20

    def test_repr(self, prot):
        assert 'Protein' in repr(prot)
        assert 'P1' in repr(prot)
