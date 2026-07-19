"""
SHEPHARD comprehensive test suite : tools

Detailed tests covering the user-facing analysis tools:

  * attribute_tools.cast_attributes
  * domain_tools.{domain_overlap, domain_overlap_fraction,
                  domain_overlap_by_position, build_missing_domains,
                  build_domains_from_track_values}
  * site_tools.build_site_density_vector
  * sequence_tools.{build_mega_string, find_string_positions}
  * track_tools.{binerize, build_track_from_domains}
  * experimental site-density helpers (smoke)

Holehouse Lab - Washington University in St. Louis
"""

import pytest

from shephard.proteome import Proteome
from shephard.exceptions import DomainException, ShephardException
from shephard.tools import (attribute_tools, domain_tools, site_tools,
                            sequence_tools, track_tools, experimental)
from shephard.interfaces import si_domains

SEQ = 'ACDEFGHIKLMNPQRSTVWY'


@pytest.fixture
def prot():
    P = Proteome([{'sequence': SEQ, 'name': 'demo',
                   'unique_ID': 'P1', 'attributes': {}}])
    return P.protein('P1')


# ===========================================================================
# attribute_tools.cast_attributes
# ===========================================================================
class TestCastAttributes:

    def test_default_casts_all(self, prot):
        prot.add_attribute('a', '1')
        prot.add_attribute('b', '2.5')
        attribute_tools.cast_attributes(prot, cast_type=float)
        assert prot.attribute('a') == 1.0
        assert prot.attribute('b') == 2.5

    def test_include_only(self, prot):
        prot.add_attribute('a', '1')
        prot.add_attribute('b', '2')
        attribute_tools.cast_attributes(prot, include=['a'], cast_type=int)
        assert prot.attribute('a') == 1
        assert prot.attribute('b') == '2'      # untouched

    def test_exclude(self, prot):
        prot.add_attribute('a', '1')
        prot.add_attribute('b', '2')
        attribute_tools.cast_attributes(prot, exclude=['b'], cast_type=int)
        assert prot.attribute('a') == 1
        assert prot.attribute('b') == '2'

    def test_include_and_exclude_raises(self, prot):
        prot.add_attribute('a', '1')
        with pytest.raises(ValueError):
            attribute_tools.cast_attributes(prot, include=['a'],
                                            exclude=['a'])

    def test_bad_cast_type_raises(self, prot):
        prot.add_attribute('a', '1')
        with pytest.raises(ValueError):
            attribute_tools.cast_attributes(prot, cast_type='float')

    def test_skip_failed(self, prot):
        prot.add_attribute('a', 'not-a-number')
        prot.add_attribute('b', '3')
        attribute_tools.cast_attributes(prot, cast_type=float,
                                        skip_failed=True)
        assert prot.attribute('a') == 'not-a-number'   # left as-is
        assert prot.attribute('b') == 3.0

    def test_failed_without_skip_raises(self, prot):
        prot.add_attribute('a', 'not-a-number')
        with pytest.raises(ValueError):
            attribute_tools.cast_attributes(prot, cast_type=float)


# ===========================================================================
# domain_tools
# ===========================================================================
class TestDomainTools:

    def test_overlap_by_position(self):
        assert domain_tools.domain_overlap_by_position(1, 10, 5, 15) is True
        assert domain_tools.domain_overlap_by_position(1, 5, 10, 15) is False
        assert domain_tools.domain_overlap_by_position(1, 10, 10, 20) is True

    def test_domain_overlap_same_protein(self, prot):
        prot.add_domain(1, 10, 'A')
        prot.add_domain(5, 15, 'B')
        a = prot.get_domains_by_type('A')[0]
        b = prot.get_domains_by_type('B')[0]
        assert domain_tools.domain_overlap(a, b) is True

    def test_domain_overlap_cross_protein_raises(self):
        P = Proteome([{'sequence': SEQ, 'name': 'n', 'unique_ID': 'P1',
                       'attributes': {}},
                      {'sequence': SEQ, 'name': 'n', 'unique_ID': 'P2',
                       'attributes': {}}])
        P.protein('P1').add_domain(1, 10, 'A')
        P.protein('P2').add_domain(1, 10, 'B')
        a = P.protein('P1').get_domains_by_type('A')[0]
        b = P.protein('P2').get_domains_by_type('B')[0]
        with pytest.raises(DomainException):
            domain_tools.domain_overlap(a, b)
        # check_origin=False suppresses the cross-protein guard
        assert domain_tools.domain_overlap(a, b, check_origin=False) is True

    def test_domain_overlap_fraction(self, prot):
        prot.add_domain(1, 10, 'A')
        prot.add_domain(5, 14, 'B')
        a = prot.get_domains_by_type('A')[0]
        b = prot.get_domains_by_type('B')[0]
        frac = domain_tools.domain_overlap_fraction(a, b)
        assert 0.0 <= frac <= 1.0

    def test_build_missing_domains(self, prot):
        prot.add_domain(5, 10, 'known')
        missing = domain_tools.build_missing_domains(prot,
                                                     new_domain_type='gap')
        assert isinstance(missing, list)
        assert all(m['domain_type'] == 'gap' for m in missing)
        # add them back and confirm full coverage of the non-domain residues
        prot.add_domains(missing)
        covered = set()
        for d in prot.domains:
            covered.update(range(d.start, d.end + 1))
        assert covered == set(range(1, len(prot) + 1))

    def test_build_domains_from_track_values(self):
        P = Proteome([{'sequence': 'A' * 60, 'name': 'n',
                       'unique_ID': 'P1', 'attributes': {}}])
        vals = [1.0] * 30 + [0.0] * 30
        P.protein('P1').add_track('t', values=vals)

        def binf(v):
            return [1 if x > 0.5 else 0 for x in v]

        dd = domain_tools.build_domains_from_track_values(
            P, 't', binf, 'hit', minimum_region_size=10, gap_closure=1,
            verbose=False)
        assert isinstance(dd, dict)
        si_domains.add_domains_from_dictionary(P, dd)
        doms = P.protein('P1').get_domains_by_type('hit')
        assert len(doms) == 1
        assert doms[0].start == 1


# ===========================================================================
# site_tools.build_site_density_vector
# ===========================================================================
class TestSiteDensityVector:

    def test_basic_density_length_matches_protein(self):
        P = Proteome([{'sequence': 'A' * 100, 'name': 'n',
                       'unique_ID': 'P1', 'attributes': {}}])
        p = P.protein('P1')
        for pos in (10, 11, 12, 50):
            p.add_site(pos, 'phos', value=1)
        vec = site_tools.build_site_density_vector(p, window_size=10)
        assert len(vec) == len(p)
        assert max(vec) > 0

    def test_site_type_filter(self):
        P = Proteome([{'sequence': 'A' * 60, 'name': 'n',
                       'unique_ID': 'P1', 'attributes': {}}])
        p = P.protein('P1')
        p.add_site(10, 'phos', value=1)
        p.add_site(20, 'acet', value=1)
        v_all = site_tools.build_site_density_vector(p, window_size=10)
        v_phos = site_tools.build_site_density_vector(
            p, site_types='phos', window_size=10)
        assert sum(v_all) > sum(v_phos)

    def test_append_leading_lagging_false_returns_shorter(self):
        # regression: this flag was previously ignored
        P = Proteome([{'sequence': 'A' * 100, 'name': 'n',
                       'unique_ID': 'P1', 'attributes': {}}])
        p = P.protein('P1')
        p.add_site(10, 'phos', value=1)
        full = site_tools.build_site_density_vector(
            p, window_size=10, append_leading_lagging=True)
        short = site_tools.build_site_density_vector(
            p, window_size=10, append_leading_lagging=False)
        assert len(full) == len(p)
        assert len(short) == len(p) - 10 + 1
        assert len(short) < len(full)

    def test_window_larger_than_protein_raises(self):
        # regression: previously an IndexError
        P = Proteome([{'sequence': 'A' * 10, 'name': 'n',
                       'unique_ID': 'P1', 'attributes': {}}])
        with pytest.raises(ShephardException):
            site_tools.build_site_density_vector(P.protein('P1'),
                                                 window_size=30)


# ===========================================================================
# sequence_tools
# ===========================================================================
class TestSequenceTools:

    def test_build_mega_string(self):
        P = Proteome([{'sequence': 'AAA', 'name': 'n', 'unique_ID': 'P1',
                       'attributes': {}},
                      {'sequence': 'CCC', 'name': 'n', 'unique_ID': 'P2',
                       'attributes': {}}])
        objs = [P.protein('P1'), P.protein('P2')]
        assert sequence_tools.build_mega_string(objs) == 'AAACCC'

    def test_build_mega_string_return_as_list(self):
        # regression: return_as_list was previously ignored
        P = Proteome([{'sequence': 'AAA', 'name': 'n', 'unique_ID': 'P1',
                       'attributes': {}},
                      {'sequence': 'CCC', 'name': 'n', 'unique_ID': 'P2',
                       'attributes': {}}])
        objs = [P.protein('P1'), P.protein('P2')]
        out = sequence_tools.build_mega_string(objs, return_as_list=True)
        assert out == ['AAA', 'CCC']

    def test_find_string_positions_protein_indexed(self):
        # protein indexing -> 1-based
        assert sequence_tools.find_string_positions('AC', 'ACXXAC') == [1, 5]

    def test_find_string_positions_zero_indexed(self):
        assert sequence_tools.find_string_positions(
            'AC', 'ACXXAC', protein_indexing=False) == [0, 4]

    def test_find_string_positions_overlapping(self):
        assert sequence_tools.find_string_positions(
            'AA', 'AAAA', protein_indexing=False) == [0, 1, 2]

    def test_find_string_positions_regex(self):
        # 'L.P' matches an L, any residue, then P
        assert sequence_tools.find_string_positions(
            'L.P', 'XLAPX', protein_indexing=False) == [1]

    def test_find_string_positions_no_match(self):
        assert sequence_tools.find_string_positions('ZZ', 'ACDE') == []


# ===========================================================================
# track_tools
# ===========================================================================
class TestTrackTools:

    def test_binerize_above(self):
        assert track_tools.binerize([0, 1, 2, 3], threshold=1) == [0, 0, 1, 1]

    def test_binerize_below(self):
        assert track_tools.binerize([0, 1, 2, 3], threshold=2,
                                    mode='below') == [1, 1, 0, 0]

    def test_binerize_bad_mode_raises(self):
        with pytest.raises(Exception):
            track_tools.binerize([1, 2], threshold=1, mode='sideways')

    def test_build_track_from_domains_all(self):
        P = Proteome([{'sequence': 'A' * 20, 'name': 'n',
                       'unique_ID': 'P1', 'attributes': {}}])
        P.protein('P1').add_domain(1, 5, 'D')
        P.protein('P1').add_domain(10, 12, 'E')
        tracks = track_tools.build_track_from_domains(P)
        assert 'P1' in tracks
        symbols = tracks['P1']
        assert len(symbols) == 20
        assert symbols[0] == '1' and symbols[5] == '0'
        assert symbols[10] == '1'

    def test_build_track_from_domains_by_type(self):
        P = Proteome([{'sequence': 'A' * 20, 'name': 'n',
                       'unique_ID': 'P1', 'attributes': {}}])
        P.protein('P1').add_domain(1, 5, 'D')
        P.protein('P1').add_domain(10, 12, 'E')
        tracks = track_tools.build_track_from_domains(P, domain_type='D')
        sym = tracks['P1']
        assert sym[0] == '1'
        assert sym[10] == '0'      # 'E' domain excluded


# ===========================================================================
# experimental (smoke)
# ===========================================================================
class TestExperimental:

    def _proteome_with_domain_and_sites(self):
        P = Proteome([{'sequence': 'A' * 100, 'name': 'n',
                       'unique_ID': 'P1', 'attributes': {}}])
        p = P.protein('P1')
        p.add_domain(20, 40, 'IDR')
        for pos in (22, 25, 30, 35):
            p.add_site(pos, 'phos', value=1)
        return p.get_domains_by_type('IDR')[0]

    def test_site_density_enrichment_runs(self):
        d = self._proteome_with_domain_and_sites()
        val = experimental.get_site_density_in_domain_normalized_by_protein(
            d, 'phos', sample_size=5)
        assert isinstance(val, (int, float))

    def test_site_density_percentile_runs(self):
        d = self._proteome_with_domain_and_sites()
        val = experimental.get_site_density_percentile_normalized_by_protein(
            d, 'phos')
        assert val is None or (0.0 <= val <= 1.0)
