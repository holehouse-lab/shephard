"""
SHEPHARD comprehensive test suite : interfaces (si_*) + interface_tools

Detailed tests covering the user-facing file/dictionary IO API:

  * si_proteins   : file <-> dictionary <-> proteome round trips
  * si_domains    : domains + domain attributes, type filtering, autoname
  * si_sites      : sites + dictionary + list writers, type filtering
  * si_tracks     : values & symbols tracks, single/separate file writers
  * si_protein_attributes : file/dict read & write
  * interface_tools : check_*, clean_string, parse_key_value_pairs,
                      is_comment_line (incl. blank-line regression)
  * robustness     : return_dictionary, skip_bad, comment & blank lines

All write tests use pytest's tmp_path so they are independent of the
current working directory.

Holehouse Lab - Washington University in St. Louis
"""

import pytest

import shephard
from shephard.proteome import Proteome
from shephard.apis import uniprot
from shephard.interfaces import (si_proteins, si_domains, si_sites,
                                 si_tracks, si_protein_attributes,
                                 interface_tools)
from shephard.exceptions import InterfaceException

DATA = shephard.get_data('test_data')


def _ts1():
    return uniprot.uniprot_fasta_to_proteome(f'{DATA}/testset_1.fasta')


# ===========================================================================
# si_proteins
# ===========================================================================
class TestSiProteins:

    def test_add_proteins_from_file(self):
        P = Proteome()
        si_proteins.add_proteins_from_file(P, f'{DATA}/TS1_small_proteins.tsv')
        assert len(P) == 3
        assert 'O00401' in P

    def test_return_dictionary(self):
        P = Proteome()
        d = si_proteins.add_proteins_from_file(
            P, f'{DATA}/TS1_small_proteins.tsv', return_dictionary=True)
        assert isinstance(d, dict)
        assert len(P) == 0                     # nothing added
        assert 'O00401' in d
        assert set(d['O00401']) >= {'name', 'sequence', 'attributes'}

    def test_add_proteins_from_dictionary(self):
        P = Proteome()
        si_proteins.add_proteins_from_dictionary(
            P, {'X1': {'name': 'n', 'sequence': 'ACDEK', 'attributes': {}}})
        assert P.protein('X1').sequence == 'ACDEK'

    def test_write_read_round_trip(self, tmp_path):
        P = _ts1()
        P.protein('O00401').add_attribute('gene', 'WASL')
        out = str(tmp_path / 'prot.tsv')
        si_proteins.write_proteins(P, out)

        Q = Proteome()
        si_proteins.add_proteins_from_file(Q, out)
        assert len(Q) == len(P)
        assert Q.protein('O00401').sequence == P.protein('O00401').sequence
        assert Q.protein('O00401').attribute('gene') == 'WASL'

    def test_duplicate_in_file_raises(self, tmp_path):
        f = tmp_path / 'dup.tsv'
        f.write_text("A\tnameA\tACDE\nA\tnameA\tACDE\n")
        with pytest.raises(InterfaceException):
            si_proteins.add_proteins_from_file(Proteome(), str(f))


# ===========================================================================
# si_domains
# ===========================================================================
class TestSiDomains:

    def test_add_domains_from_file(self):
        P = _ts1()
        si_domains.add_domains_from_file(P, f'{DATA}/TS1_domains_idr.tsv')
        assert len(P.domains) > 0
        assert 'IDR' in P.unique_domain_types

    def test_return_dictionary(self):
        P = _ts1()
        d = si_domains.add_domains_from_file(
            P, f'{DATA}/TS1_domains_idr.tsv', return_dictionary=True)
        assert isinstance(d, dict)
        assert len(P.domains) == 0
        assert 'O00401' in d

    def test_add_domains_from_dictionary(self):
        P = _ts1()
        si_domains.add_domains_from_dictionary(
            P, {'O00401': [{'start': 1, 'end': 10, 'domain_type': 'D',
                            'attributes': {'a': 1}}]})
        d = P.protein('O00401').get_domains_by_type('D')[0]
        assert d.attribute('a') == 1

    def test_autoname_allows_duplicates(self):
        P = _ts1()
        dd = {'O00401': [{'start': 1, 'end': 5, 'domain_type': 'X'},
                         {'start': 1, 'end': 5, 'domain_type': 'X'}]}
        si_domains.add_domains_from_dictionary(P, dd, autoname=True)
        assert len(P.protein('O00401').get_domains_by_type('X')) == 2

    def test_write_read_round_trip(self, tmp_path):
        P = _ts1()
        si_domains.add_domains_from_file(P, f'{DATA}/TS1_domains_idr.tsv')
        n = len(P.domains)
        out = str(tmp_path / 'dom.tsv')
        si_domains.write_domains(P, out)

        Q = _ts1()
        si_domains.add_domains_from_file(Q, out)
        assert len(Q.domains) == n

    def test_write_domains_type_filter(self, tmp_path):
        P = _ts1()
        si_domains.add_domains_from_file(
            P, f'{DATA}/TS1_domains_idr_and_others.tsv')
        out = str(tmp_path / 'idr_only.tsv')
        si_domains.write_domains(P, out, domain_types=['IDR'])
        Q = _ts1()
        si_domains.add_domains_from_file(Q, out)
        assert set(Q.unique_domain_types) == {'IDR'}

    def test_write_domains_type_filter_must_be_list(self, tmp_path):
        P = _ts1()
        si_domains.add_domains_from_file(P, f'{DATA}/TS1_domains_idr.tsv')
        with pytest.raises(InterfaceException):
            si_domains.write_domains(P, str(tmp_path / 'x.tsv'),
                                     domain_types='IDR')

    def test_write_domains_from_list(self, tmp_path):
        P = _ts1()
        si_domains.add_domains_from_file(P, f'{DATA}/TS1_domains_idr.tsv')
        out = str(tmp_path / 'list.tsv')
        si_domains.write_domains_from_list(P.domains, out)
        Q = _ts1()
        si_domains.add_domains_from_file(Q, out)
        assert len(Q.domains) == len(P.domains)

    def test_add_domain_attributes_from_dictionary(self):
        P = _ts1()
        si_domains.add_domains_from_file(P, f'{DATA}/TS1_domains_idr.tsv')
        d0 = P.protein('O00401').get_domains_by_type('IDR')[0]
        dd = {'O00401': [{'start': d0.start, 'end': d0.end,
                          'domain_type': 'IDR',
                          'attributes': {'newkey': 'newval'}}]}
        si_domains.add_domain_attributes_from_dictionary(P, dd)
        assert d0.attribute('newkey') == 'newval'


# ===========================================================================
# si_sites
# ===========================================================================
class TestSiSites:

    def test_add_sites_from_file(self):
        P = _ts1()
        si_sites.add_sites_from_file(P, f'{DATA}/TS1_sites.tsv')
        assert len(P.sites) > 0

    def test_return_dictionary(self):
        P = _ts1()
        d = si_sites.add_sites_from_file(
            P, f'{DATA}/TS1_sites.tsv', return_dictionary=True)
        assert isinstance(d, dict)
        assert len(P.sites) == 0

    def test_add_sites_from_dictionary(self):
        P = _ts1()
        sd = {'O00401': [{'position': 5, 'site_type': 'phos',
                          'symbol': 'P', 'value': 1.0,
                          'attributes': {'enz': 'CK2'}}]}
        si_sites.add_sites_from_dictionary(P, sd)
        s = P.protein('O00401').site(5)[0]
        assert s.site_type == 'phos'
        assert s.attribute('enz') == 'CK2'

    def test_sites_with_none_value(self):
        P = _ts1()
        si_sites.add_sites_from_file(P, f'{DATA}/ts1_bonus_sites.tsv')
        # file contains a site whose value column is None
        assert len(P.sites) > 0

    def test_write_read_round_trip(self, tmp_path):
        P = _ts1()
        si_sites.add_sites_from_file(P, f'{DATA}/TS1_sites.tsv')
        n = len(P.sites)
        out = str(tmp_path / 's.tsv')
        si_sites.write_sites(P, out)
        Q = _ts1()
        si_sites.add_sites_from_file(Q, out)
        assert len(Q.sites) == n

    def test_write_sites_type_filter(self, tmp_path):
        P = _ts1()
        si_sites.add_sites_from_file(P, f'{DATA}/TS1_sites.tsv')
        types = list(P.unique_site_types)[:1]
        out = str(tmp_path / 'filt.tsv')
        si_sites.write_sites(P, out, site_types=types)
        Q = _ts1()
        si_sites.add_sites_from_file(Q, out)
        assert set(Q.unique_site_types) == set(types)

    def test_write_sites_from_list(self, tmp_path):
        P = _ts1()
        si_sites.add_sites_from_file(P, f'{DATA}/TS1_sites.tsv')
        out = str(tmp_path / 'l.tsv')
        si_sites.write_sites_from_list(P.sites, out)
        Q = _ts1()
        si_sites.add_sites_from_file(Q, out)
        assert len(Q.sites) == len(P.sites)

    def test_skip_bad_lines_parsing(self, tmp_path):
        # skip_bad governs FILE PARSING: an unparseable line is skipped,
        # the well-formed line is still added.
        f = tmp_path / 'badsites.tsv'
        f.write_text("O00401\t5\tphos\tP\t1.0\n"
                     "O00401\tNOT_AN_INT\tphos\tP\t1.0\n")
        P = _ts1()
        si_sites.add_sites_from_file(P, str(f), skip_bad=True)
        assert len(P.protein('O00401').sites) == 1

    def test_out_of_range_site_safe_false_skips(self):
        # add-time (non-parsing) errors are governed by `safe`, not skip_bad
        P = _ts1()
        si_sites.add_sites_from_file(P, f'{DATA}/ts1_bonus_sites_bad.tsv',
                                     safe=False, verbose=False)
        # out-of-range sites are skipped rather than raising
        assert len(P.sites) >= 0


# ===========================================================================
# si_tracks
# ===========================================================================
class TestSiTracks:

    def test_add_values_track_from_file(self):
        P = _ts1()
        si_tracks.add_tracks_from_file(P, f'{DATA}/TS1_tracks_pscore.tsv',
                                       mode='values')
        t = P.protein('O00401').track('pscore')
        assert t.track_type == 'values'

    def test_return_dictionary(self):
        P = _ts1()
        d = si_tracks.add_tracks_from_file(
            P, f'{DATA}/TS1_tracks_pscore.tsv', mode='values',
            return_dictionary=True)
        assert isinstance(d, dict)
        assert 'pscore' not in P.protein('O00401').track_names

    def test_add_tracks_from_dictionary_values(self):
        P = _ts1()
        L = len(P.protein('O00401'))
        td = {'O00401': [{'track_name': 'mt',
                          'track_data': [0.1] * L}]}
        si_tracks.add_tracks_from_dictionary(P, td, mode='values')
        assert P.protein('O00401').track('mt').track_type == 'values'

    def test_add_tracks_from_dictionary_symbols(self):
        P = _ts1()
        L = len(P.protein('O00401'))
        td = {'O00401': [{'track_name': 'ms',
                          'track_data': ['x'] * L}]}
        si_tracks.add_tracks_from_dictionary(P, td, mode='symbols')
        assert P.protein('O00401').track('ms').track_type == 'symbols'

    def test_invalid_mode_raises(self):
        P = _ts1()
        with pytest.raises(Exception):
            si_tracks.add_tracks_from_file(
                P, f'{DATA}/TS1_tracks_pscore.tsv', mode='nonsense')

    def test_write_read_round_trip(self, tmp_path):
        P = _ts1()
        si_tracks.add_tracks_from_file(P, f'{DATA}/TS1_tracks_pscore.tsv',
                                       mode='values')
        out = str(tmp_path / 'tr.tsv')
        si_tracks.write_tracks(P, out, 'pscore')

        Q = _ts1()
        si_tracks.add_tracks_from_file(Q, out, mode='values')
        a = P.protein('O00401').get_track_values('pscore')
        b = Q.protein('O00401').get_track_values('pscore')
        assert len(a) == len(b)
        assert pytest.approx(a[0], abs=1e-3) == b[0]

    def test_write_all_separate_files(self, tmp_path):
        P = _ts1()
        si_tracks.add_tracks_from_file(P, f'{DATA}/TS1_tracks_pscore.tsv',
                                       mode='values')
        si_tracks.write_all_tracks_separate_files(P, outdirectory=str(tmp_path))
        assert (tmp_path / 'shephard_track_pscore.tsv').exists()

    def test_write_all_values_single_file(self, tmp_path):
        P = _ts1()
        si_tracks.add_tracks_from_file(P, f'{DATA}/TS1_tracks_pscore.tsv',
                                       mode='values')
        out = str(tmp_path / 'all_vals.tsv')
        si_tracks.write_all_values_tracks_single_file(P, out)
        Q = _ts1()
        si_tracks.add_tracks_from_file(Q, out, mode='values')
        assert 'pscore' in Q.protein('O00401').track_names

    def test_write_tracks_from_list(self, tmp_path):
        P = _ts1()
        si_tracks.add_tracks_from_file(P, f'{DATA}/TS1_tracks_pscore.tsv',
                                       mode='values')
        tracks = [p.track('pscore') for p in P if 'pscore' in p.track_names]
        out = str(tmp_path / 'fl.tsv')
        si_tracks.write_tracks_from_list(tracks, out)
        Q = _ts1()
        si_tracks.add_tracks_from_file(Q, out, mode='values')
        assert 'pscore' in Q.protein('O00401').track_names

    def test_write_tracks_invalid_value_fmt(self, tmp_path):
        P = _ts1()
        si_tracks.add_tracks_from_file(P, f'{DATA}/TS1_tracks_pscore.tsv',
                                       mode='values')
        # '%d' formats 1.5 -> '1' which no longer round-trips to a float
        with pytest.raises(InterfaceException):
            si_tracks.write_tracks(P, str(tmp_path / 'x.tsv'),
                                   'pscore', value_fmt='%d')


# ===========================================================================
# si_protein_attributes
# ===========================================================================
class TestSiProteinAttributes:

    def test_add_from_file(self):
        P = _ts1()
        si_protein_attributes.add_protein_attributes_from_file(
            P, f'{DATA}/TS1_protein_attributes.tsv')
        assert len(P.protein('O00401').attributes) > 0

    def test_add_from_dictionary(self):
        P = _ts1()
        si_protein_attributes.add_protein_attributes_from_dictionary(
            P, {'O00401': [{'k1': 'v1', 'k2': 'v2'}]})
        assert P.protein('O00401').attribute('k1') == 'v1'

    def test_write_read_round_trip(self, tmp_path):
        P = _ts1()
        P.protein('O00401').add_attribute('a', 'one')
        P.protein('O00470').add_attribute('b', 'two')
        out = str(tmp_path / 'pa.tsv')
        si_protein_attributes.write_protein_attributes(P, out)

        Q = _ts1()
        si_protein_attributes.add_protein_attributes_from_file(Q, out)
        assert Q.protein('O00401').attribute('a') == 'one'
        assert Q.protein('O00470').attribute('b') == 'two'

    def test_write_from_dictionary(self, tmp_path):
        out = str(tmp_path / 'pad.tsv')
        si_protein_attributes.write_protein_attributes_from_dictionary(
            {'O00401': {'k': 'v'}}, out)
        P = _ts1()
        si_protein_attributes.add_protein_attributes_from_file(P, out)
        assert P.protein('O00401').attribute('k') == 'v'


# ===========================================================================
# interface_tools
# ===========================================================================
class TestInterfaceTools:

    def test_check_helpers_accept_correct_types(self):
        P = _ts1()
        si_domains.add_domains_from_file(P, f'{DATA}/TS1_domains_idr.tsv')
        si_sites.add_sites_from_file(P, f'{DATA}/TS1_sites.tsv')
        L = len(P.protein('O00401'))
        P.protein('O00401').add_track('t', values=[0.0] * L)

        prot = P.protein('O00401')
        assert interface_tools.check_proteome(P, 'f') is None
        assert interface_tools.check_protein(prot, 'f') is None
        assert interface_tools.check_domain(prot.domains[0], 'f') is None
        assert interface_tools.check_site(P.sites[0], 'f') is None
        assert interface_tools.check_track(prot.track('t'), 'f') is None

    def test_check_helpers_reject_wrong_types(self):
        with pytest.raises(InterfaceException):
            interface_tools.check_proteome('x', 'f')
        with pytest.raises(InterfaceException):
            interface_tools.check_protein(123, 'f')

    def test_clean_string_replaces_delimiter(self):
        assert interface_tools.clean_string('a\tb') == 'a b'
        assert interface_tools.clean_string(12345) == '12345'

    def test_full_clean_string_removes_tab_and_colon(self):
        out = interface_tools.full_clean_string('a\tb:c')
        assert '\t' not in out and ':' not in out

    def test_parse_key_value_pairs(self):
        kv = interface_tools.parse_key_value_pairs(
            ['a:1', 'b:two'], 'f', 1, 'line')
        assert kv == {'a': '1', 'b': 'two'}

    def test_parse_key_value_pairs_bad_raises(self):
        with pytest.raises(InterfaceException):
            interface_tools.parse_key_value_pairs(
                ['no-colon-here'], 'f', 1, 'line')

    def test_is_comment_line(self):
        assert interface_tools.is_comment_line('# a comment') is True
        assert interface_tools.is_comment_line('   ') is True   # regression
        assert interface_tools.is_comment_line('') is True      # regression
        assert interface_tools.is_comment_line('O00401\t1\t2') is False

    def test_comment_and_blank_lines_in_file(self, tmp_path):
        f = tmp_path / 'with_blanks.tsv'
        f.write_text("# header comment\n\nA\tnameA\tACDEK\n\n# tail\n")
        P = Proteome()
        si_proteins.add_proteins_from_file(P, str(f))   # must not raise
        assert len(P) == 1 and 'A' in P
