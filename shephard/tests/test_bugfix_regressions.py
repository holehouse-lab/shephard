"""
SHEPHARD regression suite : explicit tests for every bug fixed on the
`streamline` branch.

Each test is named after the bug it guards and fails if the underlying
bug is ever reintroduced.

Bugs covered
------------
 1.  proteome.py  : `s.position. s.site_type` typo broke copying proteins
                     that have sites into a new Proteome
 2.  track.py     : `elif(symbols, str)` (always-true tuple) broke symbol
                     track validation / string-symbol handling
 3.  track.py     : non-dict attribute_dictionary raised NameError instead
                     of TrackException
 4.  domain.py    : Domain.site() used an undefined name -> NameError
                     instead of DomainException
 5.  interface_tools.py : is_comment_line() IndexError on blank lines
 6.  sequence_tools.py  : build_mega_string ignored return_as_list
 7.  site_tools.py : build_site_density_vector ignored
                     append_leading_lagging + IndexError when
                     window_size > protein length
 8.  proteome.py  : add_proteins([]) raised IndexError
 9.  proteome.py  : all-zero values track miscopied as a symbols track
 10. domain.py    : Domain.site() raised raw KeyError for in-range
                     position with no site
 11. protein.py   : convert_to_valid always raised (sentinel counted) and
                     never wrote back when copy=False
 12. fasta.py     : shephard_fasta_to_proteome attribute parse duplicated
                     the key as the value + spurious '' attribute
 13. experimental.py : missing `random` / `site_tools` imports

Holehouse Lab - Washington University in St. Louis
"""

import pytest

import shephard
from shephard.proteome import Proteome
from shephard.track import Track
from shephard.domain import Domain
from shephard.exceptions import (TrackException, DomainException,
                                 ProteinException, ShephardException,
                                 SiteException)
from shephard.interfaces import interface_tools, si_proteins
from shephard.tools import sequence_tools, site_tools
from shephard.apis import fasta, uniprot

DATA = shephard.get_data('test_data')
SEQ = 'ACDEFGHIKLMNPQRSTVWY'


def _proteome(seq=SEQ, uid='P1'):
    return Proteome([{'sequence': seq, 'name': 'demo',
                      'unique_ID': uid, 'attributes': {}}])


# --- bug 1 -----------------------------------------------------------------
def test_bug01_copy_proteome_with_sites():
    src = _proteome()
    src.protein('P1').add_site(5, 'phos', symbol='P', value=1.0,
                               attributes={'enz': 'CK2'})
    dest = Proteome([src.protein('P1')])          # must not AttributeError
    s = dest.protein('P1').site(5)[0]
    assert s.site_type == 'phos'
    assert s.symbol == 'P'
    assert s.value == 1.0
    assert s.attribute('enz') == 'CK2'


# --- bug 2 -----------------------------------------------------------------
def test_bug02_string_symbols_track():
    p = _proteome().protein('P1')
    p.add_track('s', symbols='ABCDEFGHIJKLMNOPQRST')
    assert p.track('s').symbols == list('ABCDEFGHIJKLMNOPQRST')


# --- bug 3 -----------------------------------------------------------------
def test_bug03_track_bad_attribute_dict_raises_trackexception():
    p = _proteome().protein('P1')
    with pytest.raises(TrackException):
        Track('t', p, values=[0.0] * 20, attribute_dictionary='not-a-dict')


# --- bug 4 + 10 ------------------------------------------------------------
def test_bug04_and_10_domain_site_errors():
    p = _proteome().protein('P1')
    p.add_domain(5, 12, 'IDR')
    d = p.get_domains_by_type('IDR')[0]

    # bug 4: position outside the domain -> DomainException (not NameError)
    with pytest.raises(DomainException):
        d.site(1)

    # bug 10: in-range position with no site -> DomainException (not KeyError)
    with pytest.raises(DomainException):
        d.site(7)


# --- bug 5 -----------------------------------------------------------------
def test_bug05_is_comment_line_blank():
    assert interface_tools.is_comment_line('') is True
    assert interface_tools.is_comment_line('   \t ') is True
    assert interface_tools.is_comment_line('# c') is True
    assert interface_tools.is_comment_line('data\tline') is False


def test_bug05_blank_lines_in_file(tmp_path):
    f = tmp_path / 'p.tsv'
    f.write_text("\nA\tnameA\tACDE\n\n")
    P = Proteome()
    si_proteins.add_proteins_from_file(P, str(f))   # must not IndexError
    assert len(P) == 1


# --- bug 6 -----------------------------------------------------------------
def test_bug06_build_mega_string_return_as_list():
    P = Proteome([{'sequence': 'AA', 'name': 'n', 'unique_ID': 'P1',
                   'attributes': {}},
                  {'sequence': 'CC', 'name': 'n', 'unique_ID': 'P2',
                   'attributes': {}}])
    objs = [P.protein('P1'), P.protein('P2')]
    assert sequence_tools.build_mega_string(objs) == 'AACC'
    assert sequence_tools.build_mega_string(objs, return_as_list=True) == \
        ['AA', 'CC']


# --- bug 7 -----------------------------------------------------------------
def test_bug07_append_leading_lagging_and_window_guard():
    P = Proteome([{'sequence': 'A' * 80, 'name': 'n', 'unique_ID': 'P1',
                   'attributes': {}}])
    p = P.protein('P1')
    p.add_site(10, 'phos', value=1)

    full = site_tools.build_site_density_vector(
        p, window_size=10, append_leading_lagging=True)
    short = site_tools.build_site_density_vector(
        p, window_size=10, append_leading_lagging=False)
    assert len(full) == len(p)
    assert len(short) == len(p) - 10 + 1

    tiny = Proteome([{'sequence': 'A' * 5, 'name': 'n', 'unique_ID': 'P1',
                      'attributes': {}}])
    with pytest.raises(ShephardException):
        site_tools.build_site_density_vector(tiny.protein('P1'),
                                             window_size=30)


# --- bug 8 -----------------------------------------------------------------
def test_bug08_add_proteins_empty_list():
    P = _proteome()
    P.add_proteins([])                 # must not IndexError
    assert len(P) == 1


# --- bug 9 -----------------------------------------------------------------
def test_bug09_all_zero_values_track_copies_as_values():
    src = _proteome()
    src.protein('P1').add_track('z', values=[0.0] * 20)
    dest = Proteome([src.protein('P1')])
    t = dest.protein('P1').track('z')
    assert t.track_type == 'values'
    assert t.values == [0.0] * 20
    assert t.symbols is None


# --- bug 11 ----------------------------------------------------------------
def test_bug11_convert_to_valid_substitution_only_does_not_raise():
    # B,U,X,Z are all 1:1 substitutions -> length unchanged -> no raise,
    # and (copy=False) the underlying sequence is actually updated
    p = _proteome('ACDEBUXZ').protein('P1')
    ret = p.convert_to_valid()
    assert ret is None
    assert p.sequence == 'ACDENCGQ'
    assert p.check_sequence_is_valid() is True


def test_bug11_convert_to_valid_copy_does_not_mutate():
    p = _proteome('ACDEBUXZ').protein('P1')
    out = p.convert_to_valid(copy=True)
    assert out == 'ACDENCGQ'
    assert p.sequence == 'ACDEBUXZ'


def test_bug11_convert_to_valid_length_change_safe_guard():
    p = _proteome('ACDE*-').protein('P1')
    with pytest.raises(ProteinException):
        p.convert_to_valid()                 # would shorten -> raise
    p.convert_to_valid(safe=False)
    assert p.sequence == 'ACDE'


# --- bug 12 ----------------------------------------------------------------
def test_bug12_fasta_attribute_round_trip(tmp_path):
    P = uniprot.uniprot_fasta_to_proteome(f'{DATA}/testset_1.fasta')
    P.protein('O00401').add_attribute('gene', 'WASL')
    P.protein('O00401').add_attribute('expr', 'high=very')   # value w/ '='
    out = str(tmp_path / 'a.fasta')
    fasta.proteome_to_fasta(out, P, include_attributes_in_header=True)

    Q = fasta.shephard_fasta_to_proteome(out)
    q = Q.protein('O00401')
    assert q.attribute('gene') == 'WASL'
    assert q.attribute('expr') == 'high=very'
    assert '' not in q.attributes          # no spurious empty attribute


# --- bug 13 ----------------------------------------------------------------
def test_bug13_experimental_imports():
    from shephard.tools import experimental
    P = Proteome([{'sequence': 'A' * 100, 'name': 'n', 'unique_ID': 'P1',
                   'attributes': {}}])
    p = P.protein('P1')
    p.add_domain(20, 40, 'IDR')
    for pos in (22, 25, 30, 35):
        p.add_site(pos, 'phos', value=1)
    d = p.get_domains_by_type('IDR')[0]
    # would previously NameError on `random` / `site_tools`
    v = experimental.get_site_density_in_domain_normalized_by_protein(
        d, 'phos', sample_size=5)
    assert isinstance(v, (int, float))


# --- cross-object consistency fixes ---------------------------------------
def test_site_remove_attribute_raises_siteexception():
    # bug: Site.remove_attribute raised ProteinException
    p = _proteome().protein('P1')
    p.add_site(5, 'phos', value=1)
    s = p.site(5)[0]
    with pytest.raises(SiteException):
        s.remove_attribute('absent')


def test_track_attribute_message_says_track():
    p = _proteome().protein('P1')
    p.add_track('t', values=[0.0] * 20)
    t = p.track('t')
    with pytest.raises(TrackException) as exc:
        t.attribute('absent')
    assert 'from Track' in str(exc.value)
    assert 'from protein' not in str(exc.value)


def test_site_get_track_value_safe_false_returns_none():
    # bug: None[0] -> TypeError when track missing and safe=False
    p = _proteome().protein('P1')
    p.add_site(5, 'phos', value=1)
    s = p.site(5)[0]
    assert s.get_track_value('absent', safe=False) is None
    assert s.get_track_symbol('absent', safe=False) is None


def test_remove_track_wrong_protein_raises_proteinexception():
    # bug: error path referenced self.protein (undefined) -> AttributeError
    P = Proteome([{'sequence': SEQ, 'name': 'n', 'unique_ID': 'P1',
                   'attributes': {}},
                  {'sequence': SEQ, 'name': 'n', 'unique_ID': 'P2',
                   'attributes': {}}])
    P.protein('P1').add_track('t', values=[0.0] * 20)
    foreign = P.protein('P1').track('t')
    with pytest.raises(ProteinException):
        P.protein('P2').remove_track(foreign)
