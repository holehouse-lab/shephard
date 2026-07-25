"""
SHEPHARD regression suite : explicit tests for every bug fixed for the 0.2.4
release.

Each test is named after the bug it guards and fails if the underlying bug is
ever reintroduced. Bugs fixed in earlier releases are covered in
test_bugfix_regressions.py (streamline branch) and
test_bugfix_regressions_0_2_3.py.

Bugs covered
------------
 1. proteome.py : Proteome[-1] raised a KeyError and negative slices raised a
                  ValueError, so a Proteome could not be indexed from the end
 2. protein.py  : over-writing a Track (safe=False) left the replaced Track
                  counted in Proteome.unique_track_names forever
 3. protein.py  : over-writing a Domain (safe=False) left the replaced Domain
                  counted in Proteome.unique_domain_types forever
 4. protein.py  : Protein.build_domain() accepted and documented safe and
                  autoname, then dropped both
 5. proteome.py : Track attributes were lost when Proteins were copied into a
                  new Proteome (and Protein.add_track() had no way to set them)
 6. protein.py  : Protein.sites was documented as sorted N to C terminal but
                  returned Sites in the order they happened to be added
 7. proteome.py : three internal book-keeping exceptions were missing their
                  f-string prefix, so they printed literal {} placeholders
 8. site.py     : Site.get_track_values() on a symbols track (and the reverse)
                  raised a raw TypeError rather than a SiteException, as did
                  Track.values_region()/symbols_region() themselves
 9. si_sites.py : a Site with no symbol was written as 'None' and read back in
                  as the literal string 'None'

Holehouse Lab - Washington University in St. Louis
"""

import os

import pytest

from shephard.proteome import Proteome
from shephard.exceptions import (SiteException, ProteinException,
                                 ProteomeException, TrackException,
                                 DomainException)
from shephard.interfaces import si_sites

SEQ = 'ACDEFGHIKLMNPQRSTVWY'


def _proteome(n=1):
    return Proteome([{'sequence': SEQ, 'name': f'demo{i}',
                      'unique_ID': f'P{i}', 'attributes': {}} for i in range(n)])


# --- bug 1 -----------------------------------------------------------------
def test_bug01_negative_indexing():
    P = _proteome(3)

    assert P[-1].unique_ID == 'P2'
    assert P[-3].unique_ID == 'P0'
    assert [p.unique_ID for p in P[-2:]] == ['P1', 'P2']
    assert [p.unique_ID for p in P[:-1]] == ['P0', 'P1']

    # out of range in either direction is an IndexError, as for a list
    with pytest.raises(IndexError):
        P[3]
    with pytest.raises(IndexError):
        P[-4]


# --- bug 2 -----------------------------------------------------------------
def test_bug02_track_overwrite_does_not_leak_track_name():
    P = _proteome()
    protein = P.protein('P0')

    protein.add_track('t1', values=[1]*len(SEQ))
    protein.add_track('t1', values=[2]*len(SEQ), safe=False)

    # one Track exists, so the Proteome should count exactly one
    assert protein.track('t1').values == [2.0]*len(SEQ)
    assert P.unique_track_names == ['t1']

    # ... and removing it should remove the name from the Proteome entirely
    protein.remove_track(protein.track('t1'))
    assert P.unique_track_names == []
    assert P.track_names_to_track_type == {}


def test_bug02_track_overwrite_cannot_change_track_type():
    # a track name is consistently a values or a symbols track across a whole
    # Proteome, so over-writing with the other type is an error - and must
    # leave the pre-existing Track untouched
    P = _proteome()
    protein = P.protein('P0')

    protein.add_track('t1', values=[1]*len(SEQ))

    with pytest.raises(ProteomeException):
        protein.add_track('t1', symbols='A'*len(SEQ), safe=False)

    assert protein.track('t1').track_type == 'values'
    assert P.track_names_to_track_type == {'t1': 'values'}


def test_bug02_failed_track_overwrite_leaves_original_intact():
    P = _proteome()
    protein = P.protein('P0')
    protein.add_track('t1', values=[1]*len(SEQ))

    # too short, so building the replacement fails - the original must survive
    with pytest.raises(Exception):
        protein.add_track('t1', values=[9, 9], safe=False)

    assert protein.track('t1').values == [1.0]*len(SEQ)
    assert P.unique_track_names == ['t1']


# --- bug 3 -----------------------------------------------------------------
def test_bug03_domain_overwrite_does_not_leak_domain_type():
    P = _proteome()
    protein = P.protein('P0')

    protein.add_domain(1, 5, 'dtype')
    protein.add_domain(1, 5, 'dtype', safe=False)

    assert len(protein.domains) == 1
    assert P.unique_domain_types == ['dtype']

    protein.remove_domain(protein.domain('dtype_1_5'))
    assert P.unique_domain_types == []


# --- bug 4 -----------------------------------------------------------------
def test_bug04_build_domain_honours_autoname_and_safe():

    def builder(_):
        return [{'start': 2, 'end': 6, 'domain_type': 'bd'}]

    P = _proteome()
    protein = P.protein('P0')

    protein.build_domain(None, builder)

    # with autoname the second (identical) domain gets its own name
    protein.build_domain(None, builder, autoname=True)
    assert sorted(protein.domain_names) == ['bd_2_6', 'bd_2_6_1']

    # and with safe=False a duplicate is quietly over-written rather than raising
    protein.build_domain(None, builder, safe=False)
    assert sorted(protein.domain_names) == ['bd_2_6', 'bd_2_6_1']

    # while the default (safe=True, autoname=False) still raises
    with pytest.raises(ProteinException):
        protein.build_domain(None, builder)


# --- bug 5 -----------------------------------------------------------------
def test_bug05_track_attributes_survive_a_proteome_copy():
    P = _proteome()
    protein = P.protein('P0')

    protein.add_track('t1', values=[1]*len(SEQ), attributes={'origin': 'test'})
    protein.add_track('t2', symbols='A'*len(SEQ))
    protein.track('t2').add_attribute('nested', {'a': [1, 2]})

    P2 = Proteome([protein])

    assert P2.protein('P0').track('t1').attribute('origin') == 'test'
    assert P2.protein('P0').track('t2').attribute('nested') == {'a': [1, 2]}

    # attributes must be copied by value, not shared with the original
    P2.protein('P0').track('t2').attribute('nested')['a'].append(3)
    assert protein.track('t2').attribute('nested')['a'] == [1, 2]


# --- bug 6 -----------------------------------------------------------------
def test_bug06_protein_sites_are_sorted_n_to_c():
    P = _proteome()
    protein = P.protein('P0')

    for position in [12, 3, 20, 1, 7]:
        protein.add_site(position, 'phos')

    assert [s.position for s in protein.sites] == [1, 3, 7, 12, 20]

    # and this must hold at the Proteome level too
    assert [s.position for s in P.sites] == [1, 3, 7, 12, 20]


# --- bug 7 -----------------------------------------------------------------
def test_bug07_bookkeeping_errors_are_formatted():
    P = _proteome()
    protein = P.protein('P0')
    protein.add_track('t1', values=[1]*len(SEQ))

    # force the internal book-keeping into the "this should never happen" path
    P._unique_track_names = {}

    with pytest.raises(Exception) as excinfo:
        protein.remove_track(protein.track('t1'))

    # the message must contain the actual track name, not a literal placeholder
    assert 't1' in str(excinfo.value)
    assert '{track_name}' not in str(excinfo.value)


# --- bug 8 -----------------------------------------------------------------
def test_bug08_site_track_type_mismatch_raises_site_exception():
    P = _proteome()
    protein = P.protein('P0')
    protein.add_site(5, 'phos')
    protein.add_track('vals', values=[1]*len(SEQ))
    protein.add_track('syms', symbols='A'*len(SEQ))

    site = protein.sites[0]

    with pytest.raises(SiteException):
        site.get_track_values('syms')

    with pytest.raises(SiteException):
        site.get_track_symbols('vals')

    # the correct pairing still works
    assert site.get_track_value('vals') == 1.0
    assert site.get_track_symbol('syms') == 'A'

    # ... as does asking for a track that does not exist with safe=False
    assert site.get_track_values('missing', safe=False) is None
    assert site.get_track_symbols('missing', safe=False) is None


def test_bug08_track_region_type_mismatch_raises_track_exception():
    # the same mistake made directly against a Track is also an informative
    # exception rather than a TypeError
    P = _proteome()
    protein = P.protein('P0')
    protein.add_track('vals', values=[1]*len(SEQ))
    protein.add_track('syms', symbols='A'*len(SEQ))

    with pytest.raises(TrackException):
        protein.track('syms').values_region(1, 5)

    with pytest.raises(TrackException):
        protein.track('vals').symbols_region(1, 5)

    # ... while a Domain still reports this as a DomainException
    protein.add_domain(1, 5, 'd')
    domain = protein.domains[0]

    with pytest.raises(DomainException):
        domain.get_track_values('syms')

    with pytest.raises(DomainException):
        domain.get_track_symbols('vals')


# --- bug 9 -----------------------------------------------------------------
def test_bug09_site_with_no_symbol_round_trips(tmp_path):
    P = _proteome()
    P.protein('P0').add_site(4, 'phos', symbol=None, value=None)
    P.protein('P0').add_site(6, 'phos', symbol='P', value=1.5)

    filename = os.path.join(str(tmp_path), 'sites.tsv')
    si_sites.write_sites(P, filename)

    P2 = _proteome()
    si_sites.add_sites_from_file(P2, filename)

    sites = {s.position: s for s in P2.protein('P0').sites}

    assert sites[4].symbol is None
    assert sites[4].value is None
    assert sites[6].symbol == 'P'
    assert sites[6].value == 1.5
