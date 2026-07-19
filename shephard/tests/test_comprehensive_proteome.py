"""
SHEPHARD comprehensive test suite : Proteome

Detailed tests covering the user-facing Proteome API:

  * construction (dict list / Protein list / attributes / force_overwrite / empty)
  * protein access (.protein, .proteins, check_unique_ID, safe behaviour)
  * add / remove protein(s)
  * proteome-level attributes
  * aggregate domain / site / track accessors and unique-type bookkeeping
  * dunder methods (__len__, __iter__, __contains__, __getitem__, __repr__)

Holehouse Lab - Washington University in St. Louis
"""

import pytest

import shephard
from shephard.proteome import Proteome
from shephard.protein import Protein
from shephard.exceptions import ProteomeException

test_data_dir = shephard.get_data('test_data')


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _protein_dict(uid, seq='ACDEFGHIKLMNPQRSTVWY', name=None, attributes=None):
    return {'sequence': seq,
            'name': name if name is not None else f'protein_{uid}',
            'unique_ID': uid,
            'attributes': {} if attributes is None else attributes}


def _basic_proteome(n=3):
    return Proteome([_protein_dict(f'U{i}', seq='A' * (10 + i)) for i in range(n)])


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------
class TestProteomeConstruction:

    def test_empty_constructions(self):
        assert len(Proteome()) == 0
        assert len(Proteome(None)) == 0
        assert len(Proteome([])) == 0

    def test_from_dict_list(self):
        P = Proteome([_protein_dict('A1'), _protein_dict('A2')])
        assert len(P) == 2
        assert set(P.proteins) == {'A1', 'A2'}

    def test_numeric_fields_cast_to_str(self):
        P = Proteome([{'sequence': 'ACDE', 'name': 123, 'unique_ID': 456,
                       'attributes': {}}])
        prot = P.protein('456')
        assert prot.unique_ID == '456'
        assert prot.name == '123'

    def test_from_protein_object_list_is_deep_copy(self):
        src = _basic_proteome(2)
        src.protein('U0').add_attribute('k', {'nested': [1, 2]})
        src.protein('U0').add_domain(1, 5, 'dom')
        src.protein('U0').add_site(3, 'phos', symbol='P', value=1.0)
        src.protein('U0').add_track('t', values=[0.0] * 10)

        dest = Proteome([src.protein(uid) for uid in src.proteins])
        assert len(dest) == 2

        # mutating the copy must not touch the original
        dest.protein('U0').attribute('k')['nested'].append(99)
        assert src.protein('U0').attribute('k')['nested'] == [1, 2]
        assert len(dest.protein('U0').domains) == 1
        assert len(dest.protein('U0').sites) == 1
        assert dest.protein('U0').track('t').values == [0.0] * 10

    def test_duplicate_uid_raises(self):
        with pytest.raises(ProteomeException):
            Proteome([_protein_dict('DUP'), _protein_dict('DUP')])

    def test_duplicate_uid_force_overwrite(self):
        P = Proteome()
        P.add_proteins([_protein_dict('X', seq='AAAA')])
        P.add_proteins([_protein_dict('X', seq='CCCCCC')], force_overwrite=True)
        assert P.protein('X').sequence == 'CCCCCC'

    def test_missing_key_raises(self):
        with pytest.raises(ProteomeException):
            Proteome([{'sequence': 'ACDE', 'name': 'n'}])  # no unique_ID

    def test_proteome_attributes_argument(self):
        P = Proteome([_protein_dict('A1')], attributes={'organism': 'human'})
        assert P.attribute('organism') == 'human'

    def test_bad_attributes_argument_raises(self):
        with pytest.raises(ProteomeException):
            Proteome([_protein_dict('A1')], attributes='not-a-dict')

    def test_mixed_type_input_list_raises(self):
        P = _basic_proteome(1)
        with pytest.raises(ProteomeException):
            Proteome([_protein_dict('Z'), P.protein('U0')])


# ---------------------------------------------------------------------------
# Protein access
# ---------------------------------------------------------------------------
class TestProteinAccess:

    def test_protein_safe_true_raises(self):
        P = _basic_proteome(1)
        with pytest.raises(ProteomeException):
            P.protein('NOPE')

    def test_protein_safe_false_returns_none(self):
        P = _basic_proteome(1)
        assert P.protein('NOPE', safe=False) is None

    def test_protein_uid_is_stringified(self):
        P = Proteome([_protein_dict('100')])
        assert P.protein(100).unique_ID == '100'

    def test_check_unique_ID(self):
        P = _basic_proteome(2)
        assert P.check_unique_ID('U0') is True
        assert P.check_unique_ID('absent') is False

    def test_proteins_returns_id_list(self):
        P = _basic_proteome(3)
        ids = P.proteins
        assert isinstance(ids, list)
        assert sorted(ids) == ['U0', 'U1', 'U2']


# ---------------------------------------------------------------------------
# add / remove proteins
# ---------------------------------------------------------------------------
class TestAddRemoveProteins:

    def test_add_protein(self):
        P = Proteome()
        P.add_protein('ACDEK', 'myprot', 'P1', attributes={'a': 1})
        assert len(P) == 1
        assert P.protein('P1').sequence == 'ACDEK'
        assert P.protein('P1').attribute('a') == 1

    def test_add_protein_duplicate_raises(self):
        P = Proteome()
        P.add_protein('AAAA', 'n', 'P1')
        with pytest.raises(ProteomeException):
            P.add_protein('CCCC', 'n', 'P1')

    def test_add_protein_force_overwrite(self):
        P = Proteome()
        P.add_protein('AAAA', 'n', 'P1')
        P.add_protein('CCCC', 'n', 'P1', force_overwrite=True)
        assert P.protein('P1').sequence == 'CCCC'

    def test_add_proteins_empty_list_is_noop(self):
        P = _basic_proteome(2)
        P.add_proteins([])               # regression: must not IndexError
        assert len(P) == 2

    def test_remove_protein(self):
        P = _basic_proteome(3)
        P.remove_protein('U1')
        assert len(P) == 2
        assert 'U1' not in P.proteins

    def test_remove_protein_missing_safe(self):
        P = _basic_proteome(1)
        with pytest.raises(ProteomeException):
            P.remove_protein('absent')
        P.remove_protein('absent', safe=False)  # no raise

    def test_remove_proteins_bulk(self):
        P = _basic_proteome(4)
        P.remove_proteins(['U0', 'U2'])
        assert sorted(P.proteins) == ['U1', 'U3']


# ---------------------------------------------------------------------------
# Proteome-level attributes
# ---------------------------------------------------------------------------
class TestProteomeAttributes:

    def test_add_get_attribute(self):
        P = _basic_proteome(1)
        P.add_attribute('study', 'demo')
        assert 'study' in P.attributes
        assert P.attribute('study') == 'demo'

    def test_add_attribute_safe_duplicate_raises(self):
        P = _basic_proteome(1)
        P.add_attribute('k', 1)
        with pytest.raises(ProteomeException):
            P.add_attribute('k', 2)

    def test_add_attribute_overwrite(self):
        P = _basic_proteome(1)
        P.add_attribute('k', 1)
        P.add_attribute('k', 2, safe=False)
        assert P.attribute('k') == 2

    def test_attribute_missing_safe(self):
        P = _basic_proteome(1)
        with pytest.raises(ProteomeException):
            P.attribute('absent')
        assert P.attribute('absent', safe=False) is None

    def test_remove_attribute(self):
        P = _basic_proteome(1)
        P.add_attribute('k', 1)
        P.remove_attribute('k')
        assert 'k' not in P.attributes
        with pytest.raises(ProteomeException):
            P.remove_attribute('k')
        P.remove_attribute('k', safe=False)


# ---------------------------------------------------------------------------
# aggregate accessors + unique-type bookkeeping
# ---------------------------------------------------------------------------
class TestAggregateAccessors:

    def test_domains_aggregate_and_unique_types(self):
        P = _basic_proteome(2)
        P.protein('U0').add_domain(1, 4, 'IDR')
        P.protein('U0').add_domain(5, 8, 'FD')
        P.protein('U1').add_domain(1, 4, 'IDR')
        assert len(P.domains) == 3
        assert sorted(P.unique_domain_types) == ['FD', 'IDR']
        assert len(P.get_domains_by_type('IDR')) == 2
        assert P.get_domains_by_type('ID', perfect_match=False)

    def test_unique_domain_types_decrement_on_remove(self):
        P = _basic_proteome(1)
        prot = P.protein('U0')
        prot.add_domain(1, 4, 'IDR')
        assert 'IDR' in P.unique_domain_types
        prot.remove_domain(prot.get_domains_by_type('IDR')[0])
        assert 'IDR' not in P.unique_domain_types

    def test_sites_aggregate_and_unique_types(self):
        P = _basic_proteome(2)
        P.protein('U0').add_site(2, 'phos', value=1)
        P.protein('U1').add_site(3, 'acet', value=1)
        P.protein('U1').add_site(4, 'phos', value=1)
        assert len(P.sites) == 3
        assert sorted(P.unique_site_types) == ['acet', 'phos']
        assert len(P.get_sites_by_type('phos')) == 2
        assert len(P.get_sites_by_type(['phos', 'acet'])) == 3

    def test_tracks_unique_names_and_type_map(self):
        P = _basic_proteome(2)
        P.protein('U0').add_track('vt', values=[0.0] * 10)
        P.protein('U1').add_track('vt', values=[1.0] * 11)
        P.protein('U0').add_track('st', symbols=['x'] * 10)
        assert sorted(P.unique_track_names) == ['st', 'vt']
        tmap = P.track_names_to_track_type
        assert tmap == {'vt': 'values', 'st': 'symbols'}
        # returned dict is a copy - mutating it must not corrupt internals
        tmap['vt'] = 'symbols'
        assert P.track_names_to_track_type['vt'] == 'values'

    def test_conflicting_track_type_same_name_raises(self):
        P = _basic_proteome(2)
        P.protein('U0').add_track('dup', values=[0.0] * 10)
        with pytest.raises(ProteomeException):
            P.protein('U1').add_track('dup', symbols=['a'] * 11)

    def test_tracks_aggregate_property(self):
        # Proteome.tracks parallels Proteome.domains / Proteome.sites
        P = _basic_proteome(2)
        P.protein('U0').add_track('vt', values=[0.0] * 10)
        P.protein('U0').add_track('st', symbols=['x'] * 10)
        P.protein('U1').add_track('vt', values=[1.0] * 11)
        all_tracks = P.tracks
        assert len(all_tracks) == 3
        # every element is a Track and maps back to its protein
        assert {t.name for t in all_tracks} == {'vt', 'st'}
        assert all(t.protein.unique_ID in ('U0', 'U1') for t in all_tracks)

    def test_get_tracks_by_name(self):
        P = _basic_proteome(3)
        P.protein('U0').add_track('vt', values=[0.0] * 10)
        P.protein('U1').add_track('vt', values=[1.0] * 11)
        P.protein('U2').add_track('other', values=[2.0] * 12)
        vt = P.get_tracks_by_name('vt')
        assert len(vt) == 2
        assert {t.protein.unique_ID for t in vt} == {'U0', 'U1'}
        assert P.get_tracks_by_name('absent') == []


# ---------------------------------------------------------------------------
# dunder methods
# ---------------------------------------------------------------------------
class TestProteomeDunders:

    def test_len(self):
        assert len(_basic_proteome(5)) == 5

    def test_iter_yields_protein_objects(self):
        P = _basic_proteome(3)
        objs = list(P)
        assert len(objs) == 3
        assert all(isinstance(o, Protein) for o in objs)

    def test_contains_str_and_protein(self):
        P = _basic_proteome(2)
        assert 'U0' in P
        assert 'absent' not in P
        assert P.protein('U1') in P

    def test_contains_other_type_returns_false(self):
        # __contains__ must return a bool, not implicitly None
        P = _basic_proteome(1)
        assert (12345 in P) is False
        assert (None in P) is False

    def test_getitem_int(self):
        P = _basic_proteome(3)
        assert isinstance(P[0], Protein)
        assert P[0].unique_ID == 'U0'

    def test_getitem_slice(self):
        P = _basic_proteome(4)
        sub = P[1:3]
        assert [p.unique_ID for p in sub] == ['U1', 'U2']

    def test_getitem_bad_key_raises(self):
        P = _basic_proteome(1)
        with pytest.raises(KeyError):
            P['not-an-int']
        with pytest.raises(KeyError):
            P[-1]

    def test_repr(self):
        assert 'Proteome' in repr(_basic_proteome(2))


# ---------------------------------------------------------------------------
# integration with the shared TS1 fixtures
# ---------------------------------------------------------------------------
class TestProteomeFixtureIntegration:

    def test_ts1_basic(self, TS1):
        assert len(TS1) == 9
        assert len(TS1.proteins) == 9

    def test_ts1_domains_and_sites(self, TS1_domains2_sites):
        assert len(TS1_domains2_sites.domains) > 0
        assert len(TS1_domains2_sites.sites) > 0
        assert 'IDR' in TS1_domains2_sites.unique_domain_types
