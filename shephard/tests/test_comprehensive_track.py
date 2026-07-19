"""
SHEPHARD comprehensive test suite : Track

Detailed tests covering the user-facing Track API:

  * construction (values / symbols / string-symbols / both-error /
    length-mismatch / empty / attribute-dict typing)
  * properties (name/track_type/values/symbols/protein)
  * values_region / symbols_region (range & single-position)
  * value() / symbol() (+ cross-type safe behaviour)
  * attributes (add/get/remove/safe)
  * dunders (__len__, __repr__)

Holehouse Lab - Washington University in St. Louis
"""

import pytest

from shephard.proteome import Proteome
from shephard.protein import Protein
from shephard.track import Track
from shephard.exceptions import TrackException

SEQ = 'ACDEFGHIKLMNPQRSTVWY'   # 20 residues


@pytest.fixture
def prot():
    P = Proteome([{'sequence': SEQ, 'name': 'demo',
                   'unique_ID': 'P1', 'attributes': {}}])
    return P.protein('P1')


# ---------------------------------------------------------------------------
# construction
# ---------------------------------------------------------------------------
class TestTrackConstruction:

    def test_values_track(self, prot):
        prot.add_track('v', values=list(range(20)))
        t = prot.track('v')
        assert t.track_type == 'values'
        assert t.values == [float(i) for i in range(20)]
        assert t.symbols is None

    def test_symbols_track_from_list(self, prot):
        prot.add_track('s', symbols=list('X' * 20))
        t = prot.track('s')
        assert t.track_type == 'symbols'
        assert t.symbols == ['X'] * 20
        assert t.values is None

    def test_symbols_track_from_string(self, prot):
        # regression: the string branch used a buggy `elif(symbols, str)`
        prot.add_track('s', symbols='ABCDEFGHIJKLMNOPQRST')
        t = prot.track('s')
        assert t.symbols == list('ABCDEFGHIJKLMNOPQRST')

    def test_values_and_symbols_both_raises(self, prot):
        with pytest.raises(TrackException):
            Track('bad', prot, values=[0.0] * 20, symbols=list('X' * 20))

    def test_empty_track_raises(self, prot):
        with pytest.raises(TrackException):
            Track('bad', prot, values=None, symbols=None)

    def test_values_length_mismatch_raises(self, prot):
        with pytest.raises(TrackException):
            Track('bad', prot, values=[1.0, 2.0])

    def test_symbols_length_mismatch_raises(self, prot):
        with pytest.raises(TrackException):
            Track('bad', prot, symbols=list('XY'))

    def test_non_numeric_values_raises(self, prot):
        with pytest.raises(TrackException):
            Track('bad', prot, values=['a'] * 20)

    def test_attribute_dictionary_typing(self, prot):
        t = Track('t', prot, values=[0.0] * 20,
                  attribute_dictionary={'k': 'v'})
        assert t.attribute('k') == 'v'
        # non-dict, non-None attribute dict -> TrackException
        # (regression: previously raised a bare NameError)
        with pytest.raises(TrackException):
            Track('bad', prot, values=[0.0] * 20,
                  attribute_dictionary='not-a-dict')


# ---------------------------------------------------------------------------
# properties
# ---------------------------------------------------------------------------
class TestTrackProperties:

    def test_name_and_protein(self, prot):
        prot.add_track('v', values=[0.0] * 20)
        t = prot.track('v')
        assert t.name == 'v'
        assert isinstance(t.protein, Protein)

    def test_values_track_symbol_is_none(self, prot):
        prot.add_track('v', values=[1.0] * 20)
        assert prot.track('v').symbols is None

    def test_zero_values_track_roundtrips(self, prot):
        # all-zero values must remain a values track
        prot.add_track('z', values=[0.0] * 20)
        t = prot.track('z')
        assert t.track_type == 'values'
        assert t.values == [0.0] * 20


# ---------------------------------------------------------------------------
# region access
# ---------------------------------------------------------------------------
class TestTrackRegions:

    def test_values_region_range(self, prot):
        prot.add_track('v', values=list(range(20)))
        t = prot.track('v')
        # 1-based inclusive
        assert t.values_region(1, 3) == [0.0, 1.0, 2.0]
        assert t.values_region(18, 20) == [17.0, 18.0, 19.0]

    def test_values_region_single(self, prot):
        prot.add_track('v', values=list(range(20)))
        assert prot.track('v').values_region(5) == 4.0

    def test_symbols_region_range_and_single(self, prot):
        prot.add_track('s', symbols=list('ABCDEFGHIJKLMNOPQRST'))
        t = prot.track('s')
        assert t.symbols_region(1, 3) == ['A', 'B', 'C']
        assert t.symbols_region(3) == 'C'

    def test_region_invalid_position_raises(self, prot):
        prot.add_track('v', values=[0.0] * 20)
        with pytest.raises(Exception):
            prot.track('v').values_region(0, 5)
        with pytest.raises(Exception):
            prot.track('v').values_region(1, 99)


# ---------------------------------------------------------------------------
# value() / symbol()
# ---------------------------------------------------------------------------
class TestTrackValueSymbol:

    def test_value(self, prot):
        prot.add_track('v', values=list(range(20)))
        assert prot.track('v').value(10) == 9.0

    def test_symbol(self, prot):
        prot.add_track('s', symbols=list('ABCDEFGHIJKLMNOPQRST'))
        assert prot.track('s').symbol(1) == 'A'

    def test_value_on_symbol_track_safe(self, prot):
        prot.add_track('s', symbols=list('X' * 20))
        t = prot.track('s')
        with pytest.raises(TrackException):
            t.value(1)
        assert t.value(1, safe=False) is None

    def test_symbol_on_value_track_safe(self, prot):
        prot.add_track('v', values=[1.0] * 20)
        t = prot.track('v')
        with pytest.raises(TrackException):
            t.symbol(1)
        assert t.symbol(1, safe=False) is None

    def test_value_invalid_position(self, prot):
        prot.add_track('v', values=[0.0] * 20)
        with pytest.raises(Exception):
            prot.track('v').value(0)


# ---------------------------------------------------------------------------
# attributes
# ---------------------------------------------------------------------------
class TestTrackAttributes:

    def test_add_get_duplicate_safe(self, prot):
        prot.add_track('v', values=[0.0] * 20)
        t = prot.track('v')
        t.add_attribute('k', 1)
        assert t.attribute('k') == 1
        with pytest.raises(TrackException):
            t.add_attribute('k', 2)
        t.add_attribute('k', 2, safe=False)
        assert t.attribute('k') == 2

    def test_missing_safe(self, prot):
        prot.add_track('v', values=[0.0] * 20)
        t = prot.track('v')
        with pytest.raises(TrackException):
            t.attribute('absent')
        assert t.attribute('absent', safe=False) is None

    def test_remove(self, prot):
        prot.add_track('v', values=[0.0] * 20)
        t = prot.track('v')
        t.add_attribute('k', 1)
        t.remove_attribute('k')
        assert 'k' not in t.attributes
        with pytest.raises(Exception):
            t.remove_attribute('k')
        t.remove_attribute('k', safe=False)


# ---------------------------------------------------------------------------
# dunders
# ---------------------------------------------------------------------------
class TestTrackDunders:

    def test_len(self, prot):
        prot.add_track('v', values=[0.0] * 20)
        assert len(prot.track('v')) == 20

    def test_repr(self, prot):
        prot.add_track('v', values=[0.0] * 20)
        assert 'Track' in repr(prot.track('v'))
