"""
SHEPHARD comprehensive test suite : apis.fasta + apis.uniprot

Detailed tests covering the user-facing FASTA / UniProt IO API:

  * fasta.fasta_to_proteome (auto IDs / header IDs / custom ID & attributes /
    add-to-existing / invalid sequence handling)
  * fasta.proteome_to_fasta + fasta.shephard_fasta_to_proteome round trip
  * uniprot.uniprot_accession_from_line
  * uniprot.uniprot_fasta_to_proteome + uniprot_proteome_to_fasta round trip

(metapredict / albatross APIs require ML dependencies and are covered by
their own dedicated test modules.)

Holehouse Lab - Washington University in St. Louis
"""

import pytest

import shephard
from shephard.proteome import Proteome
from shephard.apis import fasta, uniprot
from shephard.exceptions import UtilitiesException

DATA = shephard.get_data('test_data')
UNIPROT_FASTA = f'{DATA}/testset_1.fasta'
GENERIC_FASTA = f'{DATA}/seqs.fasta'


# ===========================================================================
# fasta.fasta_to_proteome
# ===========================================================================
class TestFastaToProteome:

    def test_default_auto_unique_ids(self):
        P = fasta.fasta_to_proteome(GENERIC_FASTA)
        assert isinstance(P, Proteome)
        assert len(P) == 10
        # name is the full header
        for prot in P:
            assert len(prot.sequence) > 0
            assert len(prot.name) > 0

    def test_use_header_as_unique_ID(self):
        P = fasta.fasta_to_proteome(GENERIC_FASTA,
                                    use_header_as_unique_ID=True)
        assert len(P) == 10
        for uid in P.proteins:
            assert P.protein(uid).name == uid

    def test_custom_build_unique_ID(self):
        def grab_acc(header):
            return header.split('|')[1]
        P = fasta.fasta_to_proteome(GENERIC_FASTA, build_unique_ID=grab_acc)
        assert 'P41956' in P

    def test_build_attributes(self):
        def attrs(header):
            return {'header_len': len(header)}
        P = fasta.fasta_to_proteome(GENERIC_FASTA, build_attributes=attrs)
        any_prot = P[0]
        assert 'header_len' in any_prot.attributes

    def test_add_to_existing_proteome(self):
        P = fasta.fasta_to_proteome(GENERIC_FASTA, use_header_as_unique_ID=True)
        n0 = len(P)
        fasta.fasta_to_proteome(UNIPROT_FASTA, proteome=P,
                                use_header_as_unique_ID=True)
        assert len(P) > n0

    def test_header_and_parsed_id_conflict_raises(self):
        with pytest.raises(Exception):
            fasta.fasta_to_proteome(
                GENERIC_FASTA,
                use_header_as_unique_ID=True,
                build_unique_ID=lambda h: h.split('|')[1])

    def test_invalid_sequence_action_fail(self, tmp_path):
        bad = tmp_path / 'bad.fasta'
        bad.write_text(">h1\nACDEFG123\n")
        with pytest.raises(Exception):
            fasta.fasta_to_proteome(str(bad),
                                    invalid_sequence_action='fail')

    def test_invalid_sequence_action_convert(self, tmp_path):
        bad = tmp_path / 'bad.fasta'
        bad.write_text(">h1\nACDEFBXZ\n")
        P = fasta.fasta_to_proteome(str(bad),
                                    invalid_sequence_action='convert')
        assert P[0].check_sequence_is_valid() is True


# ===========================================================================
# fasta round trip (SHEPHARD-format)
# ===========================================================================
class TestFastaRoundTrip:

    def test_proteome_to_fasta_and_back(self, tmp_path):
        P = uniprot.uniprot_fasta_to_proteome(UNIPROT_FASTA)
        out = str(tmp_path / 'out.fasta')
        fasta.proteome_to_fasta(out, P)

        Q = fasta.shephard_fasta_to_proteome(out)
        assert len(Q) == len(P)
        for uid in P.proteins:
            assert uid in Q
            assert Q.protein(uid).sequence == P.protein(uid).sequence

    def test_round_trip_with_attributes(self, tmp_path):
        P = uniprot.uniprot_fasta_to_proteome(UNIPROT_FASTA)
        P.protein('O00401').add_attribute('gene', 'WASL')
        out = str(tmp_path / 'attr.fasta')
        fasta.proteome_to_fasta(out, P, include_attributes_in_header=True)

        Q = fasta.shephard_fasta_to_proteome(out)
        assert Q.protein('O00401').attribute('gene') == 'WASL'


# ===========================================================================
# uniprot
# ===========================================================================
class TestUniprot:

    def test_accession_from_line(self):
        line = '>sp|O00401|WASL_HUMAN Neural Wiskott-Aldrich ...'
        assert uniprot.uniprot_accession_from_line(line) == 'O00401'

    def test_accession_from_line_bad_raises(self):
        with pytest.raises(UtilitiesException):
            uniprot.uniprot_accession_from_line('no-pipes-here')

    def test_uniprot_fasta_to_proteome(self):
        P = uniprot.uniprot_fasta_to_proteome(UNIPROT_FASTA)
        assert len(P) == 9
        assert 'O00401' in P
        assert 'O14786' in P

    def test_uniprot_round_trip(self, tmp_path):
        P = uniprot.uniprot_fasta_to_proteome(UNIPROT_FASTA)
        out = str(tmp_path / 'uni.fasta')
        uniprot.uniprot_proteome_to_fasta(out, P)
        Q = uniprot.uniprot_fasta_to_proteome(out)
        assert len(Q) == len(P)
        for uid in P.proteins:
            assert Q.protein(uid).sequence == P.protein(uid).sequence

    def test_uniprot_add_to_existing(self):
        P = uniprot.uniprot_fasta_to_proteome(UNIPROT_FASTA)
        with pytest.raises(Exception):
            # re-reading the same file (duplicate UIDs) must raise by default
            uniprot.uniprot_fasta_to_proteome(UNIPROT_FASTA, proteome=P)
        # force_overwrite makes it succeed
        uniprot.uniprot_fasta_to_proteome(UNIPROT_FASTA, proteome=P,
                                          force_overwrite=True)
        assert len(P) == 9
