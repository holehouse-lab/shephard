import shephard
from shephard.interfaces import si_uniprot, si_domains, si_tracks, si_sites
from shephard import quickstart
import pytest
import sys


TS1_FILE=['testset_1.fasta', 'TS1_domains_idr.tsv', 'TS1_domains_pscore.tsv', 'TS1_sites.tsv', 'TS1_tracks_pscore.tsv']
test_data_dir = shephard.get_data('test_data')

def build_proteome(fn):
    return quickstart.quickstart('%s/%s' % (test_data_dir,fn),  extract_unique_ID=si_uniprot.extract_unique_ID_uniprot)


p = build_proteome(TS1_FILE[0])
si_domains.read_in_domains(p, '%s/%s' %(test_data_dir, TS1_FILE[1]))
si_domains.read_in_domains(p, '%s/%s' %(test_data_dir, TS1_FILE[2]))
si_sites.read_in_sites(p, '%s/%s' %(test_data_dir, TS1_FILE[3]))
si_tracks.read_in_tracks(p, '%s/%s' %(test_data_dir, TS1_FILE[4]))
