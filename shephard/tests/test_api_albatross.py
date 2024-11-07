import shephard
from shephard.apis import albatross_api
from shephard.apis import uniprot
from shephard import interfaces
import numpy as np

rg_predictions = {}
rg_predictions['O00401'] = 74.3378
rg_predictions['O00470'] = 59.1429
rg_predictions['O00472'] = 82.6983
rg_predictions['O00499'] = 81.6329
rg_predictions['O00629'] = 74.7669
rg_predictions['O00712'] = 63.4530
rg_predictions['O00716'] = 74.3516
rg_predictions['O14786'] = 81.7541
rg_predictions['Q9UJX3'] = 77.4541

re_predictions = {}
re_predictions['O00401'] = 162.375
re_predictions['O00470'] = 132.6499
re_predictions['O00472'] = 184.6364
re_predictions['O00499'] = 189.2844
re_predictions['O00629'] = 172.3293
re_predictions['O00712'] = 139.4617
re_predictions['O00716'] = 177.6754
re_predictions['O14786'] = 188.8538
re_predictions['Q9UJX3'] = 178.6182

def test_proteome_wide_protein_preductions():

    # precomputed mean disorder we can compare against


    # build a proteome
    test_data_dir = shephard.get_data('test_data')
    P = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'))

    albatross_api.annotate_proteome_with_dimensions(P)

    for protein in P:
        assert np.isclose(protein.attribute('rg'), rg_predictions[protein.unique_ID])
        assert np.isclose(protein.attribute('re'), re_predictions[protein.unique_ID])


def test_proteome_wide_protein_preductions():

    # precomputed mean disorder we can compare against

    re_predictions = {}
    rg_predictions = {}

    
    re_predictions['pscore_domain_295_383'] = 83.0407
    re_predictions['pscore_domain_207_254'] = 46.1245
    
    rg_predictions['pscore_domain_295_383'] = 33.9210
    rg_predictions['pscore_domain_207_254'] = 18.8230


    # build a proteome
    test_data_dir = shephard.get_data('test_data')
    P = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'))
    interfaces.si_domains.add_domains_from_file(P, f"{test_data_dir}/TS1_domains_pscore.tsv")

    albatross_api.annotate_domains_with_dimensions(P, 'pscore_domain')
    for d in P.domains:
        assert np.isclose(d.attribute('rg'), rg_predictions[d.domain_name])
        assert np.isclose(d.attribute('re'), re_predictions[d.domain_name])

        
