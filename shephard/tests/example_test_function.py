# must import pytest
import pytest

# import additional modules needed for this code
import shephard
from shephard.apis import uniprot  
from shephard import interfaces

# import the domain_tools module which we're going to test
from shephard.tools import domain_tools

def test_first_test():

    # this line uses shephards internal get_data function to build the full
    # path to the directory where our test data is
    test_data_dir = shephard.get_data('test_data')

    # define filenames - note we use the names we defined above
    fasta_file = '%s/%s' % (test_data_dir, 'testset_1.fasta')
    domain_file = '%s/%s' % (test_data_dir, 'TS1_domains_idr.tsv')

    # create a new Proteome object and annotate with the domains
    # in the domain_file
    P = uniprot.uniprot_fasta_to_proteome(fasta_file)
    interfaces.si_domains.add_domains_from_file(P, domain_file)

    P.protein('O00401').domains

    assert domain_tools.calculate_domain_length(P.protein('O00401').domains[0]) == 20

