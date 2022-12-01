import pytest

import dommap2shephard 


# ........................................................................................
#
def test_parse_file_uniprot():

    X = dommap2shephard.read_domainmapper_file('test_data.out', uniprot=True)


    # run some tests for each
    dr = X['P04147'][0]
    
    assert dr.circularly_permuted == False
    assert dr.non_contiguous == False
    assert dr.uid == 'P04147'
    assert dr.domain_type == 'RRM_1'
    assert dr.shprd_domain_type == 'DomainMapper_1'
    assert int(dr.start) == 38
    assert int(dr.end) == 115
    assert dr.insertional == False
    assert dr.evalue == '1.00e-25'
    assert dr.ecod_fid == '304.9.1.30'
    assert dr.ecod_tgroup == 'RNA-binding domain, RBD'
    assert dr.ecod_xgroup == 'Alpha-beta plaits'
    assert dr.ecod_arch == 'a+b two layers'
    assert dr.domain_index == 1
    assert dr.domain_max_index == 1
    assert dr.shprd_domain_count == 1

    assert dr.generate_shephard_domain_line() == 'P04147	38	115	DomainMapper_1	domain_type:RRM_1	non_contiguous:False	insertional:False	circularly_permuted:False	e_value:1.00e-25	ECOD_F_ID:304.9.1.30	ECOD_T_group:RNA-binding domain, RBD	ECOD_X_group:Alpha-beta plaits	ECOD_architecture:a+b two layers	domain_index:1_of_1	domain_count:1	domain_unique_string:P04147_domain_1_RRM_1_idx_38_115_1_of_1	from_domain_mapper:True\n'

    assert dr.unique_name == 'P04147_domain_1_RRM_1_idx_38_115_1_of_1'
    assert dr.index_string == '1_of_1'


# ........................................................................................
#
def test_parse_file_all():

    X = dommap2shephard.read_domainmapper_file('test_data.out')

    # run some tests for each
    dr = X['sp|P04147|PABP_YEAST'][0]
    
    assert dr.circularly_permuted == False
    assert dr.non_contiguous == False
    assert dr.uid == 'sp|P04147|PABP_YEAST'
    assert dr.domain_type == 'RRM_1'
    assert dr.shprd_domain_type == 'DomainMapper_1'
    assert int(dr.start) == 38
    assert int(dr.end) == 115
    assert dr.insertional == False
    assert dr.evalue == '1.00e-25'
    assert dr.ecod_fid == '304.9.1.30'
    assert dr.ecod_tgroup == 'RNA-binding domain, RBD'
    assert dr.ecod_xgroup == 'Alpha-beta plaits'
    assert dr.ecod_arch == 'a+b two layers'
    assert dr.domain_index == 1
    assert dr.domain_max_index == 1
    assert dr.shprd_domain_count == 1

    assert dr.unique_name == 'sp|P04147|PABP_YEAST_domain_1_RRM_1_idx_38_115_1_of_1'
    
    assert dr.generate_shephard_domain_line() == 'sp|P04147|PABP_YEAST	38	115	DomainMapper_1	domain_type:RRM_1	non_contiguous:False	insertional:False	circularly_permuted:False	e_value:1.00e-25	ECOD_F_ID:304.9.1.30	ECOD_T_group:RNA-binding domain, RBD	ECOD_X_group:Alpha-beta plaits	ECOD_architecture:a+b two layers	domain_index:1_of_1	domain_count:1	domain_unique_string:sp|P04147|PABP_YEAST_domain_1_RRM_1_idx_38_115_1_of_1	from_domain_mapper:True\n'


    assert dr.unique_name == 'sp|P04147|PABP_YEAST_domain_1_RRM_1_idx_38_115_1_of_1'
    assert dr.index_string == '1_of_1'

    
    
    
# ........................................................................................
#
def test_parse_file_NC_protein():

    X = dommap2shephard.read_domainmapper_file('test_data.out', uniprot=True)

    # run some tests for each
    dr = X['P00549'][0]
    
    assert dr.circularly_permuted == False
    assert dr.non_contiguous == True
    assert dr.uid == 'P00549'
    assert dr.domain_type == 'PK_1st'
    assert dr.shprd_domain_type == 'DomainMapper_1'
    assert int(dr.start) == 17
    assert int(dr.end) == 120
    assert dr.insertional == False
    assert dr.evalue == '4.90e-116'
    assert dr.ecod_fid == '2002.1.1.220'
    assert dr.ecod_tgroup == 'TIM barrels'
    assert dr.ecod_xgroup == 'TIM beta/alpha-barrel'
    assert dr.ecod_arch == 'a/b barrels'
    assert dr.domain_index == 1
    assert dr.domain_max_index == 2
    assert dr.shprd_domain_count == 1

    assert dr.generate_shephard_domain_line() == 'P00549	17	120	DomainMapper_1	domain_type:PK_1st	non_contiguous:True	insertional:False	circularly_permuted:False	e_value:4.90e-116	ECOD_F_ID:2002.1.1.220	ECOD_T_group:TIM barrels	ECOD_X_group:TIM beta/alpha-barrel	ECOD_architecture:a/b barrels	domain_index:1_of_2	domain_count:1	domain_unique_string:P00549_domain_1_PK_1st_idx_17_120_1_of_2	from_domain_mapper:True\n'

    assert dr.unique_name == 'P00549_domain_1_PK_1st_idx_17_120_1_of_2'
    assert dr.index_string == '1_of_2'


    # run some tests for each
    dr = X['P00549'][1]
    
    assert dr.circularly_permuted == False
    assert dr.non_contiguous == True
    assert dr.uid == 'P00549'
    assert dr.domain_type == 'PK_1st'
    assert dr.shprd_domain_type == 'DomainMapper_1'
    assert int(dr.start) == 185
    assert int(dr.end) == 360
    assert dr.insertional == False
    assert dr.evalue == '4.90e-116'
    assert dr.ecod_fid == '2002.1.1.220'
    assert dr.ecod_tgroup == 'TIM barrels'
    assert dr.ecod_xgroup == 'TIM beta/alpha-barrel'
    assert dr.ecod_arch == 'a/b barrels'
    assert dr.domain_index == 2
    assert dr.domain_max_index == 2
    assert dr.shprd_domain_count == 1

    assert dr.generate_shephard_domain_line() == 'P00549	185	360	DomainMapper_1	domain_type:PK_1st	non_contiguous:True	insertional:False	circularly_permuted:False	e_value:4.90e-116	ECOD_F_ID:2002.1.1.220	ECOD_T_group:TIM barrels	ECOD_X_group:TIM beta/alpha-barrel	ECOD_architecture:a/b barrels	domain_index:2_of_2	domain_count:1	domain_unique_string:P00549_domain_1_PK_1st_idx_185_360_2_of_2	from_domain_mapper:True\n'

    assert dr.unique_name == 'P00549_domain_1_PK_1st_idx_185_360_2_of_2'
    assert dr.index_string == '2_of_2'
    
