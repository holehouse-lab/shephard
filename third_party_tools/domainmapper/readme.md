# dommap2shephard
##### Last updated 2022-12-01

## About
`dommap2shephard.py` is a Python script for converting output files from [DomainMapper]() to a SHEPHARD-compliant [Domains file](https://shephard.readthedocs.io/en/latest/shephard_file_types.html#domain-files). DomainMapper is a Python tool developed by the [Fried lab]() to robustly identify distinct domains, specifically annotating those with complex topologies. The resulting domains are well-suited for downstream analysis using SHEPHARD. As such, we developed `dommap2shephard.py` to make these annotations easy to obtain.

## Usage

Using `dommap2shephard.py` involves two steps:

#### Step 1: Using DomainMapper to assign domains
Instructions for download and installing DomainMapper are provided on the DomainMapper webpage. Briefly, this involves using the [ECOD HMMs](http://prodata.swmed.edu/ecod/complete/distribution) database in conjunction with [hmmer](http://hmmer.org/download.html) first to assign domains from a FASTA file. Having initially assigned those domains, DomainMapper define specific continuous topologies enabling single domains that are non-continuous or contain insertional domains to be rationally assigned, formatting domains such that they can be easily converted into a SHEPHARD-complinat format.

The output from DomainMapper is a text file with a set of comment lines followed by lines with nine tab-separated columns. This is the file that `dommap2shephard` reads in.

#### Step 2: Using dommap2shephard to build a SHEPHARD-compliant file
`dommap2shephard.py` should be run as a command-line script as follows:

	python dommap2shephard.py -f <input file> -o <output filename> [--uniprot]

This command should be executed from the directory where `dommap2shephard.py` is found and should be pretty fast. For example,  analyzing all domains annotated across the human proteome took < 3 seconds on a MacBook Pro (2022).

Here:

* The `<input file>` is file that DomainMapper generated 
* The `<output file>` is the name of the file you want dommap2shephard to write to
* The `--uniprot` flag is a flag which, if provided, means DomainMapper files generated from UniProt derived FASTA files use the UniProt ID as the protein identifier, otherwise the identifier used is extracted by DomainMapper.

#### A note on software dependencies
DomainMapper is a third-party tool – that is, it's a software tool that was not developed by the Holehouse lab. As such, we provide `dommap2shephard.py` as a convenient script to help assign SHEPHARD domains from DomainMapper-derived files, but do not formally couple the SHEPHARD code base to DomainMapper. That is, DomainMapper is not a dependency of SHEPHARD. 

If DomainMapper is updated or changes, we'll do our best to ensure `dommap2shephard.py` is updated, but this loose coupling avoids a scenario where if, DomainMapper breaks, SHEPHARD also becomes inoperable. This approach is a deliberate design principle that SHEPHARD follows to avoid chained dependencies and ensure code based on SHEPHARD will remain functional regardless of the state of other software tools.

### Output domain file format (as read by SHEPHARD)

As with all SHEPHARD domain files, the first few positions define the required columns in the [Domains file](https://shephard.readthedocs.io/en/latest/shephard_file_types.html#domain-files):

1. ID - The ID of the protein associated with the domain
2. start – The position where the domain starts (Where first position in aprotein = 1)
3. end - The position where the domain ends (inclusive)
4. type - SHEPHARD enables different types of domains, so this reflects the fact this domain came from DomainMapper (as opposed the type of domain it is in terms of structural annotation). Note that *all* domains will be of type `DomainMapper_X` where X is an ID that is unique to the domain in that protein. This enables multiple domains of the same type not to clash.

After this, we move into the Domain attributes. This is where most of the useful information is found, and attributes are simple key-value pairs that enable Domain objects in SHEPHARD to annotate with a specific type of data. The attributes and columns for dommap2shephard-generated files are shown below


5. `domain_type` - ECOD Family (F) group - a controlled vocabulary of Family (F-group) domains names.
6. `non_contiguous` - Flag which defines if the domain contains non-contigous segments or not. Will be True or False.
7. `insertional` - Flag which annotates if the domain is an insertional domain or not (e.g. a domain inserted within a loop of another domain). Will be True or False.
8. `circularly_permuted` – Flag which annotates if the domain is circularly permuted - a domain where structure is preserved but the N- and C-termini start at different positions along the domain sequence.
9. `e_value` – value (float) that gives the expected likelihood of the sequence in question being annotated by the domain of annotation by random chance (smaller means more confident the domain annotation is correct).
10. `ECOD_F_ID` – The ECOD F group ID; this is is a set of numbers separated by periods (e.g. `304.9.1.30`) that provide a controlled identifier for the domain.
11. `ECOD_T_group` – ECOD T group name - string defining a protein group based on topology name.
12. `ECOD_X_group` – ECOD X group name - string defining a protein group based on similar structure but lack convincing evidence for bona fide homology
13. `ECOD_architecture` – ECOD architecture string
14. `domain_index` – For non_contigous domains, this string is of format X_of_Y where X is the current index and Y is the total number of non_contigous domains. For 'normal' domains that are contigous this is just 1_of_1
15. `domain_count` - The Y described in the domain index above (i.e. how many separate records parts are associated with the underlying biological domain that this SHEPHARD domain object is associated with). Will be 1 most of the time, but for non-contiguous domains will report the number of non-contigous segments. 
16. `domain_unique_string` - A unique string that defines this domain for easy external referencing. The string has the format: `<UID>_domain_<domain_count>_<domain_type>_idx_<start>_<end>_<domain_index>`, e.g. this might look like `P04147_domain_1_RRM_1_idx_38_115_1_of_1`.
17. `from_domain_mapper` : Flag which is always set to True defining this domain as from DomainMapper. The main reason to include this as an attribute is it enables SHEPHARD analysis to filter DomainMapper domains using the following pattern

		if domain.attribute('from_domain_mapper',safe=False) is 'True':
		    # do something

## Precomputed domains
To make life easy, we have added precomputed ECOD domains using the dommap2shephard pipeline for human and yeast proteomes, with the remaining proteomes coming soon. All data are shared on our [SHEPHARD-data GitHub page](https://github.com/holehouse-lab/shephard-data).

## Example analysis
To see an example of the type of analysis one could do, see the `demo_domainmapper.ipynb` notebook in this directory.

## Tests

To run the dommap2shephard tests run

	pytest -v test_dommap2shephard.py
	
Note this does obviously require `pytest` to be installed, which if it's not can be installed using pip:

	pip install pytest	
	
Note that the tests should be updated if/when DomainMapper is updated, because they use a precomputed output file from DomainMapper by default.	
## Help
If you encounter any problems, please raise an issue on the [SHEPHARD GitHub page](https://github.com/holehouse-lab/shephard/issues). Alternatively contact Alex at alex dot holehouse 

## Changelog

#### Version 1.0 (Nov 2022)
Initial release.
	