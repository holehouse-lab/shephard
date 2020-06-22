Examples
=================
Most examples will be provided as stand-alone runable code, but some simple examples are included here too

The following code would print the set of proteins that have a domain at their immediate N termini


.. code-block:: Python

    from shephard.apis import uniprot
    from shephard.interfaces import si_domains, si_sites

    # read in a UniProt based FASTA file
    P = uniprot.uniprot_fasta_to_proteome('fasta_file.fasta')

    si_domains.add_domains_from_file(proteome, 'domains_file.tsv')

    for protein in P:
        for domain_idx in protein.domains:
            domain = protein.domain(domain_idx)
            if domain.start ==  1:
                print('Protein %s has an N-terminal domain: %s' % (str(protein), str(domain)))

