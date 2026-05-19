Examples
=================

Most examples are provided as stand-alone runnable notebooks on the `SHEPHARD colab repository <https://github.com/holehouse-lab/shephard-colab>`_, but a few representative snippets are included here. They illustrate the kind of integrative, proteome-wide questions SHEPHARD is designed to make easy to ask.

Example 1: domains at the N-terminus
--------------------------------------

Print the set of proteins that have a domain at their immediate N-terminus.

.. code-block:: python

    from shephard.apis import uniprot
    from shephard.interfaces import si_domains

    # read a UniProt-based FASTA file
    P = uniprot.uniprot_fasta_to_proteome('fasta_file.fasta')

    si_domains.add_domains_from_file(P, 'domains_file.tsv')

    for protein in P:
        for domain in protein.domains:
            if domain.start == 1:
                print(f'Protein {protein} has an N-terminal domain: {domain}')


Example 2: are PTMs enriched in IDRs?
---------------------------------------

A recurring observation in the literature is that intrinsically disordered regions (IDRs) are enriched for post-translational modification (PTM) sites. With a Proteome annotated with IDR ``Domains`` and PTM ``Sites``, asking this question is just a traversal of the hierarchy: for every IDR domain, count the sites that fall inside it (``domain.sites`` returns exactly the Sites within the domain boundaries).

.. code-block:: python

    from shephard.apis import uniprot
    from shephard.interfaces import si_domains, si_sites

    # load the proteome and annotate it with IDRs and PTM sites
    human = uniprot.uniprot_fasta_to_proteome('human_proteome.fasta')
    si_domains.add_domains_from_file(human, 'idrs.tsv')
    si_sites.add_sites_from_file(human, 'ptms.tsv')

    ptm_residues_in_idrs = 0
    idr_residues = 0

    # for every protein in the human proteome
    for protein in human:

        # for every domain in that protein
        for d in protein.domains:

            # if that domain is an IDR ...
            if d.domain_type == 'IDR':
                idr_residues += len(d)

                # ... count the PTM sites inside the IDR
                ptm_residues_in_idrs += len(d.sites)

    frac = ptm_residues_in_idrs / idr_residues
    print(f'{frac:.3%} of IDR residues carry a PTM site')

This is essentially the analysis that, after correcting for sequence composition and solvent accessibility, showed that many classes of PTM remain genuinely enriched in IDRs (Ginell *et al.*, *Bioinformatics*, 2023).


Example 3: building a Track on the fly and discretizing it into Domains
------------------------------------------------------------------------

Annotations do not have to come from a file. Here we compute a per-residue charge track from the sequence, then use the ``tools`` modules to turn contiguous high-signal stretches into Domains.

.. code-block:: python

    from shephard.apis import uniprot
    from shephard.tools import domain_tools
    from shephard.interfaces import si_domains

    P = uniprot.uniprot_fasta_to_proteome('fasta_file.fasta')

    # add an |charge| values track to every protein
    def abs_charge(seq):
        return [1.0 if r in 'KRED' else 0.0 for r in seq]

    for protein in P:
        protein.build_track_values_from_sequence('charged', abs_charge)

    # discretize the 'charged' track into 'charged_region' domains
    def binarize(values):
        return [1 if v > 0.5 else 0 for v in values]

    domain_dict = domain_tools.build_domains_from_track_values(
        P, 'charged', binarize, 'charged_region',
        minimum_region_size=15, gap_closure=2, verbose=False)

    si_domains.add_domains_from_dictionary(P, domain_dict)

The same analysis code can be re-pointed at a completely different dataset without modification — a core design goal of SHEPHARD.
