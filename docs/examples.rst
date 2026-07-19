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


Example 4: annotating a proteome directly from UniProt
--------------------------------------------------------

Annotations do not have to come from a file you already have. The :code:`uniprot` API can fetch them live from UniProt and convert them into Domains, Sites, and Tracks (see :doc:`apis`). The example below is fully runnable — it uses the small test proteome bundled with SHEPHARD, so it needs nothing but a network connection.

.. code-block:: python

    import shephard
    from shephard.apis import uniprot

    # the nine-protein test set that ships with SHEPHARD
    test_data = shephard.get_data('test_data')
    P = uniprot.uniprot_fasta_to_proteome(f'{test_data}/testset_1.fasta')

    # fetch annotations for all nine proteins in a single batched query.
    # safe=False means a protein we cannot annotate is reported rather than
    # aborting the whole run
    summary = uniprot.annotate_proteome_with_uniprot(P, safe=False, verbose=False)

    for unique_ID in sorted(summary):
        result = summary[unique_ID]

        if 'error' in result:
            print(f"{unique_ID}: not annotated ({result['error']})")
        else:
            print(f"{unique_ID}: {result['domains']} domains, {result['sites']} sites")

Running this prints something close to::

    O00401: 12 domains, 6 sites
    O00470: 7 domains, 0 sites
    O00472: 10 domains, 2 sites
    ...
    O14786: 9 domains, 24 sites
    Q9UJX3: not annotated (sequence mismatch)

(the exact counts will drift as UniProt is updated)

That last line is the sequence check doing its job. The bundled FASTA file predates a UniProt revision of ``Q9UJX3``, which went from 599 to 565 residues. Because UniProt annotations are positional, applying today's annotations to yesterday's sequence would have produced plausible-looking but incorrect Domains, so SHEPHARD refuses and tells you why.

Every annotation carries its provenance, so you can always ask where something came from:

.. code-block:: python

    for domain in P.protein('O00401').domains:
        print(domain.start, domain.end,
              domain.attribute('uniprot_feature_type'),
              domain.attribute('description'),
              domain.attribute('evidence'))

Because UniProt is a moving target — and because it is a free public resource — it is usually best to fetch annotations once and then save them, rather than re-querying on every analysis run:

.. code-block:: python

    from shephard.interfaces import si_domains, si_sites

    si_domains.write_domains(P, 'uniprot_domains.tsv')
    si_sites.write_sites(P, 'uniprot_sites.tsv')

Subsequent analyses can then read those files back, which is faster and — importantly — reproducible.


Example 5: are UniProt PTM sites enriched in predicted IDRs?
--------------------------------------------------------------

This is Example 2 again, but with both halves of the annotation obtained programmatically: PTM Sites come from UniProt, and IDR Domains from `metapredict <https://metapredict.readthedocs.io/>`_. Nothing here is read from disk, and the two annotation sources know nothing about each other — SHEPHARD is what makes them comparable.

.. code-block:: python

    import shephard
    from shephard.apis import uniprot, metapredict_api

    test_data = shephard.get_data('test_data')
    P = uniprot.uniprot_fasta_to_proteome(f'{test_data}/testset_1.fasta')

    # (1) modified residues (phosphorylation, acetylation, ...) from UniProt.
    # Everything else is switched off so we download only what we need
    uniprot.annotate_proteome_with_uniprot(P,
                                           modified_residues=True,
                                           domains=False, regions=False,
                                           motifs=False, repeats=False,
                                           zinc_fingers=False, coiled_coils=False,
                                           compositional_bias=False,
                                           dna_binding=False, transmembrane=False,
                                           binding_sites=False, active_sites=False,
                                           other_sites=False, glycosylation=False,
                                           lipidation=False, disulfide_bonds=False,
                                           cross_links=False,
                                           safe=False, verbose=False)

    # (2) IDRs from metapredict
    metapredict_api.annotate_proteome_with_disordered_domains(P, safe=False)

    # (3) ask the question
    ptms_in_idrs = 0
    idr_residues = 0

    for protein in P:
        for domain in protein.get_domains_by_type('IDR'):
            idr_residues += len(domain)
            ptms_in_idrs += len(domain.get_sites_by_type('uniprot_modified_residue',
                                                         return_list=True))

    total_ptms = len(P.get_sites_by_type('uniprot_modified_residue'))
    total_residues = sum([len(protein) for protein in P])

    print(f'{ptms_in_idrs}/{total_ptms} PTMs fall in IDRs')
    print(f'IDRs are {idr_residues/total_residues:.1%} of all residues')

On a nine-protein test set the numbers are not meaningful, but the code is identical for a whole proteome — which is the point.


Example 6: which parts of a protein have no experimental structure?
---------------------------------------------------------------------

UniProt records which regions of a protein are covered by an experimentally-determined structure via its PDB cross-references. SHEPHARD can turn this into a per-residue Track, which makes "what is structurally uncharacterized?" a straightforward question — and one that pairs naturally with disorder prediction.

.. code-block:: python

    import shephard
    from shephard.apis import uniprot

    test_data = shephard.get_data('test_data')
    P = uniprot.uniprot_fasta_to_proteome(f'{test_data}/testset_1.fasta')

    # experimental_structures adds a Domain per resolved region (carrying the
    # PDB ID, method, and resolution), while structure_coverage_track adds a
    # 1/0 per-residue Track
    uniprot.annotate_proteome_with_uniprot(P,
                                           experimental_structures=True,
                                           structure_coverage_track=True,
                                           safe=False, verbose=False)

    for protein in P:

        track = protein.track('uniprot_structure_coverage', safe=False)
        if track is None:
            continue

        covered = sum(track.values)
        print(f'{protein.unique_ID}: {covered/len(protein):.0%} structurally covered')

        # what structures contribute, and how good are they?
        for structure in protein.get_domains_by_type('uniprot_experimental_structure'):
            print(f"    {structure.attribute('pdb_id')} "
                  f"({structure.attribute('method')}, {structure.attribute('resolution')}) "
                  f"covers {structure.start}-{structure.end}")

Because the coverage Track is just a values Track like any other, it can be fed straight into the ``domain_tools`` machinery shown in Example 3 — for instance to define "structurally dark" Domains as contiguous stretches of at least 30 residues with no structural coverage:

.. code-block:: python

    from shephard.tools import domain_tools
    from shephard.interfaces import si_domains

    def uncovered(values):
        return [1 if v == 0 else 0 for v in values]

    dark_regions = domain_tools.build_domains_from_track_values(
        P, 'uniprot_structure_coverage', uncovered, 'structurally_dark',
        minimum_region_size=30, verbose=False)

    si_domains.add_domains_from_dictionary(P, dark_regions)
