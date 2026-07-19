Getting Started
=================

The typical SHEPHARD workflow has three steps: (1) build a ``Proteome`` from a FASTA file, (2) layer on annotations from SHEPHARD files (or compute them on the fly), and (3) query the relationship between sequence and those annotations.

A minimal end-to-end example
-------------------------------

The following code reads a UniProt FASTA file, annotates it with domains, and then prints every protein that has a domain starting at its immediate N-terminus.

.. code-block:: python

    from shephard.apis import uniprot
    from shephard.interfaces import si_domains

    # read a UniProt-formatted FASTA file into a new Proteome
    P = uniprot.uniprot_fasta_to_proteome('fasta_file.fasta')

    # annotate every protein with domains from a SHEPHARD domains file
    si_domains.add_domains_from_file(P, 'domains_file.tsv')

    # iterate the Proteome -> Protein -> Domain hierarchy
    for protein in P:
        for domain in protein.domains:
            if domain.start == 1:
                print(f'Protein {protein} has an N-terminal domain: {domain}')

A few things to note:

* Iterating a ``Proteome`` (``for protein in P``) yields ``Protein`` objects directly.
* ``protein.domains`` returns a list of ``Domain`` objects (sorted by start position) — you iterate the objects directly; you do **not** index back through the protein.
* All positions use biology-style indexing (1-based, inclusive), so ``domain.start == 1`` really is the first residue.

Building a Proteome by hand
----------------------------

You do not need a file to get started — Proteins can be added programmatically:

.. code-block:: python

    from shephard import Proteome

    P = Proteome()
    P.add_protein('MAEPQRDGSALLK', 'my cool protein', 'test1')

    prot = P.protein('test1')
    prot.add_domain(1, 5, 'IDR')
    prot.add_site(4, 'phosphorylation', symbol='P', value=1.0)

    print(prot.get_sequence_region(1, 5))      # 'MAEPQ'
    print(prot.get_domains_by_type('IDR'))     # [ |Domain: IDR ... ]
    print(prot.site(4))                        # [ |Site: phosphorylation ... ]

For more worked examples, see :doc:`examples`.
