# UniProt API demos

### `uniprot_api_examples.ipynb`

A walk through SHEPHARD's UniProt REST API functionality — annotating proteins and proteomes with domains, sites, and tracks pulled live from UniProt.

The notebook covers:

* discovering the available annotation groups
* annotating a single protein (`annotate_protein_with_uniprot`)
* how UniProt features map onto SHEPHARD Domains and Sites, including the disulfide bond / cross-link case
* the metadata and evidence codes attached to every annotation
* fetching only the annotation classes you need
* experimental structures and per-residue Tracks (with a figure)
* batch-annotating a whole proteome (`annotate_proteome_with_uniprot`)
* error handling, and why the sequence check matters
* saving annotations to SHEPHARD files for reproducibility
* a worked example asking whether PTMs are enriched in disordered regions

The notebook is saved with its outputs, so it can be read without running anything. To run it you need network access — every example makes live calls to UniProt. Note that the exact annotation counts will drift over time as UniProt is updated.

Reference documentation for these functions is in the [SHEPHARD API docs](https://shephard.readthedocs.io/en/latest/apis.html).
