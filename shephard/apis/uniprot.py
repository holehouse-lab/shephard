"""
uniprot.py

From the SHEPHARD package
Sequence-based Hierachical and Extendable Platform for High-throughput Analysis of Region of Disorder
Ginell & Holehouse, 2020-2022

Handles all I/O associated with uniprot-derived files.

"""

import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request

from shephard.exceptions import UtilitiesException, APIException, ShephardException
from shephard import exceptions as shephard_exceptions
import protfasta
from . import fasta

## ------------------------------------------------------------------------
##
def uniprot_accession_from_line(line):
    """
    Function that converts a header from a uniprot fasta file
    to extract the uniprot ID. This an example of the type of function
    that can be passed to quickstart using the extract_unique_id argument.

    This function assumes the uniprot-standard format for the header
    file has been maintained - i.e.

    >>> >xx|ACCESSION|xxxx

    where ACCESSION is the uniprot accession. 

    Parameters
    -----------

    line : string
        String where we expect the uniprot ID to be contained within two 'pipe' 
        characters ('|'). 

    Returns
    -----------
    string
        Returns the uniprot ID, although this is not formally validated. 
        However, assuming the string follows standard uniprot fasta header 
        conventions this should be true!

    """
    try:
        return line.split('|')[1].strip()
    except (IndexError, AttributeError):
        raise UtilitiesException(f'Unable to parse string [{line}] to identify uniprot ID')

        
        
## ------------------------------------------------------------------------
##
def uniprot_fasta_to_proteome(filename, 
                              proteome = None,
                              force_overwrite=False,
                              invalid_sequence_action='fail'):
                              
    """
    Stand alone function that allows the user to build a proteome from a 
    standard FASTA file downloaded from UniProt

    This function assumes the uniprot-standard format for the header
    file has been maintained - i.e.

    >>> >xx|ACCESSION|xxxx

    Where ACCESSION is the uniprot accession and will be used as the 
    unique_ID
    
    Parameters
    ------------

    filename : string
        Name of the FASTA file we're going to parse in. Note the protein 
        name will be defined as the full FASTA header for each entry.
        
    proteome : Proteome
        If a Proteome object is provided the FASTA file will be read and 
        added to the existing proteome, whereas if set to None a new 
        Proteome will be generated.

    force_overwrite : bool (default  = False)
        If this flag is set to true  and we encounter a unique_ID that is 
        already in the proteome the newer value overwrites the older one. 
        This is mostly useful if you are adding in a file with known 
        duplicate entries OR combining multiple FASTA files where you know 
        there's some duplications. Important - if we're building unique IDs
        based on numerical record indices then EVERY FASTA entry will be given 
        a unique_ID (meaning force_overwrite is irrelevant in this case).

    invalid_sequence_action : str (default = 'fail')
        Selector which defines the behaviour if a sequence with a non-
        standard amino acid is encountered. Valid options and their meaning
        are listed below:

            * ``ignore``  - invalid sequences are completely ignored

            * ``fail``    - invalid sequence cause parsing to fail and throw an exception
  
            * ``remove`` -  invalid sequences are removed

            * ``convert`` - invalid residues are converted to valid residues                            

            * ``convert-ignore`` - invalid sequences are converted to valid sequences and any remaining invalid residues are ignored.
    
    Returns 
    --------
    Proteome
        Returns an initialized Proteome object 
    
    """
    
    return fasta.fasta_to_proteome(filename, proteome=proteome, build_unique_ID=uniprot_accession_from_line, force_overwrite=force_overwrite, invalid_sequence_action=invalid_sequence_action)


## ------------------------------------------------------------------------
##
def uniprot_proteome_to_fasta(filename, proteome):                              
    """
    Stand alone function that allows the user to write a FASTA file from
    a Proteome under the assumption that the Proteome was built from a 
    uniprot FASTA.

    Practically, this just means that the Protein.name variable is used
    for the FASTA header, although the function will fail if duplicate
    headers are found.

    
    Parameters
    ------------

    filename : string
        Name of the FASTA file we're going to write sequences to

    proteome : Proteome
        The Proteome object from which FASTA file will be generated

    
    Returns 
    --------
    None
        No return variable but wll write to file 
    """

    out_dict = {}
    for p in proteome:

        header = p.name

        if header in out_dict:
            raise UtilitiesException(f'Duplicate name entries found in Proteome ({header}). Should not happen for UniProt headers')
        
        out_dict[header] = p.sequence

    protfasta.write_fasta(out_dict, filename, linelength=80)


## ###########################################################################
##
## UNIPROT REST API - LIVE ANNOTATION LOOKUP
##
## The functions below query the UniProt REST API and convert UniProt
## annotations into SHEPHARD Domains, Sites and Tracks. Note we deliberately
## use urllib from the standard library rather than requests, so that querying
## UniProt does not introduce a new runtime dependency for SHEPHARD.
##
## ###########################################################################

# base URL for the UniProtKB REST API
UNIPROT_REST_URL = 'https://rest.uniprot.org/uniprotkb'

# number of accessions requested in a single batch call. The /accessions
# endpoint tolerates more than this, but URLs get unwieldy and a smaller
# chunk means a transient failure costs us less work
MAX_ACCESSIONS_PER_REQUEST = 100

# prefix applied to every domain_type, site_type and track name generated here,
# so UniProt-derived annotations are always distinguishable from your own
DEFAULT_ANNOTATION_PREFIX = 'uniprot_'

# UniProt accession format, as defined by UniProt themselves. We allow an
# optional isoform suffix (e.g. P04637-2)
_ACCESSION_PATTERN = re.compile(r'^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[0-9]+)?$')

# HTTP status codes where retrying makes sense (rate limiting and transient
# server-side failures)
_RETRY_STATUS_CODES = (429, 500, 502, 503, 504)


##
## Definition of every annotation group that can be requested. Each group maps
## to:
##
##   fields : the UniProt REST 'fields' values needed to retrieve it
##   types  : the UniProt feature types the group selects
##   mode   : how the feature is converted into a SHEPHARD annotation
##
##             'region' - always added as a Domain
##             'auto'   - added as a Site if it covers a single residue,
##                        otherwise as a Domain
##             'paired' - the two positions are bonded partners rather than
##                        the ends of a region, so each is added as a Site
##                        (this is how UniProt encodes disulfide bonds and
##                        cross-links)
##
UNIPROT_ANNOTATION_GROUPS = {

    # ---- region-like features -------------------------------------------
    'domains':             {'fields': ('ft_domain',),   'types': ('Domain',),            'mode': 'region'},
    'regions':             {'fields': ('ft_region',),   'types': ('Region',),            'mode': 'region'},
    'motifs':              {'fields': ('ft_motif',),    'types': ('Motif',),             'mode': 'region'},
    'repeats':             {'fields': ('ft_repeat',),   'types': ('Repeat',),            'mode': 'region'},
    'zinc_fingers':        {'fields': ('ft_zn_fing',),  'types': ('Zinc finger',),       'mode': 'region'},
    'coiled_coils':        {'fields': ('ft_coiled',),   'types': ('Coiled coil',),       'mode': 'region'},
    'compositional_bias':  {'fields': ('ft_compbias',), 'types': ('Compositional bias',),'mode': 'region'},
    'dna_binding':         {'fields': ('ft_dna_bind',), 'types': ('DNA binding',),       'mode': 'region'},

    'transmembrane':       {'fields': ('ft_transmem', 'ft_intramem', 'ft_topo_dom'),
                            'types':  ('Transmembrane', 'Intramembrane', 'Topological domain'),
                            'mode':   'region'},

    'molecule_processing': {'fields': ('ft_chain', 'ft_peptide', 'ft_propep', 'ft_signal', 'ft_transit', 'ft_init_met'),
                            'types':  ('Chain', 'Peptide', 'Propeptide', 'Signal', 'Transit peptide', 'Initiator methionine'),
                            'mode':   'region'},

    'secondary_structure': {'fields': ('ft_helix', 'ft_strand', 'ft_turn'),
                            'types':  ('Helix', 'Beta strand', 'Turn'),
                            'mode':   'region'},

    # ---- site-like features ----------------------------------------------
    'binding_sites':       {'fields': ('ft_binding',),  'types': ('Binding site',),      'mode': 'auto'},
    'active_sites':        {'fields': ('ft_act_site',), 'types': ('Active site',),       'mode': 'auto'},
    'other_sites':         {'fields': ('ft_site',),     'types': ('Site',),              'mode': 'auto'},
    'modified_residues':   {'fields': ('ft_mod_res',),  'types': ('Modified residue',),  'mode': 'auto'},
    'glycosylation':       {'fields': ('ft_carbohyd',), 'types': ('Glycosylation',),     'mode': 'auto'},
    'lipidation':          {'fields': ('ft_lipid',),    'types': ('Lipidation',),        'mode': 'auto'},
    'mutagenesis':         {'fields': ('ft_mutagen',),  'types': ('Mutagenesis',),       'mode': 'auto'},
    'natural_variants':    {'fields': ('ft_variant',),  'types': ('Natural variant',),   'mode': 'auto'},

    # ---- bonded pairs ------------------------------------------------------
    'disulfide_bonds':     {'fields': ('ft_disulfid',), 'types': ('Disulfide bond',),    'mode': 'paired'},
    'cross_links':         {'fields': ('ft_crosslnk',), 'types': ('Cross-link',),        'mode': 'paired'},
}

# secondary structure symbols used by the secondary structure Track
_SECONDARY_STRUCTURE_SYMBOLS = {'Helix': 'H', 'Beta strand': 'E', 'Turn': 'T'}
_SECONDARY_STRUCTURE_BLANK = '-'


## ------------------------------------------------------------------------
##
def uniprot_annotation_groups():
    """
    Function that returns the names of every UniProt annotation group that
    can be requested by annotate_protein_with_uniprot() and
    annotate_proteome_with_uniprot().

    Returns
    -----------
    list of str
        Alphabetically sorted list of the valid annotation group names.

    """

    return sorted(UNIPROT_ANNOTATION_GROUPS.keys())


## ------------------------------------------------------------------------
##
def is_valid_uniprot_accession(accession):
    """
    Function that checks if a passed string looks like a valid UniProt
    accession. Note this checks the FORMAT only - it makes no statement as to
    whether the accession actually exists in UniProt.

    Parameters
    -----------
    accession : str
        String to be checked.

    Returns
    -----------
    bool
        Returns True if the passed string is a validly-formatted UniProt
        accession (with or without an isoform suffix), else False.

    """

    if not isinstance(accession, str):
        return False

    return _ACCESSION_PATTERN.match(accession.strip()) is not None


## ------------------------------------------------------------------------
##
def _http_get_json(url, timeout=30, n_retries=3, _sleep=time.sleep):
    """
    INTERNAL FUNCTION (not for public API use)

    Function that performs an HTTP GET and returns the decoded JSON response.
    This is the ONLY place in SHEPHARD where a network request is made, which
    means it is also the single seam that tests monkeypatch to run offline.

    Transient failures (rate limiting, 5xx responses, network errors) are
    retried with exponential backoff. Anything else raises immediately.

    Parameters
    ----------------
    url : str
        Fully-formed URL to request.

    timeout : int or float (default = 30)
        Per-request timeout in seconds.

    n_retries : int (default = 3)
        Number of times a transient failure is retried before giving up. Set
        to 0 to disable retrying.

    _sleep : function (default = time.sleep)
        Function used to wait between retries. Exposed so tests do not have
        to actually sleep.

    Returns
    ----------------
    dict
        Returns the decoded JSON response.

    Raises
    ----------------
    shephard.exceptions.APIException
        Raised if the request fails, the server returns an error status, or
        the response cannot be decoded as JSON.

    """

    if n_retries < 0:
        raise APIException(f'n_retries must be 0 or greater (got {n_retries})')

    # build a user agent that identifies SHEPHARD - UniProt ask for this so they
    # can contact us if we are misbehaving
    try:
        import shephard
        user_agent = f'SHEPHARD/{shephard.get_version()} (https://github.com/holehouse-lab/shephard)'
    except Exception:
        user_agent = 'SHEPHARD (https://github.com/holehouse-lab/shephard)'

    headers = {'Accept': 'application/json', 'User-Agent': user_agent}

    attempt = 0
    while True:

        try:
            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(request, timeout=timeout) as fh:
                raw = fh.read().decode('utf-8')

            try:
                return json.loads(raw)
            except ValueError as e:
                raise APIException(f'Could not decode the response from UniProt as JSON [{url}].\n\nError was: {e}')

        except urllib.error.HTTPError as e:

            # transient server-side problems and rate limiting are worth retrying
            if e.code in _RETRY_STATUS_CODES and attempt < n_retries:
                attempt = attempt + 1
                _sleep(_retry_delay(e, attempt))
                continue

            # try and surface UniProt's own error message, which is usually
            # much more useful than the status code alone
            detail = _extract_error_message(e)
            raise APIException(f'UniProt request failed with HTTP status {e.code} [{url}].\n\n{detail}')

        except urllib.error.URLError as e:

            if attempt < n_retries:
                attempt = attempt + 1
                _sleep(_retry_delay(None, attempt))
                continue

            raise APIException(f'Could not reach UniProt [{url}].\n\nError was: {e.reason}\n\nIs this machine online?')


## ------------------------------------------------------------------------
##
def _retry_delay(error, attempt):
    """
    INTERNAL FUNCTION (not for public API use)

    Function that computes how long to wait before retrying a failed request.
    If the server sent a Retry-After header we respect it, otherwise we back
    off exponentially (1s, 2s, 4s, ...).

    Parameters
    ----------------
    error : urllib.error.HTTPError or None
        The error that triggered the retry, if there was one.

    attempt : int
        The retry attempt we are about to make (indexed from 1).

    Returns
    ----------------
    float
        Number of seconds to wait.

    """

    if error is not None:
        try:
            retry_after = error.headers.get('Retry-After')
            if retry_after is not None:
                return float(retry_after)
        except (AttributeError, TypeError, ValueError):
            pass

    return float(2**(attempt - 1))


## ------------------------------------------------------------------------
##
def _extract_error_message(error):
    """
    INTERNAL FUNCTION (not for public API use)

    Function that pulls the human-readable error message out of a UniProt
    error response, falling back to the raw body if it is not JSON.

    Parameters
    ----------------
    error : urllib.error.HTTPError
        The error raised by urllib.

    Returns
    ----------------
    str
        Returns a message suitable for including in an exception.

    """

    try:
        body = error.read().decode('utf-8')
    except Exception:
        return 'No further information was returned by UniProt.'

    try:
        parsed = json.loads(body)
        messages = parsed.get('messages')
        if messages:
            return 'UniProt says: ' + '; '.join([str(m) for m in messages])
    except ValueError:
        pass

    return f'UniProt says: {body[:500]}'


## ------------------------------------------------------------------------
##
def _build_fields(annotation_groups, include_structures, verify_sequence):
    """
    INTERNAL FUNCTION (not for public API use)

    Function that builds the minimal set of UniProt REST 'fields' values
    needed to satisfy the requested annotation groups. Requesting only what
    we need keeps responses small - a full UniProtKB entry can be megabytes
    once natural variants are included.

    Parameters
    ----------------
    annotation_groups : list of str
        Names of the annotation groups that were requested.

    include_structures : bool
        Whether PDB cross-references are needed.

    verify_sequence : bool
        Whether the sequence is needed so it can be checked against the
        local Protein.

    Returns
    ----------------
    list of str
        Returns the ordered, de-duplicated list of field names.

    """

    # accession is always required so we can map responses back to proteins
    fields = ['accession']

    if verify_sequence:
        fields.append('sequence')

    for group in annotation_groups:
        for field in UNIPROT_ANNOTATION_GROUPS[group]['fields']:
            fields.append(field)

    if include_structures:
        fields.append('xref_pdb')

    # de-duplicate while preserving order
    return list(dict.fromkeys(fields))


## ------------------------------------------------------------------------
##
def _fetch_uniprot_records(accessions,
                           fields,
                           base_url=UNIPROT_REST_URL,
                           timeout=30,
                           n_retries=3,
                           chunk_size=MAX_ACCESSIONS_PER_REQUEST):
    """
    INTERNAL FUNCTION (not for public API use)

    Function that retrieves UniProt records for a list of accessions, batching
    the request into chunks so that a proteome-scale annotation is a handful
    of API calls rather than one call per protein.

    Parameters
    ----------------
    accessions : list of str
        Accessions to retrieve. Duplicates are collapsed.

    fields : list of str
        UniProt REST field names to retrieve.

    base_url : str (default = UNIPROT_REST_URL)
        Base URL of the UniProtKB REST API.

    timeout : int or float (default = 30)
        Per-request timeout in seconds.

    n_retries : int (default = 3)
        Number of retries on a transient failure.

    chunk_size : int (default = MAX_ACCESSIONS_PER_REQUEST)
        Number of accessions requested per API call.

    Returns
    ----------------
    dict
        Returns a dictionary that maps accession to the UniProt record. Note
        records are also indexed under any secondary accessions they carry,
        so a query using an obsolete accession still resolves.

    """

    if chunk_size < 1:
        raise APIException(f'chunk_size must be 1 or greater (got {chunk_size})')

    # collapse duplicates but keep the order stable so behaviour is predictable
    unique_accessions = list(dict.fromkeys(accessions))

    records = {}

    for idx in range(0, len(unique_accessions), chunk_size):

        chunk = unique_accessions[idx:idx + chunk_size]

        query = urllib.parse.urlencode({'accessions': ','.join(chunk),
                                        'fields': ','.join(fields)})

        response = _http_get_json(f'{base_url}/accessions?{query}',
                                  timeout=timeout,
                                  n_retries=n_retries)

        if not isinstance(response, dict):
            raise APIException(f'Unexpected response from UniProt - expected a JSON object but got {type(response)}')

        for record in response.get('results', []):

            primary = record.get('primaryAccession')
            if primary is None:
                continue

            records[primary] = record

            # index under secondary accessions too, but never let a secondary
            # accession displace a primary one
            for secondary in record.get('secondaryAccessions', []):
                if secondary not in records:
                    records[secondary] = record

    return records


## ------------------------------------------------------------------------
##
def _feature_positions(feature):
    """
    INTERNAL FUNCTION (not for public API use)

    Function that extracts the start and end positions from a UniProt feature.

    Note that UniProt positions can be unknown (a null value) or fuzzy (a
    modifier of 'OUTSIDE' or 'UNKNOWN'). We take the position at face value
    when it is present and return None when it is not, which lets the caller
    skip annotations that cannot be placed on a sequence.

    Parameters
    ----------------
    feature : dict
        A single UniProt feature dictionary.

    Returns
    ----------------
    tuple or None
        Returns a (start, end) tuple of ints, or None if the feature has no
        usable position.

    """

    location = feature.get('location')
    if not isinstance(location, dict):
        return None

    start = location.get('start', {})
    end = location.get('end', {})

    if not isinstance(start, dict) or not isinstance(end, dict):
        return None

    start_value = start.get('value')
    end_value = end.get('value')

    # a null value means UniProt does not know where this feature is
    if start_value is None or end_value is None:
        return None

    try:
        start_value = int(start_value)
        end_value = int(end_value)
    except (TypeError, ValueError):
        return None

    # UniProt should never do this, but if it did we would silently build a
    # backwards domain, so guard against it
    if end_value < start_value:
        return None

    return (start_value, end_value)


## ------------------------------------------------------------------------
##
def _format_evidence(feature):
    """
    INTERNAL FUNCTION (not for public API use)

    Function that flattens a UniProt feature's evidence list into a single
    string, so it can be stored as a SHEPHARD attribute and written to a
    SHEPHARD file without further processing.

    Evidence is formatted as ECO code, optionally qualified by its source,
    with multiple pieces of evidence separated by a comma - for example::

        ECO:0000269|PubMed:25732823,ECO:0000255

    Parameters
    ----------------
    feature : dict
        A single UniProt feature dictionary.

    Returns
    ----------------
    str
        Returns the formatted evidence string, or an empty string if the
        feature carries no evidence.

    """

    evidences = feature.get('evidences')
    if not evidences:
        return ''

    formatted = []
    for evidence in evidences:

        if not isinstance(evidence, dict):
            continue

        code = str(evidence.get('evidenceCode', '')).strip()
        source = str(evidence.get('source', '')).strip()
        identifier = str(evidence.get('id', '')).strip()

        if source and identifier:
            formatted.append(f'{code}|{source}:{identifier}')
        elif code:
            formatted.append(code)

    return ','.join(formatted)


## ------------------------------------------------------------------------
##
def _format_cross_references(feature):
    """
    INTERNAL FUNCTION (not for public API use)

    Function that flattens a UniProt feature's cross-references (for example
    the ChEBI identifier of a bound ligand) into a single string.

    Parameters
    ----------------
    feature : dict
        A single UniProt feature dictionary.

    Returns
    ----------------
    str
        Returns a string of the form 'ChEBI:CHEBI:29105', with multiple
        cross-references separated by a comma, or an empty string if there
        are none.

    """

    cross_references = feature.get('featureCrossReferences')
    if not cross_references:
        return ''

    formatted = []
    for reference in cross_references:
        if not isinstance(reference, dict):
            continue

        database = str(reference.get('database', '')).strip()
        identifier = str(reference.get('id', '')).strip()

        if database and identifier:
            formatted.append(f'{database}:{identifier}')

    return ','.join(formatted)


## ------------------------------------------------------------------------
##
def _feature_attributes(feature, accession, extra=None):
    """
    INTERNAL FUNCTION (not for public API use)

    Function that builds the SHEPHARD attribute dictionary attached to every
    annotation this module creates. Every value is a string so that
    annotations round-trip cleanly through SHEPHARD files.

    Parameters
    ----------------
    feature : dict
        A single UniProt feature dictionary.

    accession : str
        The UniProt accession the feature came from.

    extra : dict (default = None)
        Additional key-value pairs to merge into the attribute dictionary.

    Returns
    ----------------
    dict
        Returns the attribute dictionary.

    """

    attributes = {'source': 'uniprot',
                  'uniprot_accession': str(accession),
                  'uniprot_feature_type': str(feature.get('type', '')),
                  'description': str(feature.get('description', ''))}

    evidence = _format_evidence(feature)
    if evidence:
        attributes['evidence'] = evidence

    cross_references = _format_cross_references(feature)
    if cross_references:
        attributes['cross_references'] = cross_references

    feature_id = feature.get('featureId')
    if feature_id:
        attributes['uniprot_feature_id'] = str(feature_id)

    if extra:
        for k in extra:
            attributes[k] = str(extra[k])

    return attributes


## ------------------------------------------------------------------------
##
def _annotation_name(feature_type, prefix):
    """
    INTERNAL FUNCTION (not for public API use)

    Function that converts a UniProt feature type into the domain_type or
    site_type used in SHEPHARD - e.g. 'Beta strand' becomes
    'uniprot_beta_strand'.

    Parameters
    ----------------
    feature_type : str
        The UniProt feature type.

    prefix : str
        Prefix prepended to the generated name.

    Returns
    ----------------
    str
        Returns the SHEPHARD annotation name.

    """

    slug = str(feature_type).strip().lower()

    # collapse anything that is not alphanumeric into single underscores, so
    # names are always safe to write into a SHEPHARD file
    slug = re.sub(r'[^a-z0-9]+', '_', slug).strip('_')

    return f'{prefix}{slug}'


## ------------------------------------------------------------------------
##
def _parse_pdb_chain_ranges(chains_string):
    """
    INTERNAL FUNCTION (not for public API use)

    Function that parses the 'Chains' property of a UniProt PDB
    cross-reference into a list of (chains, start, end) tuples.

    The property looks like 'A/C=324-358' - i.e. one or more chain
    identifiers, then the range of the UniProt sequence those chains cover.
    Multiple ranges are separated by commas, e.g. 'A=1-100, B=200-300'.

    Parameters
    ----------------
    chains_string : str
        The raw 'Chains' property value.

    Returns
    ----------------
    list of tuple
        Returns a list of (chains, start, end) tuples. Entries that cannot be
        parsed are silently skipped, since a malformed chain string should
        not stop the rest of an entry being annotated.

    """

    parsed = []

    if not chains_string:
        return parsed

    for block in str(chains_string).split(','):

        block = block.strip()
        if '=' not in block:
            continue

        chains, _, position_range = block.partition('=')

        if '-' not in position_range:
            continue

        start_string, _, end_string = position_range.partition('-')

        try:
            start = int(start_string.strip())
            end = int(end_string.strip())
        except ValueError:
            continue

        if end < start:
            continue

        parsed.append((chains.strip(), start, end))

    return parsed


## ------------------------------------------------------------------------
##
def _pdb_structure_regions(record):
    """
    INTERNAL FUNCTION (not for public API use)

    Function that extracts the regions of a protein covered by an
    experimentally-determined structure, using the PDB cross-references in a
    UniProt record.

    Parameters
    ----------------
    record : dict
        A UniProt record.

    Returns
    ----------------
    list of dict
        Returns a list of dictionaries, each with 'start', 'end' and
        'attributes' keys, where the attributes carry the PDB identifier,
        experimental method, resolution and chains.

    """

    regions = []

    for cross_reference in record.get('uniProtKBCrossReferences', []):

        if not isinstance(cross_reference, dict):
            continue

        if cross_reference.get('database') != 'PDB':
            continue

        pdb_id = str(cross_reference.get('id', ''))

        # flatten the property list into something we can index into
        properties = {}
        for entry in cross_reference.get('properties', []):
            if isinstance(entry, dict) and 'key' in entry:
                properties[entry['key']] = entry.get('value', '')

        method = str(properties.get('Method', ''))
        resolution = str(properties.get('Resolution', ''))

        for (chains, start, end) in _parse_pdb_chain_ranges(properties.get('Chains', '')):

            regions.append({'start': start,
                            'end': end,
                            'attributes': {'source': 'uniprot',
                                           'uniprot_accession': str(record.get('primaryAccession', '')),
                                           'description': f'Experimental structure {pdb_id}',
                                           'pdb_id': pdb_id,
                                           'method': method,
                                           'resolution': resolution,
                                           'chains': chains}})

    return regions


## ------------------------------------------------------------------------
##
def _collect_features(record, annotation_groups):
    """
    INTERNAL FUNCTION (not for public API use)

    Function that walks the features in a UniProt record and returns only
    those belonging to the requested annotation groups, paired with the mode
    that says how each should be converted into a SHEPHARD annotation.

    Parameters
    ----------------
    record : dict
        A UniProt record.

    annotation_groups : list of str
        Names of the requested annotation groups.

    Returns
    ----------------
    list of tuple
        Returns a list of (feature, mode) tuples.

    """

    # build a lookup of feature type -> mode for the requested groups only
    type_to_mode = {}
    for group in annotation_groups:
        definition = UNIPROT_ANNOTATION_GROUPS[group]
        for feature_type in definition['types']:
            type_to_mode[feature_type] = definition['mode']

    collected = []
    for feature in record.get('features', []):

        if not isinstance(feature, dict):
            continue

        feature_type = feature.get('type')
        if feature_type in type_to_mode:
            collected.append((feature, type_to_mode[feature_type]))

    return collected


## ------------------------------------------------------------------------
##
def _annotate_protein_from_record(protein,
                                  record,
                                  annotation_groups,
                                  include_structures,
                                  secondary_structure_track,
                                  structure_coverage_track,
                                  prefix,
                                  safe,
                                  autoname,
                                  verbose):
    """
    INTERNAL FUNCTION (not for public API use)

    Function that applies a single UniProt record to a single Protein. This is
    the shared engine behind both annotate_protein_with_uniprot() and
    annotate_proteome_with_uniprot(), which is what guarantees that annotating
    one protein and annotating a whole proteome give identical results.

    Parameters
    ----------------
    protein : shephard.protein.Protein
        The Protein to annotate.

    record : dict
        The UniProt record for this Protein.

    annotation_groups : list of str
        Names of the requested annotation groups.

    include_structures : bool
        Whether to add Domains for regions with an experimental structure.

    secondary_structure_track : bool
        Whether to add a symbols Track of secondary structure.

    structure_coverage_track : bool
        Whether to add a values Track of experimental structure coverage.

    prefix : str
        Prefix applied to every generated domain_type, site_type and Track
        name.

    safe : bool
        If True, problems (out-of-range annotations, clashes with existing
        annotations) raise an exception. If False they are skipped, with a
        warning printed when verbose is True.

    autoname : bool
        Passed to Protein.add_domain(). Should normally be left True, since
        UniProt regularly reports several distinct features that share their
        coordinates.

    verbose : bool
        Whether to print warnings about skipped annotations.

    Returns
    ----------------
    dict
        Returns a summary dictionary with the counts of added 'domains',
        'sites' and 'tracks', and of 'skipped' annotations.

    """

    summary = {'domains': 0, 'sites': 0, 'tracks': 0, 'skipped': 0}

    accession = record.get('primaryAccession', protein.unique_ID)
    protein_length = len(protein)

    def _in_range(position):
        return 1 <= position <= protein_length

    def _handle_problem(message):
        """
        Deals with an annotation we cannot add - either by raising or by
        warning and moving on.
        """
        summary['skipped'] = summary['skipped'] + 1
        if safe:
            raise APIException(message)
        if verbose:
            shephard_exceptions.print_warning(message)

    ## ---- features -------------------------------------------------------
    for (feature, mode) in _collect_features(record, annotation_groups):

        feature_type = feature.get('type')
        positions = _feature_positions(feature)

        if positions is None:
            _handle_problem(f'Skipping a [{feature_type}] feature on {accession} because UniProt does not give it a usable position')
            continue

        (start, end) = positions

        if not _in_range(start) or not _in_range(end):
            _handle_problem(f'Skipping a [{feature_type}] feature on {accession} at {start}-{end} because it falls outside the protein (length {protein_length}). This usually means the local sequence and the UniProt sequence differ')
            continue

        try:
            # disulfide bonds and cross-links record two bonded residues rather
            # than the ends of a region, so each partner becomes its own Site
            if mode == 'paired':

                site_type = _annotation_name(feature_type, prefix)

                if start == end:
                    protein.add_site(start, site_type, attributes=_feature_attributes(feature, accession))
                    summary['sites'] = summary['sites'] + 1
                else:
                    protein.add_site(start, site_type, attributes=_feature_attributes(feature, accession, extra={'partner_position': end}))
                    protein.add_site(end, site_type, attributes=_feature_attributes(feature, accession, extra={'partner_position': start}))
                    summary['sites'] = summary['sites'] + 2

            # a single-residue feature becomes a Site, anything longer a Domain
            elif mode == 'auto' and start == end:

                protein.add_site(start, _annotation_name(feature_type, prefix), attributes=_feature_attributes(feature, accession))
                summary['sites'] = summary['sites'] + 1

            else:

                protein.add_domain(start, end, _annotation_name(feature_type, prefix),
                                   attributes=_feature_attributes(feature, accession),
                                   safe=safe,
                                   autoname=autoname)
                summary['domains'] = summary['domains'] + 1

        except ShephardException as e:
            _handle_problem(f'Could not add a [{feature_type}] annotation on {accession} at {start}-{end}: {e}')

    ## ---- regions covered by an experimental structure --------------------
    structure_regions = []
    if include_structures or structure_coverage_track:
        structure_regions = _pdb_structure_regions(record)

    if include_structures:

        domain_type = f'{prefix}experimental_structure'

        for region in structure_regions:

            start = region['start']
            end = region['end']

            if not _in_range(start) or not _in_range(end):
                _handle_problem(f"Skipping experimental structure {region['attributes']['pdb_id']} on {accession} at {start}-{end} because it falls outside the protein (length {protein_length})")
                continue

            try:
                protein.add_domain(start, end, domain_type,
                                   attributes=dict(region['attributes']),
                                   safe=safe,
                                   autoname=autoname)
                summary['domains'] = summary['domains'] + 1

            except ShephardException as e:
                _handle_problem(f"Could not add experimental structure {region['attributes']['pdb_id']} on {accession}: {e}")

    ## ---- secondary structure track ---------------------------------------
    if secondary_structure_track:

        symbols = [_SECONDARY_STRUCTURE_BLANK]*protein_length

        for feature in record.get('features', []):

            if not isinstance(feature, dict):
                continue

            symbol = _SECONDARY_STRUCTURE_SYMBOLS.get(feature.get('type'))
            if symbol is None:
                continue

            positions = _feature_positions(feature)
            if positions is None:
                continue

            (start, end) = positions
            if not _in_range(start) or not _in_range(end):
                continue

            for position in range(start, end + 1):
                symbols[position - 1] = symbol

        try:
            protein.add_track(f'{prefix}secondary_structure', symbols=symbols, safe=safe)
            summary['tracks'] = summary['tracks'] + 1
        except ShephardException as e:
            _handle_problem(f'Could not add the secondary structure track on {accession}: {e}')

    ## ---- experimental structure coverage track ---------------------------
    if structure_coverage_track:

        values = [0.0]*protein_length

        for region in structure_regions:

            start = max(1, region['start'])
            end = min(protein_length, region['end'])

            for position in range(start, end + 1):
                values[position - 1] = 1.0

        try:
            protein.add_track(f'{prefix}structure_coverage', values=values, safe=safe)
            summary['tracks'] = summary['tracks'] + 1
        except ShephardException as e:
            _handle_problem(f'Could not add the structure coverage track on {accession}: {e}')

    return summary


## ------------------------------------------------------------------------
##
def _resolve_annotation_groups(requested):
    """
    INTERNAL FUNCTION (not for public API use)

    Function that converts the boolean keyword arguments of the public
    functions into a validated list of annotation group names.

    Parameters
    ----------------
    requested : dict
        Dictionary that maps annotation group name to a bool.

    Returns
    ----------------
    list of str
        Returns the names of the requested groups, in the order they are
        defined in UNIPROT_ANNOTATION_GROUPS.

    """

    for name in requested:
        if name not in UNIPROT_ANNOTATION_GROUPS:
            raise APIException(f'Unknown UniProt annotation group [{name}]. Valid options are: {uniprot_annotation_groups()!s}')

    return [name for name in UNIPROT_ANNOTATION_GROUPS if requested.get(name, False)]


## ------------------------------------------------------------------------
##
def _check_sequence_matches(protein, record, safe, verbose):
    """
    INTERNAL FUNCTION (not for public API use)

    Function that checks the local Protein sequence against the sequence
    UniProt holds. This matters because every UniProt annotation is a set of
    sequence positions - if the two sequences differ (a different isoform, an
    older release, a sequence that has been edited locally) then those
    positions describe a different protein and the resulting annotations
    would be quietly wrong.

    Parameters
    ----------------
    protein : shephard.protein.Protein
        The Protein being annotated.

    record : dict
        The UniProt record.

    safe : bool
        If True a mismatch raises, if False it warns and reports failure.

    verbose : bool
        Whether to print a warning when safe is False.

    Returns
    ----------------
    bool
        Returns True if the sequences match (or if UniProt did not return a
        sequence to compare against), else False.

    """

    sequence_block = record.get('sequence')

    # nothing to compare against - this happens if the caller turned off
    # sequence verification and we therefore never requested the sequence
    if not isinstance(sequence_block, dict):
        return True

    uniprot_sequence = sequence_block.get('value')
    if uniprot_sequence is None:
        return True

    if protein.sequence == uniprot_sequence:
        return True

    accession = record.get('primaryAccession', protein.unique_ID)
    message = (f'The sequence for {protein.unique_ID} does not match the UniProt sequence for {accession} '
               f'(local length {len(protein)}, UniProt length {len(uniprot_sequence)}). UniProt annotations are '
               f'positional, so applying them to a different sequence would give incorrect annotations. Set '
               f'verify_sequence=False to annotate anyway.')

    if safe:
        raise APIException(message)

    if verbose:
        shephard_exceptions.print_warning(message)

    return False


## ------------------------------------------------------------------------
##
def annotate_proteome_with_uniprot(proteome,
                                   domains=True,
                                   regions=True,
                                   motifs=True,
                                   repeats=True,
                                   zinc_fingers=True,
                                   coiled_coils=True,
                                   compositional_bias=True,
                                   dna_binding=True,
                                   transmembrane=True,
                                   molecule_processing=False,
                                   secondary_structure=False,
                                   binding_sites=True,
                                   active_sites=True,
                                   other_sites=True,
                                   modified_residues=True,
                                   glycosylation=True,
                                   lipidation=True,
                                   mutagenesis=False,
                                   natural_variants=False,
                                   disulfide_bonds=True,
                                   cross_links=True,
                                   experimental_structures=False,
                                   secondary_structure_track=False,
                                   structure_coverage_track=False,
                                   accession_from_protein=None,
                                   verify_sequence=True,
                                   prefix=DEFAULT_ANNOTATION_PREFIX,
                                   safe=True,
                                   autoname=True,
                                   verbose=True,
                                   base_url=UNIPROT_REST_URL,
                                   timeout=30,
                                   n_retries=3,
                                   chunk_size=MAX_ACCESSIONS_PER_REQUEST):
    """
    Function that annotates every Protein in a Proteome with annotations
    pulled live from the UniProt REST API.

    This is the batch equivalent of annotate_protein_with_uniprot(). Rather
    than making one API call per protein, accessions are requested in batches
    (of chunk_size), so annotating a proteome of a few thousand proteins is a
    few tens of API calls. The records that come back are then processed
    exactly as they are for a single protein, so the two functions always
    give the same answer.

    By default the unique_ID of each Protein is used as its UniProt
    accession, which is why we recommend using UniProt accessions as
    unique_IDs. If your unique_IDs are something else, pass a function via
    accession_from_protein.

    Which annotations are retrieved is controlled by the boolean flags below.
    Only the data needed for the requested annotations is downloaded, so
    turning flags off makes the query faster and the response smaller.

    UniProt features are converted to SHEPHARD annotations as follows:

    * Features that span a region become **Domains**

    * Features on a single residue become **Sites**

    * Disulfide bonds and cross-links become a pair of **Sites**, one on each
      bonded residue, each recording its partner's position

    * Optionally, secondary structure and experimental structure coverage
      become **Tracks**

    Every annotation carries an attribute dictionary with the UniProt feature
    type, its description, its evidence codes, any cross-references (e.g. the
    ChEBI ID of a bound ligand), the UniProt feature ID where one exists, and
    the accession the annotation came from.

    Parameters
    -----------------
    proteome : shephard.proteome.Proteome
        Proteome object to be annotated.

    domains : bool (default = True)
        Annotate UniProt 'Domain' features.

    regions : bool (default = True)
        Annotate 'Region' features - the free-text regions of interest that
        cover things like interaction interfaces.

    motifs : bool (default = True)
        Annotate 'Motif' features.

    repeats : bool (default = True)
        Annotate 'Repeat' features.

    zinc_fingers : bool (default = True)
        Annotate 'Zinc finger' features.

    coiled_coils : bool (default = True)
        Annotate 'Coiled coil' features.

    compositional_bias : bool (default = True)
        Annotate 'Compositional bias' features.

    dna_binding : bool (default = True)
        Annotate 'DNA binding' features.

    transmembrane : bool (default = True)
        Annotate 'Transmembrane', 'Intramembrane' and 'Topological domain'
        features.

    molecule_processing : bool (default = False)
        Annotate 'Chain', 'Peptide', 'Propeptide', 'Signal', 'Transit
        peptide' and 'Initiator methionine' features. Off by default because
        the 'Chain' feature usually spans the entire protein.

    secondary_structure : bool (default = False)
        Annotate 'Helix', 'Beta strand' and 'Turn' features as Domains. Note
        that for most purposes secondary_structure_track is more useful.

    binding_sites : bool (default = True)
        Annotate 'Binding site' features.

    active_sites : bool (default = True)
        Annotate 'Active site' features.

    other_sites : bool (default = True)
        Annotate general 'Site' features.

    modified_residues : bool (default = True)
        Annotate 'Modified residue' features (phosphorylation, acetylation
        and so on).

    glycosylation : bool (default = True)
        Annotate 'Glycosylation' features.

    lipidation : bool (default = True)
        Annotate 'Lipidation' features.

    mutagenesis : bool (default = False)
        Annotate 'Mutagenesis' features. Off by default as these describe
        experiments rather than the native protein.

    natural_variants : bool (default = False)
        Annotate 'Natural variant' features. Off by default because
        well-studied proteins can carry thousands of these.

    disulfide_bonds : bool (default = True)
        Annotate 'Disulfide bond' features as paired Sites.

    cross_links : bool (default = True)
        Annotate 'Cross-link' features as paired Sites.

    experimental_structures : bool (default = False)
        Add a Domain for every region of the protein covered by an
        experimentally-determined structure, taken from the PDB
        cross-references. Each Domain records the PDB identifier, the
        experimental method, the resolution and the chains.

    secondary_structure_track : bool (default = False)
        Add a symbols Track of the UniProt secondary structure, where each
        residue is 'H' (helix), 'E' (beta strand), 'T' (turn) or '-'.

    structure_coverage_track : bool (default = False)
        Add a values Track where each residue is 1 if it is covered by an
        experimental structure and 0 if it is not.

    accession_from_protein : function (default = None)
        Optional function that takes a Protein and returns the UniProt
        accession to look it up by. If None the Protein's unique_ID is used.

    verify_sequence : bool (default = True)
        If True, the local sequence is checked against the UniProt sequence
        and proteins whose sequences differ are not annotated. Strongly
        recommended - UniProt annotations are positional, so applying them to
        a different sequence gives wrong answers rather than no answer.

    prefix : str (default = 'uniprot_')
        Prefix applied to every generated domain_type, site_type and Track
        name.

    safe : bool (default = True)
        If True, any problem (an accession that is missing from UniProt, a
        sequence mismatch, an annotation that falls outside the protein, a
        clash with an existing annotation) raises an exception. If False
        these are skipped and reported in the returned summary.

    autoname : bool (default = True)
        Passed through to Protein.add_domain(). Recommended to leave True,
        since UniProt regularly reports several distinct features that share
        the same coordinates.

    verbose : bool (default = True)
        Whether to print warnings about skipped annotations.

    base_url : str (default = UNIPROT_REST_URL)
        Base URL of the UniProtKB REST API.

    timeout : int or float (default = 30)
        Per-request timeout in seconds.

    n_retries : int (default = 3)
        Number of times a transient failure (rate limiting, a 5xx response, a
        network blip) is retried before giving up.

    chunk_size : int (default = 100)
        Number of accessions requested per API call.

    Returns
    -----------------
    dict
        Returns a summary dictionary keyed by unique_ID, where each value is
        a dictionary with the counts of 'domains', 'sites' and 'tracks' added
        and of 'skipped' annotations. Proteins that could not be annotated at
        all appear with an 'error' key explaining why.

    Raises
    -----------------
    shephard.exceptions.APIException
        Raised if the UniProt request fails, or - when safe is True - if any
        protein cannot be annotated.

    """

    from shephard.interfaces import interface_tools
    interface_tools.check_proteome(proteome, 'annotate_proteome_with_uniprot (uniprot)')

    annotation_groups = _resolve_annotation_groups(
        {'domains': domains, 'regions': regions, 'motifs': motifs,
         'repeats': repeats, 'zinc_fingers': zinc_fingers,
         'coiled_coils': coiled_coils, 'compositional_bias': compositional_bias,
         'dna_binding': dna_binding, 'transmembrane': transmembrane,
         'molecule_processing': molecule_processing,
         'secondary_structure': secondary_structure,
         'binding_sites': binding_sites, 'active_sites': active_sites,
         'other_sites': other_sites, 'modified_residues': modified_residues,
         'glycosylation': glycosylation, 'lipidation': lipidation,
         'mutagenesis': mutagenesis, 'natural_variants': natural_variants,
         'disulfide_bonds': disulfide_bonds, 'cross_links': cross_links})

    # the secondary structure track needs the same features as the secondary
    # structure group, so make sure they are requested even if the Domains
    # themselves were not asked for
    if secondary_structure_track and 'secondary_structure' not in annotation_groups:
        fields_groups = annotation_groups + ['secondary_structure']
    else:
        fields_groups = annotation_groups

    need_structures = experimental_structures or structure_coverage_track

    if len(annotation_groups) == 0 and not need_structures and not secondary_structure_track:
        raise APIException(f'No annotations were requested, so there is nothing to do. Valid annotation groups are: {uniprot_annotation_groups()!s}')

    fields = _build_fields(fields_groups, need_structures, verify_sequence)

    # work out which accession to look each protein up by
    summary = {}
    accession_to_proteins = {}

    for protein in proteome:

        if accession_from_protein is None:
            accession = protein.unique_ID
        else:
            accession = accession_from_protein(protein)

        if not is_valid_uniprot_accession(accession):
            message = f'[{accession}] (from protein {protein.unique_ID}) is not a validly-formatted UniProt accession'
            if safe:
                raise APIException(message)
            if verbose:
                shephard_exceptions.print_warning(message)
            summary[protein.unique_ID] = {'domains': 0, 'sites': 0, 'tracks': 0,
                                          'skipped': 0, 'error': 'invalid accession'}
            continue

        accession = accession.strip()
        accession_to_proteins.setdefault(accession, []).append(protein)

    # nothing left to look up
    if len(accession_to_proteins) == 0:
        return summary

    records = _fetch_uniprot_records(list(accession_to_proteins.keys()),
                                     fields,
                                     base_url=base_url,
                                     timeout=timeout,
                                     n_retries=n_retries,
                                     chunk_size=chunk_size)

    for accession in accession_to_proteins:

        record = records.get(accession)

        for protein in accession_to_proteins[accession]:

            if record is None:
                message = f'UniProt returned no record for accession [{accession}] (protein {protein.unique_ID})'
                if safe:
                    raise APIException(message)
                if verbose:
                    shephard_exceptions.print_warning(message)
                summary[protein.unique_ID] = {'domains': 0, 'sites': 0, 'tracks': 0,
                                              'skipped': 0, 'error': 'no record returned'}
                continue

            if verify_sequence and not _check_sequence_matches(protein, record, safe, verbose):
                summary[protein.unique_ID] = {'domains': 0, 'sites': 0, 'tracks': 0,
                                              'skipped': 0, 'error': 'sequence mismatch'}
                continue

            summary[protein.unique_ID] = _annotate_protein_from_record(
                protein,
                record,
                annotation_groups,
                include_structures=experimental_structures,
                secondary_structure_track=secondary_structure_track,
                structure_coverage_track=structure_coverage_track,
                prefix=prefix,
                safe=safe,
                autoname=autoname,
                verbose=verbose)

    return summary


## ------------------------------------------------------------------------
##
def annotate_protein_with_uniprot(protein,
                                  domains=True,
                                  regions=True,
                                  motifs=True,
                                  repeats=True,
                                  zinc_fingers=True,
                                  coiled_coils=True,
                                  compositional_bias=True,
                                  dna_binding=True,
                                  transmembrane=True,
                                  molecule_processing=False,
                                  secondary_structure=False,
                                  binding_sites=True,
                                  active_sites=True,
                                  other_sites=True,
                                  modified_residues=True,
                                  glycosylation=True,
                                  lipidation=True,
                                  mutagenesis=False,
                                  natural_variants=False,
                                  disulfide_bonds=True,
                                  cross_links=True,
                                  experimental_structures=False,
                                  secondary_structure_track=False,
                                  structure_coverage_track=False,
                                  accession=None,
                                  verify_sequence=True,
                                  prefix=DEFAULT_ANNOTATION_PREFIX,
                                  safe=True,
                                  autoname=True,
                                  verbose=True,
                                  base_url=UNIPROT_REST_URL,
                                  timeout=30,
                                  n_retries=3):
    """
    Function that annotates a single Protein with annotations pulled live
    from the UniProt REST API.

    The function is stateless - it makes one API call, converts the returned
    UniProt features into SHEPHARD annotations, and applies them to the
    passed Protein.

    By default the Protein's unique_ID is used as the UniProt accession,
    which is why we recommend using UniProt accessions as unique_IDs. If your
    unique_IDs are something else, pass the accession explicitly.

    If you are annotating more than a handful of proteins, use
    annotate_proteome_with_uniprot() instead - it batches accessions into a
    small number of API calls rather than making one call per protein, and
    processes the results identically.

    UniProt features are converted to SHEPHARD annotations as follows:

    * Features that span a region become **Domains**

    * Features on a single residue become **Sites**

    * Disulfide bonds and cross-links become a pair of **Sites**, one on each
      bonded residue, each recording its partner's position

    * Optionally, secondary structure and experimental structure coverage
      become **Tracks**

    Every annotation carries an attribute dictionary with the UniProt feature
    type, its description, its evidence codes, any cross-references (e.g. the
    ChEBI ID of a bound ligand), the UniProt feature ID where one exists, and
    the accession the annotation came from.

    Parameters
    -----------------
    protein : shephard.protein.Protein
        Protein object to be annotated.

    domains : bool (default = True)
        Annotate UniProt 'Domain' features.

    regions : bool (default = True)
        Annotate 'Region' features.

    motifs : bool (default = True)
        Annotate 'Motif' features.

    repeats : bool (default = True)
        Annotate 'Repeat' features.

    zinc_fingers : bool (default = True)
        Annotate 'Zinc finger' features.

    coiled_coils : bool (default = True)
        Annotate 'Coiled coil' features.

    compositional_bias : bool (default = True)
        Annotate 'Compositional bias' features.

    dna_binding : bool (default = True)
        Annotate 'DNA binding' features.

    transmembrane : bool (default = True)
        Annotate 'Transmembrane', 'Intramembrane' and 'Topological domain'
        features.

    molecule_processing : bool (default = False)
        Annotate 'Chain', 'Peptide', 'Propeptide', 'Signal', 'Transit
        peptide' and 'Initiator methionine' features.

    secondary_structure : bool (default = False)
        Annotate 'Helix', 'Beta strand' and 'Turn' features as Domains.

    binding_sites : bool (default = True)
        Annotate 'Binding site' features.

    active_sites : bool (default = True)
        Annotate 'Active site' features.

    other_sites : bool (default = True)
        Annotate general 'Site' features.

    modified_residues : bool (default = True)
        Annotate 'Modified residue' features.

    glycosylation : bool (default = True)
        Annotate 'Glycosylation' features.

    lipidation : bool (default = True)
        Annotate 'Lipidation' features.

    mutagenesis : bool (default = False)
        Annotate 'Mutagenesis' features.

    natural_variants : bool (default = False)
        Annotate 'Natural variant' features. Off by default because
        well-studied proteins can carry thousands of these.

    disulfide_bonds : bool (default = True)
        Annotate 'Disulfide bond' features as paired Sites.

    cross_links : bool (default = True)
        Annotate 'Cross-link' features as paired Sites.

    experimental_structures : bool (default = False)
        Add a Domain for every region covered by an experimentally-determined
        structure, taken from the PDB cross-references, recording the PDB
        identifier, experimental method, resolution and chains.

    secondary_structure_track : bool (default = False)
        Add a symbols Track of secondary structure ('H', 'E', 'T' or '-').

    structure_coverage_track : bool (default = False)
        Add a values Track that is 1 for residues covered by an experimental
        structure and 0 elsewhere.

    accession : str (default = None)
        UniProt accession to look this Protein up by. If None the Protein's
        unique_ID is used.

    verify_sequence : bool (default = True)
        If True, the local sequence is checked against the UniProt sequence
        and the protein is not annotated if they differ. Strongly
        recommended - UniProt annotations are positional.

    prefix : str (default = 'uniprot_')
        Prefix applied to every generated domain_type, site_type and Track
        name.

    safe : bool (default = True)
        If True, any problem raises an exception. If False, problems are
        skipped and counted in the returned summary.

    autoname : bool (default = True)
        Passed through to Protein.add_domain().

    verbose : bool (default = True)
        Whether to print warnings about skipped annotations.

    base_url : str (default = UNIPROT_REST_URL)
        Base URL of the UniProtKB REST API.

    timeout : int or float (default = 30)
        Per-request timeout in seconds.

    n_retries : int (default = 3)
        Number of times a transient failure is retried before giving up.

    Returns
    -----------------
    dict
        Returns a summary dictionary with the counts of 'domains', 'sites'
        and 'tracks' added, and of 'skipped' annotations. If the protein
        could not be annotated the dictionary also carries an 'error' key.

    Raises
    -----------------
    shephard.exceptions.APIException
        Raised if the UniProt request fails, or - when safe is True - if the
        protein cannot be annotated.

    """

    from shephard.interfaces import interface_tools
    interface_tools.check_protein(protein, 'annotate_protein_with_uniprot (uniprot)')

    if accession is None:
        accession = protein.unique_ID

    if not is_valid_uniprot_accession(accession):
        message = f'[{accession}] is not a validly-formatted UniProt accession'
        if safe:
            raise APIException(message)
        if verbose:
            shephard_exceptions.print_warning(message)
        return {'domains': 0, 'sites': 0, 'tracks': 0, 'skipped': 0,
                'error': 'invalid accession'}

    accession = accession.strip()

    annotation_groups = _resolve_annotation_groups(
        {'domains': domains, 'regions': regions, 'motifs': motifs,
         'repeats': repeats, 'zinc_fingers': zinc_fingers,
         'coiled_coils': coiled_coils, 'compositional_bias': compositional_bias,
         'dna_binding': dna_binding, 'transmembrane': transmembrane,
         'molecule_processing': molecule_processing,
         'secondary_structure': secondary_structure,
         'binding_sites': binding_sites, 'active_sites': active_sites,
         'other_sites': other_sites, 'modified_residues': modified_residues,
         'glycosylation': glycosylation, 'lipidation': lipidation,
         'mutagenesis': mutagenesis, 'natural_variants': natural_variants,
         'disulfide_bonds': disulfide_bonds, 'cross_links': cross_links})

    if secondary_structure_track and 'secondary_structure' not in annotation_groups:
        fields_groups = annotation_groups + ['secondary_structure']
    else:
        fields_groups = annotation_groups

    need_structures = experimental_structures or structure_coverage_track

    if len(annotation_groups) == 0 and not need_structures and not secondary_structure_track:
        raise APIException(f'No annotations were requested, so there is nothing to do. Valid annotation groups are: {uniprot_annotation_groups()!s}')

    fields = _build_fields(fields_groups, need_structures, verify_sequence)

    records = _fetch_uniprot_records([accession],
                                     fields,
                                     base_url=base_url,
                                     timeout=timeout,
                                     n_retries=n_retries,
                                     chunk_size=1)

    record = records.get(accession)

    if record is None:
        message = f'UniProt returned no record for accession [{accession}]'
        if safe:
            raise APIException(message)
        if verbose:
            shephard_exceptions.print_warning(message)
        return {'domains': 0, 'sites': 0, 'tracks': 0, 'skipped': 0,
                'error': 'no record returned'}

    if verify_sequence and not _check_sequence_matches(protein, record, safe, verbose):
        return {'domains': 0, 'sites': 0, 'tracks': 0, 'skipped': 0,
                'error': 'sequence mismatch'}

    return _annotate_protein_from_record(protein,
                                         record,
                                         annotation_groups,
                                         include_structures=experimental_structures,
                                         secondary_structure_track=secondary_structure_track,
                                         structure_coverage_track=structure_coverage_track,
                                         prefix=prefix,
                                         safe=safe,
                                         autoname=autoname,
                                         verbose=verbose)
