from localcider.sequenceParameters import SequenceParameters
from shephard.exceptions import InterfaceException
from shephard.interfaces import interface_tools 
from localcider import SequenceException
from shephard import general_utilities

"""
Functions associated with interfacing with the localcider package 

"""


def apply_track_NCPR(proteome, blobsize=5):    
    """
    
    """

    interface_tools.check_proteome(proteome, 'apply_track_NCPR (si_localcider)')

    for protein in proteome:

        seq = protein.sequence

        # skip short sequences
        if len(seq) < blobsize:
            continue

        ncpr=[]
        for i in seq:

            if i in ['D','E']:
                ncpr.append(-1)
            elif i in ['K','R']:
                ncpr.append(+1)
            else:
                ncpr.append(0)

        averaged = []
        for i in range(0,(len(ncpr)-blobsize)):
            averaged.append(general_utilities.numerical_average(ncpr[i:i+blobsize]))
            

        ntd_mod = int(blobsize/2)
        ctd_mod = blobsize - ntd_mod
        final_averaged = [averaged[0]]*ntd_mod + averaged + [averaged[-1]]*ctd_mod
        

        """
        # this is where we convert sequence into NCPR
        try:
            ncpr = SequenceParameters(seq).get_linear_NCPR()[1]
        except SequenceException:
            seq = seq.replace('U','C')
            ncpr = SequenceParameters(seq).get_linear_NCPR()[1]
        """

        # this is where we add the NCPR track
        protein.add_track('NCPR', values=final_averaged)


def apply_track_residue_density(proteome, residue_set, name='residue_density', block_size=30, safe=True):    
    """
    
    This is an apply_track_* function that generates a values-track that describes the local density of amino 
    acids along the sequence. The set of amino acids to be used is defined by the residue_set, 
    while name and block_size are optional parameters that define different things.
    
    Density is calculated by computing the number of residues within a block_size subsequence and dividing by
    the block_size. This yeilds a value between 0 and 1.

    NB: If a protein sequence is shorter than the block_size then that protein gets skipped. 

    Parameters
    ------------------

    proteome : proteome object
        This is the proteome object, and for each protein the residue_density function will be applied. Note
        that if this is not a proteome object this function will throw an exception.

    residue_set : list of single-character strings representing amino acids
        This is a list of amino acids for which the local density will be computed. This could be a single 
        residue, or multiple residues. Duplicates are removed. Note that beyond removing duplictes and ensuring
        residues are upper case no sanity check is done here.

    name : string 
        The name defines the track name and allows the user to customize what they want to call this set of
        tracks. If safe is set to true (default) and a name is proposed that already and safe=True then 
        the function will throw an exception. If safe=False then an existing track will be overwritten.

    block_size : int (default = 30)
        This is the size of the subsequence over which density is calculated. 

    safe : bool (default = True)
        Boolean that defines if track names should be overwritten (if False) or throw an error (if True)
        when a track name is proposed that already exists in the proteome.


    Return 
    ----------

    No return value, but the proteome will have the set of tracks added to all proteins of sufficient
    length

    """

    # validate that a proteome was actually passed
    interface_tools.check_proteome(proteome, 'apply_track_residue_density (si_localcider)')
    
    # first remove any duplicates
    try:        
        residue_set = list(set(residue_set))
    except Exception:
        raise InterfaceException('Error when applying "apply_track_residue_density" from the localcider interface. Could not convert the residue_set to a set. residue_set should be a list/tuple/set of amino acid residues which will be used to compute local density')

    # now for each protein
    for protein in proteome:

        # get the protein sequence and number of residues. If there are fewer residues in
        # the protein than the block size then we skip this protein.
        seq = protein.sequence
        nres = len(seq)
        if nres < block_size:
            continue

        # this is where we convert sequence into a local density by cylcing through
        # each $block_size subsequence, computing the local density, and appending 
        # that density to the ever-growing density vector
        density_vector = []
        for pos in range(0, (nres - block_size)+1):

            total = 0            
            subseq = seq[pos:pos+block_size]            
            for i in residue_set:
                total = total + subseq.count(i.upper())

            density_vector.append(total/block_size)

        # this code creates leading/lagging values to fill in the missing ones such that the actual value
        # of the density vector reports on the density half-way across the blocksize and the 
        # track length = nres
        leading_values = [density_vector[0]]*int(block_size/2)
        lagging_values = [density_vector[-1]]*(nres - (len(density_vector) + len(leading_values)))
          
        # this line then combines the leading, desnity and lagging lists into a single list, and we
        # then add this numerical list as a track
        final = leading_values + density_vector + lagging_values    
        protein.add_track(name, values=final, safe=safe)


def apply_attribute_kappa(proteome):    

    interface_tools.check_proteome(proteome, 'apply_attribute_kappa (si_localcider)')

    for protein in proteome:

        # get the protein sequence
        seq = protein.sequence

        # this is where we convert sequence into NCPR
        kappa = SequenceParameters(protein.sequence).get_kappa()

        # this is where we add the NCPR track
        protein.add_attribute('kappa', kappa)


