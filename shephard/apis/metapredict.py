

##
## API into metapredict
##
## For all APIS we do do not make hard dependencies, safe informative
## import checking should be done.
##


try:
    import metapredict as meta
    
except ModuleNotFoundError:
    print('Unable to import metapredict')
    print('To use the metapredict API, make sure metapredict is installed')
    print('This can be done as follows:')
    print('pip install metapredict')

def annotate_proteome_with_disorder_track(proteome,                                        
                                          name='disorder',
                                          gpuid=00,
                                          show_progress_bar=True,
                                          batch_mode=None,
                                          safe=True):

    """
    Function that annotates a proteome with disorder Tracks 
    for every protein. 

    By default, disorder Tracks are named 'disorder', although
    this can be changed by setting the `track_name` parameter.

    Disorder prediction uses the batch mode in metapredict, which 
    leverages parallel predictions automatically on GPUs or CPUs.
    However, if a specific device is requested, this can be passed

    Parameters
    -----------------
    proteome : shephard.proteome.Proteome 
        Proteome object to be annotated.

    track_name : str
        Name of the Track added to each Protein. 
        Default = 'disorder'

    gpuid : int 
        Identifier for the GPU being requested. Note that if
        this is left unset the code will use the first GPU available
        and if none is available will default back to CPU; in 
        general, it is recommended not to try and set this unless
        there's a specific reason why a specific GPU should be
        used. Default = 0.

    show_progress_bar : bool
        Flag which, if set to True, means a progress bar is printed as 
        predictions are made, while if False no progress bar is printed.
        Default  =  True

    batch_mode : int
        Indictora which, if set to 1 or 2 will FORCE the batch 
        algorithm to use mode 1 or mode 2 for batch 
        decomposition.

        Mode 1 means we pre-filter sequences into groups where 
        they're all the same length, avoiding padding/packing. 
        This works in all versions of torch, and will be faster
        if you have very large datasets or have many copies of 
        the same sequence.

        Mode 2 involves padding/packing the sequences so that 
        all sequences can be passed in a batchsize of 32. This 
        is only available if pytorch 1.11 or higher is available, 
        but for small sets of sequences 1-10,000 will be much 
        faster than mode 1. We default to mode 2 if available, 
        but in special cases you may want to force mode 1.

        Default = None, which means dynamic selection occurs (2
        if available, fall-back to 1). However 1 may often actually
        be more efficient, so it's worth testing modes to see if
        there's any change in perforance for a given dataset.


    safe : bool
        Flag which, if set to False, means the function overwrites 
        existing tracks and domains if present. If True, overwriting
        will trigger an exception.
        Default = True.
        
    Returns
    -----------------
    None
        No return type, but the Protein objects in the Proteome 
        will be annotated with per-residue disorder Tracks.

    """
    uid2seq = {}
    for p in proteome:
        uid2seq[p.unique_ID] = p.sequence

    # batch predict disorder
    D = meta.predict_disorder_batch(uid2seq, gpuid=gpuid, show_progress_bar=show_progress_bar, batch_mode=batch_mode)

    for k in uid2seq:
        proteome.protein(k).add_track(name, values=D[k][1], safe=safe)

    


def annotate_proteome_with_disordered_domains(proteome,
                                              name='IDR',
                                              disorder_threshold=0.5,
                                              annotate_folded_domains=False,
                                              folded_domain_name = 'FD',
                                              gpuid=00,
                                              show_progress_bar=True,
                                              batch_mode=None,                                              
                                              safe=True):
    """
    Function that annotates a proteome with disordered
    Domains (IDRs) for every protein. 

    By default, disordered Domains are named as 'IDR's, although
    this can be changed by setting the `name` parameter.

    In addition, if requested, folded domains can also be annotated
    as those domains which are not IDRs. These folded domains are
    named 'FD's by default, although this can be changed by setting
    the `folded_domain_name` parameter.
    
    Disorder prediction uses the batch mode in metapredict, which 
    leverages parallel predictions automatically on GPUs or CPUs.
    However, if a specific device is requested this can be passed

    Parameters
    -----------------
    proteome : shephard.proteome.Proteome 
        Proteome object to be annotated.

    disorder_threshold : float
        Threshold to be used to define IDRs by the metapredict
        domain decomposition algorithm. The default is 0.5, 
        and we strongly recommend sticking with this value.

    track_name : str
        Name of the Track added to each Protein. 
        Default = 'disorder'

    domain_name : str
        Name of the Domain added to each Protein. 
        Default = 'IDR'

    annotate_folded_domains : bool
        Flag which, if included, means we ALSO annotate 
        the regions that are not IDRs as 'FD' (folded
        domains), where the name can be changed using
        the folded_domain_name variable.
        Default = False

    folded_domain_name : str
        String used to name Folded Domains. Only relevant
        if annotate_folded_domains is set to True.
        Default = 'FD'

    gpuid : int 
        Identifier for the GPU being requested. Note that if
        this is left unset the code will use the first GPU available
        and if none is available will default back to CPU; in 
        general, it is recommended not to try and set this unless
        there's a specific reason why a specific GPU should be
        used. Default = 0.

    show_progress_bar : bool
        Flag which, if set to True, means a progress bar is printed as 
        predictions are made, while if False no progress bar is printed.
        Default  =  True

    batch_mode : int
        Indictora which, if set to 1 or 2 will FORCE the batch 
        algorithm to use mode 1 or mode 2 for batch 
        decomposition.

        Mode 1 means we pre-filter sequences into groups where 
        they're all the same length, avoiding padding/packing. 
        This works in all versions of torch, and will be faster
        if you have very large datasets or have many copies of 
        the same sequence.

        Mode 2 involves padding/packing the sequences so that 
        all sequences can be passed in a batchsize of 32. This 
        is only available if pytorch 1.11 or higher is available, 
        but for small sets of sequences 1-10,000 will be much 
        faster than mode 1. We default to mode 2 if available, 
        but in special cases you may want to force mode 1.

        Default = None, which means dynamic selection occurs (2
        if available, fall-back to 1). However 1 may often actually
        be more efficient, so it's worth testing modes to see if
        there's any change in perforance for a given dataset.

    safe : bool
        Flag which, if set to False, means the function overwrites 
        existing tracks and domains if present. If True, overwriting
        will trigger an exception.
        Default = True.
        
    Returns
    -----------------
    None
        No return type, but the Protein objects in the Proteome 
        will be annotated with disordered Domain annotations.
        
    """
    
    uid2seq = {}
    for p in proteome:
        uid2seq[p.unique_ID] = p.sequence

    # batch predict disorder
    D = meta.predict_disorder_batch(uid2seq, return_domains=True, gpuid=gpuid, show_progress_bar=show_progress_bar, batch_mode=batch_mode)

    for k in uid2seq:

        X = D[k]

        for boundaries in X.disordered_domain_boundaries:
            proteome.protein(k).add_domain(boundaries[0]+1, boundaries[1], name, safe=safe)

        if annotate_folded_domains:
            for boundaries in X.folded_domain_boundaries:
                proteome.protein(k).add_domain(boundaries[0]+1, boundaries[1], folded_domain_name, safe=safe)
            



def annotate_proteome_with_disorder_tracks_and_disordered_domains(proteome,
                                                                  disorder_threshold=0.5,
                                                                  track_name='disorder',
                                                                  domain_name='IDR',
                                                                  annotate_folded_domains=False,
                                                                  folded_domain_name = 'FD',
                                                                  gpuid=00,
                                                                  show_progress_bar=True,
                                                                  batch_mode=None,
                                                                  safe=True):
    """
    Function that annotates a proteome with disorder Tracks and
    disorder Domains for every protein. 

    By default, disorder Tracks are named 'disoder', although
    this can be changed by setting the `track_name` parameter.

    By default, disordered Domains are named as 'IDR's, although
    this can be changed by setting the `name` parameter.

    In addition, if requested, folded domains can also be annotated
    as those domains which are not IDRs. These folded domains are
    named 'FD's by default, although this can be changed by setting
    the `folded_domain_name` parameter.

    Disorder prediction uses the batch mode in metapredict, which 
    leverages parallel predictions automatically on GPUs or CPUs.
    However, if a specific device is requested this can be passed

    Parameters
    -----------------
    proteome : shephard.proteome.Proteome 
        Proteome object to be annotated.

    disorder_threshold : float
        Threshold to be used to define IDRs by the metapredict
        domain decomposition algorithm. Default is 0.5 and strongly
        recommend sticking with this value.

    track_name : str
        Name of the Track added to each Protein. 
        Default = 'disorder'

    domain_name : str
        Name of the Domain added to each Protein. 
        Default = 'IDR'

    annotate_folded_domains : bool
        Flag which, if included, means we ALSO annotate 
        the regions that are not IDRs as 'FD' (folded
        domains), where the name can be changed using
        the folded_domain_name variable.
        Default = False

    folded_domain_name : str
        String used to name Folded Domains. Only relevant
        if annotate_folded_domains is set to True.
        Default = 'FD'

    gpuid : int 
        Identifier for the GPU being requested. Note that if
        this is left unset the code will use the first GPU available
        and if none is available will default back to CPU; in 
        general, it is recommended not to try and set this unless
        there's a specific reason why a specific GPU should be
        used. Default = 0.

    show_progress_bar : bool
        Flag which, if set to True, means a progress bar is printed as 
        predictions are made, while if False no progress bar is printed.
        Default  =  True

    batch_mode : int
        Indictora which, if set to 1 or 2 will FORCE the batch 
        algorithm to use mode 1 or mode 2 for batch 
        decomposition.

        Mode 1 means we pre-filter sequences into groups where 
        they're all the same length, avoiding padding/packing. 
        This works in all versions of torch, and will be faster
        if you have very large datasets or have many copies of 
        the same sequence.

        Mode 2 involves padding/packing the sequences so that 
        all sequences can be passed in a batchsize of 32. This 
        is only available if pytorch 1.11 or higher is available, 
        but for small sets of sequences 1-10,000 will be much 
        faster than mode 1. We default to mode 2 if available, 
        but in special cases you may want to force mode 1.

        Default = None, which means dynamic selection occurs (2
        if available, fall-back to 1). However 1 may often actually
        be more efficient, so it's worth testing modes to see if
        there's any change in perforance for a given dataset.

    safe : bool
        Flag which, if set to False, means the function overwrites 
        existing tracks and domains if present. If True, overwriting
        will trigger an exception.
        Default = True.
        
    Returns
    -----------------
    None
        No return type, but the Protein objects in the Proteome 
        will be annotated with per-residue disorder Tracks and
        disordered Domain annotations.

    """
    
    uid2seq = {}
    for p in proteome:
        uid2seq[p.unique_ID] = p.sequence

    # batch predict disorder annotations/scores
    D = meta.predict_disorder_batch(uid2seq, return_domains=True, gpuid=gpuid, show_progress_bar=show_progress_bar, batch_mode=batch_mode)

    # for each unique ID
    for k in uid2seq:

        # X = DisorderObject
        X = D[k]

        # cycle through IDR boundaries
        for boundaries in X.disordered_domain_boundaries:
            proteome.protein(k).add_domain(boundaries[0]+1, boundaries[1], domain_name, safe=safe)

        if annotate_folded_domains:
            for boundaries in X.folded_domain_boundaries:
                proteome.protein(k).add_domain(boundaries[0]+1, boundaries[1], folded_domain_name, safe=safe)


        proteome.protein(k).add_track(track_name, values=X.disorder, safe=safe)

