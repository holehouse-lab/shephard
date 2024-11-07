##
## API into metapredict
##
## For all APIS we do do not make hard dependencies, safe informative
## import checking should be done.
##


## Check metapredict is installed
try:
    import metapredict as meta

    try:
        meta.predict_disorder('AAPAPA',version=3)
    except TypeError:
        print('shephard requires metapredict V3 or higher. Please upgrade:')
        print('pip install --upgrade metapredict')
    
except ModuleNotFoundError:
    print('Unable to import metapredict')
    print('To use the metapredict API, make sure metapredict is installed')
    print('This can be done as follows:')
    print('pip install metapredict')

## then check batch mode is available

    
## ------------------------------------------------------------------------
##
def annotate_proteome_with_disorder_track(proteome,                                       
                                          name='disorder',
                                          device=None,
                                          version=3,                                          
                                          show_progress_bar=True,
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

    device : int or str 
        Identifier for the device to be used for predictions. 
        Possible inputs: 'cpu', 'mps', 'cuda', or an int that corresponds to
        the index of a specific cuda-enabled GPU. If 'cuda' is specified and
        cuda.is_available() returns False, instead of falling back to CPU, 
        metapredict will raise an Exception so you know that you are not
        using CUDA as you were expecting. 
        Default: None
            When set to None, we will check if there is a cuda-enabled
            GPU. If there is, we will try to use that GPU. 
            If you set the value to be an int, we will use cuda:int as the device
            where int is the int you specify. The GPU numbering is 0 indexed, so 0 
            corresponds to the first GPU, and so on. Only specify this if you
            know which GPU you want to use. 
            * Note: MPS is only supported in Pytorch 2.1 or later. 
            MPS is still fairly new, so use it at your own risk.  
                                          
    version : int
        Defines the metapredict version to use (must be one of 1, 2 
        or 3).

    show_progress_bar : bool
        Flag which, if set to True, means a progress bar is printed as 
        predictions are made, while if False no progress bar is printed.
        Default  =  True

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
    try:
        D = meta.predict_disorder(uid2seq, device=device, show_progress_bar=show_progress_bar, version=version)
    except:
        D = meta.predict_disorder(uid2seq, device=device, show_progress_bar=show_progress_bar, version=version)

    for k in uid2seq:
        proteome.protein(k).add_track(name, values=D[k][1], safe=safe)

    

## ------------------------------------------------------------------------
##
def annotate_proteome_with_disordered_domains(proteome,
                                              name='IDR',
                                              disorder_threshold=0.5,
                                              annotate_folded_domains=False,
                                              folded_domain_name = 'FD',
                                              device=None,
                                              version=3,                                                                                        
                                              show_progress_bar=True,
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

    name : str
        Name to give IDR domains.

    disorder_threshold : float
        Threshold to be used to define IDRs by the metapredict
        domain decomposition algorithm. The default is 0.5, 
        and we strongly recommend sticking with this value.

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

    device : int or str 
        Identifier for the device to be used for predictions. 
        Possible inputs: 'cpu', 'mps', 'cuda', or an int that corresponds to
        the index of a specific cuda-enabled GPU. If 'cuda' is specified and
        cuda.is_available() returns False, instead of falling back to CPU, 
        metapredict will raise an Exception so you know that you are not
        using CUDA as you were expecting. 
        Default: None
            When set to None, we will check if there is a cuda-enabled
            GPU. If there is, we will try to use that GPU. 
            If you set the value to be an int, we will use cuda:int as the device
            where int is the int you specify. The GPU numbering is 0 indexed, so 0 
            corresponds to the first GPU and so on. Only specify this if you
            know which GPU you want to use. 
            * Note: MPS is only supported in Pytorch 2.1 or later. 
            MPS is still fairly new, so use it at your own risk.  

    version : int
        Defines the metapredict version to use (must be one of 1, 2 
        or 3).

    show_progress_bar : bool
        Flag which, if set to True, means a progress bar is printed as 
        predictions are made, while if False no progress bar is printed.
        Default  =  True

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
    D = meta.predict_disorder(uid2seq, device=device, show_progress_bar=show_progress_bar, version=version, return_domains=True)

    for k in uid2seq:

        X = D[k]

        for boundaries in X.disordered_domain_boundaries:
            proteome.protein(k).add_domain(boundaries[0]+1, boundaries[1], name, safe=safe)

        if annotate_folded_domains:
            for boundaries in X.folded_domain_boundaries:
                proteome.protein(k).add_domain(boundaries[0]+1, boundaries[1], folded_domain_name, safe=safe)



## ------------------------------------------------------------------------
##
def annotate_proteome_with_disorder_tracks_and_disordered_domains(proteome,
                                                                  track_name='disorder',
                                                                  domain_name='IDR',
                                                                  disorder_threshold=0.5,
                                                                  annotate_folded_domains=False,
                                                                  folded_domain_name = 'FD',
                                                                  device=None,
                                                                  version=3,                                                                                        
                                                                  show_progress_bar=True,
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

    track_name : str
        Name of the Track added to each Protein. 
        Default = 'disorder'

    domain_name : str
        Name of the Domain added to each Protein. 
        Default = 'IDR'

    disorder_threshold : float
        Threshold to be used to define IDRs by the metapredict
        domain decomposition algorithm. Default is 0.5 and strongly
        recommend sticking with this value.

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

    device : int or str 
        Identifier for the device to be used for predictions. 
        Possible inputs: 'cpu', 'mps', 'cuda', or an int that corresponds to
        the index of a specific cuda-enabled GPU. If 'cuda' is specified and
        cuda.is_available() returns False, instead of falling back to CPU, 
        metapredict will raise an Exception so you know that you are not
        using CUDA as you were expecting. 
        Default: None
            When set to None, we will check if there is a cuda-enabled
            GPU. If there is, we will try to use that GPU. 
            If you set the value to be an int, we will use cuda:int as the device
            where int is the int you specify. The GPU numbering is 0 indexed, so 0 
            corresponds to the first GPU and so on. Only specify this if you
            know which GPU you want to use. 
            * Note: MPS is only supported in Pytorch 2.1 or later. 
            MPS is still fairly new, so use it at your own risk.             

    version : int
        Defines the metapredict version to use (must be one of 1, 2 
        or 3).

    show_progress_bar : bool
        Flag which, if set to True, means a progress bar is printed as 
        predictions are made, while if False no progress bar is printed.
        Default  =  True

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
    D = meta.predict_disorder(uid2seq, device=device, show_progress_bar=show_progress_bar, version=version, return_domains=True)

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

