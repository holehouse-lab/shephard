##
## API into ALBATROSS
##
## For all APIS we do do not make hard dependencies, safe informative
## import checking should be done.
##


## Check metapredict is installed
try:
    from sparrow.predictors.batch_predict import batch_predict
    
except ModuleNotFoundError:
    print('Unable to import sparrow (the package where ALBATROSS is implemeted)')
    print('To use the ALBATROSS, make sure sparrow is installed')
    print('This can be done as follows:')
    print('pip install git+https://git@github.com/idptools/sparrow.git')

    
## ------------------------------------------------------------------------
##
def annotate_proteome_with_dimensions(proteome,                                        
                                      rg_name = 'rg',
                                      re_name = 're',
                                      gpuid=00,
                                      show_progress_bar=True,
                                      batch_mode=None,
                                      safe=True):
    """
    Function that annotates a proteome with it's predicted radius
    of gyration (rg) and end-to-end distance (re) for every protein.

    By default, rg and re are added as attributes to each Protein, with
    the names 'rg' and 're' respectively. However, this can be changed
    by setting the `rg_name` and `re_name` parameters.

    Dimension prediction uses the batch mode in sparrow, which
    leverages parallel predictions automatically on GPUs or CPUs.
    However, if a specific device is requested, this can be passed via
    the `gpuid` parameter.

    Parameters
    -----------------
    proteome : shephard.proteome.Proteome 
        Proteome object to be annotated.

    rg_name : str
        Name of the rg attribute added to each Protein.

    re_name : str
        Name of the re attribute added to each Protein.

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

    # batch predict dimensions for all proteins
    rg = batch_predict(uid2seq, network='scaled_rg', gpuid=gpuid, show_progress_bar=show_progress_bar)
    re = batch_predict(uid2seq, network='scaled_re', gpuid=gpuid, show_progress_bar=show_progress_bar)

    # add as an attribute to the proteins
    for k in rg:
        proteome.protein(k).add_attribute(rg_name, rg[k][1], safe=safe)
        
    for k in re:
        proteome.protein(k).add_attribute(re_name, re[k][1], safe=safe)


        
## ------------------------------------------------------------------------
##
def annotate_domains_with_dimensions(proteome,
                                     domain_type,
                                     rg_name = 'rg',
                                     re_name = 're',
                                     gpuid=00,
                                     show_progress_bar=True,
                                     batch_mode=None,
                                     safe=True):
    """
    Function that annotates every domain matching the domain_name in
    a proteome with it's predicted radius of gyration (rg) and end-to-end
    distance (re).

    By default, rg and re are added as attributes to each Domain, with
    the names 'rg' and 're' respectively. However, this can be changed
    by setting the `rg_name` and `re_name` parameters.

    Dimension prediction uses the batch mode in sparrow, which
    leverages parallel predictions automatically on GPUs or CPUs.
    However, if a specific device is requested, this can be passed via
    the `gpuid` parameter.

    Parameters
    -----------------
    proteome : shephard.proteome.Proteome 
        Proteome object to be annotated.

    domain_type : str
        Type of the domain to be annotated.

    rg_name : str
        Name of the rg attribute added to each Protein.

    re_name : str
        Name of the re attribute added to each Protein.

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

    # build the dictionary of unique IDs to sequences
    uid2seq = {}
    for d in proteome.domains:
        if d.domain_type == domain_type:
            unique_name = d.protein.unique_ID + '_' + d.domain_name
            uid2seq[unique_name] = d.sequence

    # batch predict dimensions for all proteins
    rg = batch_predict(uid2seq, network='scaled_rg', gpuid=gpuid, show_progress_bar=show_progress_bar)
    re = batch_predict(uid2seq, network='scaled_re', gpuid=gpuid, show_progress_bar=show_progress_bar)

    for d in proteome.domains:
        if d.domain_type == domain_type:
            
            unique_name = d.protein.unique_ID + '_' + d.domain_name
            rg_val = rg[unique_name][1]
            re_val = re[unique_name][1]

            d.add_attribute(rg_name, rg_val, safe=safe)
            d.add_attribute(re_name, re_val, safe=safe)
    
        
