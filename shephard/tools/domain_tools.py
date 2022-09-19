"""
SHEPHARD: 
Sequence-based Hierarchical and Extendable Platform for High-throughput Analysis of Region of Disorder

Authors: Garrett M. Ginell & Alex S. Holehouse
Contact: (g.ginell@wustl.edu)

Holehouse Lab - Washington University in St. Louis
"""

import numpy as np
import random
from . import site_tools
from shephard import general_utilities
from shephard import exceptions
from shephard.interfaces.interface_tools import check_protein

## ------------------------------------------------------------------------
##
def domain_overlap(domain_1, domain_2, check_origin=True):
    """
    Given two domains asks if their boundaries overlap. By default
    this expects the two domains in question to be from the same protein
    and checks this. If we dont want to enforce this assumption set
    check_origin to False.

    Parameters
    -----------
    domain_1 : shephard.domain.Domain
        The first domain object of interest

    domain_2 : shephard.domain.Domain
        The first domain object of interest

    check_origin : bool
        Flag that if set to True will cause an exception if domain_1 and
        domain_2 are from different proteins. If set to false, no such
        sanity checks are performed.

    Returns
    ----------
    boolean
        Returns true if the two domains overlap, else returns false
      
    """

    if check_origin:
        if domain_1.protein.unique_ID != domain_2.protein.unique_ID:
            raise exceptions.DomainException('Examining overlap of %s and %s but these are from different proteins' % (str(domain_1), str(domain_2)))
            
    return domain_overlap_by_position(domain_1.start, domain_1.end, domain_2.start, domain_2.end)



## ------------------------------------------------------------------------
##
def domain_overlap_fraction(domain_1, domain_2, check_origin=True):
    """
    Given two domains asks what fraction the shorter domain 
    overlaps the longer one with.
    Parameters
    -----------
    domain_1 : shephard.domain.Domain
        The first domain object of interest

    domain_2 : shephard.domain.Domain
        The first domain object of interest

    check_origin : bool (default = True)
        Flag that if set to True will cause an exception if domain_1 and
        domain_2 are from different proteins. If set to false, no such
        sanity checks are performed.

    Returns
    ----------
    float
        Returns a float between 0 and 1 that corresponds to what 
        fraction of the shorter domain overlaps with the longer domain.
      
    """

    if check_origin:
        if domain_1.protein.unique_ID != domain_2.protein.unique_ID:
            raise exceptions.DomainException('Examining overlap of %s and %s but these are from different proteins' % (str(domain_1), str(domain_2)))


    if len(domain_1) < len(domain_2):
        d_short = domain_1
        d_long = domain_2
    else:
        d_short = domain_2
        d_long = domain_1

    
    # ......OOOOOOOOOOOOOOOOOOOOOOO............       long
    #...XXXXXXXXXXXX?????                             short
    #
    if d_short.start < d_long.start:
        if d_short.end < d_long.start:
            return 0.0
        else:
            return ((d_short.end - d_long.start)+1)/len(d_short)

    # ......OOOOOOOOOOOOOOOOOOOOOOO.............     end
    #.                   ?????XXXXXXXXXXXX......     short
    #
    if d_short.end > d_long.end:
        if d_short.start > d_long.end:
            return 0.0
        else:

            # note - the +1 accounts for the fact we're using 
            # an inclusive counting 
            return ((d_long.end - d_short.start)+1)/len(d_short)

    # if we get here d_short.start equal to or larger than d_long start
    # and d_short end equal to or smaller than d_long end, so 100% overlap
    return 1.0
    


## ------------------------------------------------------------------------
##
def domain_overlap_by_position(boundary_start1, boundary_end1, boundary_start2, boundary_end2):
    """
    Given four sets of starting/ending positions, this function asks if 
    their boundaries overlap. 

    Parameters
    -----------
    boundary_start1 : int
        Position of domain 1 start

    boundary_end : int
        Position of domain 1 end

    boundary_start2 : int
        Position of domain 2 start

    boundary_end : int
        Position of domain 2 end

    Returns
    ----------
    boolean
        Returns true if the two domains overlap, else returns false
    """

    # note we do this swapping around which of the two domains passed as boundaries is 'a' and 'b' in 
    # the visual schematic below
    for x in [[boundary_start1, boundary_end1, boundary_start2, boundary_end2], [boundary_start2, boundary_end2,boundary_start1, boundary_end1]]:
        
        a_start = x[0]
        a_end   = x[1]
        b_start = x[2]
        b_end   = x[3]


        # scenario 1
        # AAAAAAAA
        #     BBBBBBBBB
        if a_start <= b_start and a_end >= b_start:
            return True

        # scenario 2
        #           AAAAAAAA
        #     BBBBBBBBB
        if a_start >= b_start and a_start <= b_end:
            return True

    return False
    



## ------------------------------------------------------------------------
##
def build_missing_domains(protein, new_domain_type = 'missing'):
    """
    Function which takes a protein and builds a set of domains that represent 
    the "empty spaces". Domains are returned as a list of domain dictionaries
    which can be added to a protein via the add_domains() function.

    This tool is stateless - i.e. it does not alter the passed protein
    but instead only generates a numerical list which could be added as 
    a track.

    One could always combine this directly with the add_domains() function
    into a single line - e.g.
    
    # this line will automatically add all the missing regions as 
    # 'missing' proteins to the domain
    protein.add_domains(build_missing_domains(protein))

    Parameters
    -----------
    protein : shephard.protein.Protein object
        Protein object over which sites are identified

    new_domain_type : str (default = 'missing')
        Name to assign to the 'empty' domains.

    Returns
    -------------
    list of domain dictionaries

        Returns a list of domain dictionaries which can be then parsed or
        added to a protein via the add_domains() function.
    
    """

    check_protein(protein, 'build_missing_domains()')

    ##
    ## NOTE -this function uses i0 indexing to keep things simple and then
    ## corrects at the end
    ##
    
    # first constrct an empty vector of 0s. We're going to build '1s' into
    # the positions in the sequence occupied by domains
    all_res = [0]*len(protein)
    
    # for each domain in the protein
    for d in protein.domains:

        # for each position in each domain (note this will overwrite
        # which is fine). NOTE that we cyle from d.start-1 to d.end because
        # start and end index from real-world positions, but range is a non
        # inclusive function (while the domain boundaries are inclusive). 
        for i in range(d.start-1, d.end):
            all_res[i] = 1
                
    # we are now left with an all_res list where the 0s represent regions
    # not assigned to any domain

    # if we start inside an empty domain set inside to true
    if all_res[0] == 0:
        start = 0
        inside=True

    # else set start to -1  and inside to false
    else:
        start = -1
        inside=False

    all_domains=[]

    # for each position in the sequence (starting at 1 because 
    # we already assessed 0 above. 
    for i in range(1,len(all_res)):

        # if we find ourselves in an empty space
        if all_res[i] == 0:
            if not inside:
                inside=True
                start = i
        # if we fine ourselves in a domain of some kind
        else:

            # if we were previously in an empty space
            if inside:
                inside = False

                # not here we're adding in i0 indexing
                all_domains.append([start, i-1])

    # if we ended and were inside when we ended then the last
    # empty domain spans to the final residue in the sequence
    if inside:
        all_domains.append([start,len(all_res)-1])

    new_domain_list =[]
    for d in all_domains:

        # note the +1 here is so we insert at residue positions
        # that correctly map to real-space positions (i.e.i1 indexing)

        # construct a domain dictionary 
        local_domain={}        
        local_domain['start'] = d[0]+1
        local_domain['end'] = d[1]+1
        local_domain['domain_type'] = new_domain_type

        new_domain_list.append(local_domain)

    return new_domain_list
            
        

## ------------------------------------------------------------------------
##
def build_domains_from_track_values(proteome, 
                                    track_name, 
                                    binerize_function, 
                                    domain_type, 
                                    gap_closure = 3, 
                                    minimum_region_size = 20, 
                                    extend_ends = None, 
                                    verbose = True):

    
    """
    Function which takes a Proteome and builds a set of domains based on
    values tracks in each Protein in that Proteome. This effectively 
    allows you to discretize some continous variable into distinct local 
    domains, which can often facilitate specific types of analysis. This
    conversion is done using a custom-passed binerize function which 
    converts a normal track into a track of 0s and 1s. Residues that
    are assigned a value of 1 will be included in a domain assuming they
    fall within a contigous region of sufficient size, as defined by
    the parameters gap_closure and minimum_region_size, as discussed 
    below.

    This function operates on an entire Proteome-level, and is stateless
    (i.e. does not directly alter the passed proteome). Instead, the 
    function dictionary where keys are unique_IDs of proteins and values 
    is a list of one or more Domain dictioinaries (with a start, end,
    and domain_type key:value pair).

    The domains dictionary can be added to a proteome using the 
    si_domains.add_domains_from_dictionary(). As an example, as possible 
    workflow is as follows:

    >>> d = build_domains_from_track_values(proteome, 'cool_track', trackfx)
    >>> si_domains.add_domains_from_dictionary(proteome, d)
    
    Under the hood, the function works by cycling through each protein, 
    extracting the track, and converting into domains.
        
    If a protein is too short or it lacks a given track, the protein is 
    skipped.
    
    Parameters
    -----------
    proteome : shephard.proteome.Proteome
        The Proteome which is going to be scanned for each track. Note that
        the underlying Proteome is not altered by this function

    track_name : string
        Name of the track to convert. If the track name does not exist in
        a given protein that protein is skipped. In this way, a Proteome
        where only a subset of Proteins have tracks can be parsed without
        issue. The track must be a values track - symbols tracks should be
        converted to a values track first to avoid issue.

    binerize_function : function
        A function which takes a track and converts it to 0 or 1 (binerize, 
        as in, make binary). This enables a complex and continous 
        track to be converted into a binary classification, which is 
        practically what a domain-assigment needs (yes/no inside domain).
        This function must take in a single variable (the track values) 
        and return a new list or numpy array that is the same length as 
        the track values but possesses only 0 and 1 in each element.

    domain_type : str
        String that defines the name of the new domains to create. Can in
        principle be anything.

    gap_closure : int (default = 3)
        Defines spacing between 1s or 0s that will be filled in to generate
        contigous stretches of 0s or 1s. This helps avoid a scenario where
        breaks in contigous stretches impede the definition of a domain 
        above a certain size, as defined by minimum_region_size. In general
        a value of 3 works reasonably well in most scenarios.

    minimum_region_size : int (default = 20)
        Defines the smallest size for a domain allowed. This can be varied
        depending on the question or data, and it may make sense to have 
        corresponding changes in gap_closure if this value becomes 
        substantially larger than a gap_closure of 3.

    extend_ends : int (default = None)
        This is a somewhat niche feature which, if set to a number, means 
        that we check the extend_ends-th value at the N- and C-terminus 
        of the binarized track, and if 1 set all values from that position
        to the N and/or C terminus to 1. This is provided because sometimes
        binerize functions will inherently struggle with the very ends of 
        sequences, so this provides a way to cast the first and last 
        extend_ends values to be 1. This is fairly specific and probably
        only worth using in a scenario where there is a clear issue

    verbose : bool (default = True)
        This flag enables the function to print statues every 500 proteins.
        If the binerize function is expensive this can be good to ensure 
        progress is proceeding.

    Returns
    -------------
    dict
        Returns a dictionary of key-value pairs, where each key is a unique
        ID and each value is a list of 1 or more domain dictionaries. This
        return dictionary can be directly added to a Proteome using the 
        Proteome.add_domains_from_dictonary() function.

    """
    
    new_domains = {}

    c = 0
    for protein in proteome:

        # if our protein is too short do not try and generate a domain
        if len(protein) < 3*gap_closure + 1:
            continue
        
        # this is the counter of proteins we've actually scanned - if
        # we hit a 500-protein milestone print status if verbose is true
        c = c + 1
        if verbose and c % 500 == 0:
            print('On %i of %i' %(c, len(proteome)))
        
        # safe = false so will return None if no track of that name
        # found, and if so we continue to next protein
        t = protein.track(track_name, safe=False)
        if t is None:
            continue

        # extract out the values
        t = t.values
        
        # and apply the binerize function to the values
        B = binerize_function(t)

        # add to help debugging
        if len(B) != len(protein):
            raise exceptions.DomainException('Binerize function failed to return an array or list that matches the protein sequence length')
            

        ## Part 1 - remove gapes
        for g in range(1, gap_closure+1):
            # first fill in 


            # for each position
            i = 0
            finished = False
            while not finished:

                p1 = i
                p2 = i + g
                p3 = i + 2*g
                p4 = i + 3*g

                # if the complete set of smaller regions ahead is empty or 
                # fully assigned skip ahead because nothing to do here...
                if np.sum(B[p1:p4]) == 0:
                    i = p4

                elif np.sum(B[p1:p4]) == 3*g:

                    # we jump to the p3 position (and NOT p4) as this allows us to skip along without
                    # discarding positions we need for filling. Note if we know everything is empty
                    # it doesn't matter and we can jump to p4
                    i = p3

                else:
                    # if we have gapsize number of hits
                    if np.sum(B[p1:p2])  == g:

                        # and if a gap away there is another gapsize 
                        if general_utilities.numerical_sum(B[p3:p4])  == g:                        
                            B[p2:p3] = [1]*g
                    
                    i = i + 1

                if i + 3*g >= len(B):
                    finished = True


        ## Part 2 - remove domains that are too small - we adde the '-' caps so we can use
        # replace and distinguish c/n terminal values
        B_string = '-'
        for i in B:
            if i == 1:
                B_string = B_string + "1"
            else:
                B_string = B_string + "0"

        B_string=B_string+'-'
        
        # for sizes of contigous stretches that are up minimum_region_size + 1
        # replace with empty ('0') strings
        for i in range(1, minimum_region_size + 1):      

            # 011110 -> 000000
            B_string = B_string.replace('0' + i*'1' + '0', '0' + i*'0' + '0')

            # -11110 -> -00000
            B_string = B_string.replace('-'+i*'1' + '0', '-'+i*'0' + '0')

            # 01111- -> 00000-
            B_string = B_string.replace('0' + i*'1'+'-' , '0'+ i*'0'+'-' )

            # -1111- -> -0000-
            B_string = B_string.replace('-' + i*'1'+'-' , '-'+ i*'0'+'-' )


        # 1 to -1 to cut off the artifical caps we added
        for i in range(1, len(B_string)-1):
            B[i-1] = int(B_string[i])

        ## Part 3 - extend ends if required
        if extend_ends:
                
            # note we don't check length of extend ends initially, so if they're too big for the sequence
            # dynamicaly resize them so they're 1/2 of the sequence length
            if extend_ends + 1 >= len(B):
                extend_ends_val = len(B)/2
            else:
                extend_ends_val = extend_ends
                
            if B[extend_ends_val + 1] == 1:
                B[0:extend_ends_val] = [1]*len(B[0:extend_ends_val])

            if B[:-(extend_ends_val + 1)] == 1:
                B[-extend_ends_val:] = [1]*len(B[0:extend_ends_val])
            
        ## Part 4 - extract domain boundaires
        local_domains=[]
        if B[0] == 1:
            inside = True
            start = 0
        else:
            inside = False
            
        # for each position
        for idx in range(0,len(B)):

            i = B[idx]
                
            if i == 1:
                if inside:
                    continue
                else:
                    inside = True
                    start = idx

            if i == 0:
                if inside:
                    inside = False
                    end = idx-1
                    local_domains.append({'start':start+1,'end':end+1,'domain_type':domain_type})

        # if we finished inside
        if inside:
            local_domains.append({'start':start+1,'end':len(B),'domain_type':domain_type})

        # if we found any local domains...
        if len(local_domains) > 0:
            new_domains[protein.unique_ID] = local_domains

    return new_domains


   
# used for test example - kept in case we want to expand or update this in          

# def calculate_domain_length(domain):
#    """
#    Function that returns a domain's length
#    
#    Parameters
#    -------------
#    domain : Domain
#    Domain object in question
#
#    Returns
#    ------------
#    int
#    Returns the length of the domain
#    """
#    
#    return len(domain)


                
            


