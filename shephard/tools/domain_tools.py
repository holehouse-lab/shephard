import numpy as np
import random
from . import site_tools
from shephard import general_utilities

def build_missing_domains(protein, new_domain_type = 'missing'):
    """
    Function which takes a protein an builds a set of domains that represent 
    the "empty spaces". Domains are returned as a list of domain dictionaries
    which can be added to a protein via the add_domains() function.

    One could always combine this directly with the add_domains() function
    into a single line - e.g.
    
    # this line will automatically add all the missing regions as 
    # 'missing' proteins to the domain
    protein.add_domains(build_missing_domains(protein))

    Parameters
    -----------
    protein : Protein object
        Protein object over which sites are identified

    new_domain_type : string (default 'missing')
        Name to assign to the 'empty' domains.

    Returns
    -------------
    list of domain dictionaries

        Returns a list of domain dictionaries which can be then parsed or
        added to a protein via the add_domains() function.
    
    """
    
    # first constrct an empty vector of 0s. We're going to build '1s' into
    # the positions in the sequence occupied by domains
    all_res = [0]*len(protein)
    
    # for each domain in the protein
    for d_idx in protein.domains:
        d = protein.domain(d_idx)

        # for each position in each domain (not this will overwrite
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
            
        


def build_domains_from_track_values(proteome, track_name, binerize_function, domain_type, gap_closure=3, minimum_region_size=1, extend_ends=None, print_progress=False):
    
    new_domains={}

    c=0
    for protein in proteome:

        if len(protein) < 3*gap_closure+1:
            continue
        
        c=c+1
        if print_progress and c % 500 == 0:
            print('On %i of %i' %(c, len(proteome)))
        

        # safe = false so will return None if no track of that name
        # found
        t = protein.track(track_name, safe=False)
        
        # if no track found continue to next protein
        if t is None:
            continue

        # extract out the values
        t = t.values
        
        # and apply the binerize function to the values
        B = binerize_function(t)

        ## Part 1 - remove gapes
        for g in range(1,gap_closure+1):
            # first fill in 


            # for each position
            i=0
            finished = False
            while not finished:


                p1 = i
                p2 = i + g
                p3 = i + 2*g
                p4 = i + 3*g

                # if the complete set of smaller regions ahead is empty or 
                # fully assigned skip ahead because nothing to do here...
                #if general_utilities.numerical_sum(B[p1:p4]) == 0:
                if np.sum(B[p1:p4]) == 0:
                    i = p4

                #elif general_utilities.numerical_sum(B[p1:p4]) == 3*g:
                elif np.sum(B[p1:p4]) == 3*g:

                    # we jump to the p3 position (and NOT p4) as this allows us to skip along without
                    # discarding positions we need for filling. Note if we know everything is empty
                    # it doesn't matter and we can jump to p4
                    i = p3

                else:
                    # if we have gapsize number of hits
                    #if general_utilities.numerical_sum(B[p1:p2])  == g:
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
                B_string=B_string+"1"
            else:
                B_string=B_string+"0"

        B_string=B_string+'-'
        
        for i in range(1, minimum_region_size):                    
            B_string = B_string.replace('0' + i*'1' + '0', '0' + i*'0' + '0')
            B_string = B_string.replace('-'+i*'1' + '0', '-'+i*'0' + '0')
            B_string = B_string.replace('0' + i*'1'+'-' , '0'+ i*'1'+'-' )


        # 1 to -1 to cut off the artifical caps we added
        for i in range(1, len(B_string)-1):
            B[i-1] = int(B_string[i])



        ## Part 3 - extend ends if required
        if extend_ends:
                
            # note we don't check length of extend ends initially, so if they're too big for the sequence
            # dynamicaly resize them so they're 1/2 of the sequence length
            if extend_ends+1 >= len(B):
                extend_ends = len(B)/2
                
            if B[extend_ends+1] == 1:
                B[0:extend_ends] = [1]*len(B[0:extend_ends])

            if B[:-(extend_ends+1)] == 1:
                B[-extend_ends:] = [1]*len(B[0:extend_ends])
            
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

        new_domains[protein.unique_ID] = local_domains

    return new_domains


            
                    
                
                
            
