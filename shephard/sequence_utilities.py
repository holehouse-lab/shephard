

def inside_region(boundary_start, boundary_end, position):
    """
    Given a position, does the position fall inside the defined region 
    """
    if position >= boundary_start and position <= boundary_end:
        return True
    else:
        return False



def update_from_realworld(vals, use_realworld_numbering=True):
    """
    Takes in a list of values and IF use_realworld_numbering is
    true will subtract 1 from each (i.e. moving from a 1-indexed
    numbering to a 0-indexed numbering).

    If use_realword is false then the same numbers are returned
    and no alterations are made.

    """

    if use_realworld_numbering:
        return_vals = []
        for i in vals:
            return_vals.append(i - 1)
        return return_vals

    else:
        return vals


def update_to_realworld(vals, use_realworld_numbering=True):
    """
    Takes in a list of values and IF use_realworld_numbering is
    true will add 1 from each (i.e. moving from a 0-indexed
    numbering to a 1-indexed numbering).

    If use_realword is false then the same numbers are returned
    and no alterations are made.

    """

    if use_realworld_numbering:
        return_vals = []
        for i in vals:
            return_vals.append(i + 1)
        return tuple(return_vals)

    else:
        return tuple(vals)

            
