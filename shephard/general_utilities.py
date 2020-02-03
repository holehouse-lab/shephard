from . import exceptions

def cast_or_none(value, cast_type):
    """
    Functo

    """
    if value is None:
        return None
    else:
        return cast_type(value)


def string_to_list_of_strings(inval):
    """
    Function that converts a string or list of
    strings to a list of strings. 

    """

    if type(inval) == str:
        return [inval]
    elif type(inval) == list:
        return inval

    else:
        raise exceptions.ShephardException('Expected either a single string or a list of strings')


def numerical_average(l):
    return numerical_sum(l)/len(l)


def numerical_sum(l):

    t = 0

    for i in l:
        t = t + i

    return t


