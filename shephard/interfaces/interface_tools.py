from shephard.exceptions import InterfaceException


def check_proteome(p, function_name):
    if "<class 'shephard.proteome.Proteome'>" == str(p.__class__):
        return None
    else:
        raise InterfaceException('First argument passed to function [%s] was not a proteome' %(function_name))
