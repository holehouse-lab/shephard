def cast_attributes(obj, include=[], exclude=[], cast_type=float, skip_failed=False):
    """
    Function which takes a SHEPHARD object with attributes and then casts 
    key associated with the attributes to the type defined by cast_type.

    Note this works by reading the value associated with the attribute
    key, casting to the associated cast type, and then re-defining the
    attribute key with the new cast value.

    Note you should define EITHER include or exclude, not both. Which 
    you choose will depend on the sitiation. 

    The default behavior is to cast all attributes to the defined type.

    Parameters
    ----------
    obj : SHEPHARD object
        Object with attributes to cast. Must have an .attributes variable

    include : list
        List of attribute keys to explicitly include in the casting. Keys
        should be strings.

    exclude : list
        List of attribute keys to explicitly exclude in the casting. Keys
        should be strings.

    cast_type : type
        Type to cast the attribute values to. Default is float.

    skip_failed : bool
        If True, we will try and cast attributes and if they throw and 
        exception we'll just move on. If False, we'll raise an exception
        if we try and cast an attribute and it cannot be cast.

    Returns
    -------
    None
        No return type, but the object will have its attributes cast to 
        the defined type.

    Raises
    ------
    ValueError
        If both include and exclude are defined, if the object does not
        have attributes, or if the cast_type is not a valid type.

    """

    # validate only one of include or exclude is defined
    if len(include) > 0 and len(exclude) > 0:
        raise ValueError('Cannot have both include and exclude keys')

    # validate cast_type is a type
    if not isinstance(cast_type, type):
        raise ValueError('The cast_type must be a valid type')

    # confirm object has attributes
    if not hasattr(obj, 'attributes'):
        raise ValueError('Passed object does not have attributes')

    # if neither include or exclude are defined, means we cast
    # everything
    if len(include) == 0 and len(exclude) == 0:
        exclude_override = True
    else:
        exclude_override = False
        
    # if we had attributes keys to explicitly include 
    if len(include) > 0:
        for a in include:
            if a in obj.attributes:
                try:
                    # we are just overwriting the attribute in the obj after casting it
                    obj.add_attribute(a, cast_type(obj.attribute(a)), safe=False)
                    
                except ValueError as e:
                    if skip_failed:
                        continue
                    else:
                        print('Failed to cast attribute: {}'.format(a))
                        raise e

    # if we had attributes keys to explicitly exclude (note the default behavior
    # is to cast everything)
    if len(exclude) > 0 or exclude_override:
        for a in obj.attributes:
            if a not in exclude:
                try:
                    # we are just overwriting the attribute in the obj after casting it
                    obj.add_attribute(a, cast_type(obj.attribute(a)), safe=False)
                    
                except ValueError as e:
                    if skip_failed:
                        continue
                    else:
                        print('Failed to cast attribute: {}'.format(a))
                        raise e
                    
