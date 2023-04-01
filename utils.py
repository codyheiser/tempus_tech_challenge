import vcf


def split_vals(d, key, list_delim=";"):
    """
    Given a dictionary and key within it, split values into multiple keys

    Parameters
    ----------
    d : `dict`
        Dictionary
    key : `str`
        Key within `d` to split
    list_delim : `str`, optional (default=";")
        String to delimit list values in `d[key]` with

    Returns
    -------
    d : `dict`
        Dictionary with split values
    """
    assert list_delim in d[key], "Key {} does not contain delimiter {}".format(key, list_delim)
    vals = d[key].split(list_delim)  # list of split values
    d[key] = vals[0]  # original key gets first value
    for i, val in enumerate(vals[1:]):
        # add key1, key2, etc. and corresponding values
        d["{}{}".format(val, i+1)] = val


def unpack_dict(d, list_delim=";", verbose=False, indent=""):
    """
    Given a dictionary, unpack values such that lists become singular.
    Intended to prep dictionaries for pandas compatibility.

    Parameters
    ----------
    d : `dict`
        Dictionary to unpack
    list_delim : `str`, optional (default=";")
        String to delimit list values in `d` with
    verbose : `bool`, optional (default=`False`)
        If `True`, print information when values are unpacked

    Returns
    -------
    d_unpacked : `dict`
        Dictionary with unpacked values
    """
    d_unpacked = d.copy()
    for key in d:
        if isinstance(d[key], list):
            # if key is a list, do some unpacking
            if len(d[key]) > 1:
                # if value is list, join by list_delim
                d_unpacked[key] = [str(elem) for elem in d_unpacked[key]]
                d_unpacked[key] = list_delim.join(d_unpacked[key])
                if verbose:
                    print("{}Joining key {} using delimiter {}".format(indent, key, list_delim))
            elif len(d_unpacked[key]) == 1:
                # if value is list with one element, unpack
                d_unpacked[key] = d_unpacked[key][0]
                if verbose:
                    print("{}Unpacked 1 value from key {}: {}".format(indent, key, d_unpacked[key]))
            else:
                # if empty list (for some reason), remove the key
                del d_unpacked[key]
        elif isinstance(d[key], dict):
            # if key is another dict, do full unpacking
            if verbose:
                print("{}Unpacking key {} as dictionary with {} keys".format(indent, key, len(d_unpacked[key])))
            # generate new keys from dict
            unpacked_key = unpack_dict(d_unpacked[key], list_delim=list_delim, verbose=verbose, indent="{}\t".format(indent))
            # delete original key-value pair
            del d_unpacked[key]
            # merge the two dicts
            d_unpacked = {**d_unpacked, **unpacked_key}
            if verbose:
                print("{}\tAdded {} keys to dict".format(indent, len(unpacked_key)))
        elif d[key] is None:
            # if empty value, just remove it
            del d_unpacked[key]

    return d_unpacked


def extract_call(call_obj, format_str="GT:GL:GOF:GQ:NR:NV", format_delim=":"):
    """
    Given a `vcf.model._Call` object, unpack values into dictionary

    Parameters
    ----------
    call_obj : `vcf.model._Call`
        PyVCF `model._Call` object to extract information from
    format_str : `str`, optional (default="GT:GL:GOF:GQ:NR:NV")
        String describing format of `vcf.model._Call` object, and keys present
    format_delim : `str`, optional (default=":")
        String to delimit list values in `format_str` with

    Returns
    -------
    call_dict : `dict`
        Dictionary with unpacked values from `call_obj`
    """
    assert isinstance(call_obj, vcf.model._Call), "call_obj must be a vcf.model._Call object"
    call_dict = {}  # initialize empty dictionary
    keys = format_str.split(format_delim)  # list of expected keys in call_obj
    for key in keys:
        call_dict[key] = getattr(call_obj.data, key)

    return call_dict

