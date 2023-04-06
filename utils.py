import requests, sys
import vcf

import pandas as pd


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
    assert list_delim in d[key], "Key {} does not contain delimiter {}".format(
        key, list_delim
    )
    vals = d[key].split(list_delim)  # list of split values
    d[key] = vals[0]  # original key gets first value
    for i, val in enumerate(vals[1:]):
        # add key1, key2, etc. and corresponding values
        d["{}{}".format(key, i + 1)] = val


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
    indent : `str`, optional (default="")
        For pretty traceback. Ignored if `verbose==False`.

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
                    print(
                        "{}Joining key {} using delimiter {}".format(
                            indent, key, list_delim
                        )
                    )
            elif len(d_unpacked[key]) == 1:
                # if value is list with one element, unpack
                d_unpacked[key] = d_unpacked[key][0]
                if verbose:
                    print(
                        "{}Unpacked 1 value from key {}: {}".format(
                            indent, key, d_unpacked[key]
                        )
                    )
            else:
                # if empty list (for some reason), remove the key
                del d_unpacked[key]
        elif isinstance(d[key], dict):
            # if key is another dict, do full unpacking
            if verbose:
                print(
                    "{}Unpacking key {} as dictionary with {} keys".format(
                        indent, key, len(d_unpacked[key])
                    )
                )
            # generate new keys from dict
            unpacked_key = unpack_dict(
                d_unpacked[key],
                list_delim=list_delim,
                verbose=verbose,
                indent="{}\t".format(indent),
            )
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
    assert isinstance(
        call_obj, vcf.model._Call
    ), "call_obj must be a vcf.model._Call object"
    call_dict = {}  # initialize empty dictionary
    keys = format_str.split(format_delim)  # list of expected keys in call_obj
    for key in keys:
        call_dict[key] = getattr(call_obj.data, key)

    return call_dict


def extract_alt(d, alt_obj, alt_key="ALT"):
    """
    Given a dictionary containing `vcf.model._Substitution` object, unpack
    values into dictionary

    Parameters
    ----------
    d : `dict`
        Dictionary with "ALT" key to unpack
    alt_obj : `vcf.model._Substitution` or `list`
        PyVCF `vcf.model._Substitution` object to extract information from
        (i.e. `record.ALT`)
    alt_key : `str`, optional (default="ALT")
        key within `d` containing PyVCF `model._Substitution` object

    Returns
    -------
    d_unpacked : `dict`
        Dictionary with unpacked values
    """
    d_unpacked = d.copy()

    if isinstance(alt_obj, vcf.model._Substitution):
        # get variant type from alt_key
        d_unpacked["{}_type".format(alt_key)] = alt_obj.type
        # get variant sequence as string and replace alt_key
        d_unpacked[alt_key] = alt_obj.sequence
    elif isinstance(alt_obj, list):
        # unpack list of ALTs
        d_unpacked["{}_type".format(alt_key)] = ";".join([x.type for x in alt_obj])
        d_unpacked[alt_key] = ";".join([x.sequence for x in alt_obj])
    else:
        print("alt_obj is not a list or vcf.model._Substitution object; skipping.")

    return d_unpacked


def vep_region_post(
    record,
    server="https://grch37.rest.ensembl.org",
    ext="/vep/homo_sapiens/region",
    **kwargs,
):
    """
    Given a `vcf.model._Record` object, send POST request to Ensembl
    vep/:species/region

    Parameters
    ----------
    record : `vcf.model._Record`
        PyVCF `model._Record` object to send to Ensembl API
    server : `str`, optional (default="https://grch37.rest.ensembl.org")
        URL to Ensembl server
    ext : `str`, optional (default="/vep/homo_sapiens/region")
        URL extension for region annotation endpoint
    **kwargs
        Keyword arguments for VEP endpoint

    Returns
    -------
    vep_dict : `dict`
        Dictionary with returned values from Ensembl API call
    """
    assert isinstance(
        record, vcf.model._Record
    ), "record must be a vcf.model._Record object"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    r = requests.post(
        server + ext,
        headers=headers,
        data='{{ "variants" : ["{} {} {} {} {} . . ."], {} }}'.format(
            record.CHROM,
            record.POS,
            record.ID if record.ID is not None else ".",
            record.REF,
            record.ALT[0].sequence,
            str(kwargs).replace("{", "").replace("}", "").replace("'", '"'),
        ),
    )

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    vep_dict = r.json()
    assert len(vep_dict) == 1, "More than one VEP result returned!"
    return vep_dict[0]


def overlap_get(record, server="https://grch37.rest.ensembl.org"):
    """
    Given a `vcf.model._Record` object, send POST request to Ensembl
    vep/:species/region

    Parameters
    ----------
    record : `vcf.model._Record`
        PyVCF `model._Record` object to send to Ensembl API
    server : `str`, optional (default="https://grch37.rest.ensembl.org")
        URL to Ensembl server

    Returns
    -------
    overlap_dict : `dict`
        Dictionary with returned values from Ensembl API call
    """
    assert isinstance(
        record, vcf.model._Record
    ), "record must be a vcf.model._Record object"
    ext = "/overlap/region/homo_sapiens/{}:{}-{}?feature=gene".format(
        record.CHROM,
        record.start,
        record.end,
    )
    headers = {"Content-Type": "application/json"}
    r = requests.get(server + ext, headers=headers)

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    overlap_dict = r.json()
    return overlap_dict


def process_record(record, index=1):
    """
    Given a `vcf.model._Record` object, annotate it and return as a dataframe

    Parameters
    ----------
    record : `vcf.model._Record`
        PyVCF `model._Record` object to annotate
    index : `int`, optional (default=1)
        Index of resulting `pd.DataFrame`

    Returns
    -------
    record_unpacked : `pd.DataFrame`
        Annotated dataframe row corresponding to `record`
    """
    # unpack record into dictionary of depth 1
    record_unpacked = unpack_dict(record.__dict__)
    # unpack ALT allele objects
    record_unpacked = extract_alt(record_unpacked, alt_obj=record.ALT, alt_key="ALT")
    # unpack model._Call object
    call_dict = unpack_dict(
        extract_call(record_unpacked["samples"], format_str=record_unpacked["FORMAT"]),
    )
    # remove old keys
    del record_unpacked["samples"]
    del record_unpacked["FORMAT"]
    # merge dicts
    record_unpacked = {**record_unpacked, **call_dict}

    # get gene where mutation is
    overlap_dict = overlap_get(record)
    try:
        # first record
        record_unpacked["gene_id"] = overlap_dict[0]["gene_id"]
        record_unpacked["gene_symbol"] = overlap_dict[0]["external_name"]
        record_unpacked["gene_biotype"] = overlap_dict[0]["biotype"]
        record_unpacked["gene_description"] = overlap_dict[0]["description"]
        if len(overlap_dict) > 1:
            for i in range(1, len(overlap_dict)):
                record_unpacked["gene_id{}".format(i)] = overlap_dict[i]["gene_id"]
                record_unpacked["gene_symbol{}".format(i)] = overlap_dict[i][
                    "external_name"
                ]
                record_unpacked["gene_biotype{}".format(i)] = overlap_dict[i]["biotype"]
                record_unpacked["gene_description{}".format(i)] = overlap_dict[i][
                    "description"
                ]
    except:
        print("No gene overlap found for {}".format(record))

    # call Ensembl VEP endpoint
    vep_dict = vep_region_post(
        record,
        hgvs=1,
        # protein=1,
        # uniprot=1,
        # GO=1,
        # vcf_string=1,
        # domains=1,
        distance=0,  # don't look for upstream/downstream effects for simplicity
    )
    # extract VEP info
    record_unpacked["most_severe_consequence"] = vep_dict["most_severe_consequence"]

    # get MAF and snp info if available
    if "colocated_variants" in vep_dict:
        for d in vep_dict["colocated_variants"]:
            if "minor_allele" in d:
                record_unpacked["minor_allele"] = d["minor_allele"]
                record_unpacked["minor_allele_freq"] = d["minor_allele_freq"]
                record_unpacked["snp_id"] = d["id"]

    return pd.DataFrame(record_unpacked, index=[index])


def calc_frequencies(out_df):
    """
    Calculate variant and reference allele frequencies (VAF and RAF)

    Parameters
    ----------
    out_df : `pd.DataFrame`
        DataFrame containing annotations and columns `["NR","NV"]`

    Returns
    -------
    `out_df` is edited in place, adding "NRef", "VAF" and "RAF" columns
    """
    out_df[["NR", "NR1"]] = out_df["NR"].astype(str).str.split(";", expand=True)
    out_df[["NV", "NV1"]] = out_df["NV"].astype(str).str.split(";", expand=True)
    out_df["VAF"] = out_df["NV"].astype(int) / out_df["NR"].astype(int)
    out_df["NRef"] = out_df["NR"].astype(int) - out_df["NV"].astype(int)
    out_df["RAF"] = out_df["NRef"] / out_df["NR"].astype(int)

    # now in rows with multiple values delimited by ;
    out_df.loc[~out_df.NR1.isnull(), "VAF1"] = out_df.loc[
        ~out_df.NR1.isnull(), "NV1"
    ].astype(int) / out_df.loc[~out_df.NR1.isnull(), "NR1"].astype(int)
    out_df.loc[~out_df.NR1.isnull(), "NRef1"] = out_df.loc[
        ~out_df.NR1.isnull(), "NR1"
    ].astype(int) - out_df.loc[~out_df.NR1.isnull(), "NV1"].astype(int)
    out_df.loc[~out_df.NR1.isnull(), "RAF1"] = out_df.loc[
        ~out_df.NR1.isnull(), "NRef1"
    ].astype(int) / out_df.loc[~out_df.NR1.isnull(), "NR1"].astype(int)

    # recombine columns
    for col in ["VAF", "NRef", "RAF", "NR", "NV"]:
        out_df.loc[~out_df.NR1.isnull(), col] = (
            out_df.loc[~out_df.NR1.isnull(), col].astype(str)
            + ";"
            + out_df.loc[~out_df.NR1.isnull(), "{}1".format(col)].astype(str)
        )

    out_df.drop(columns=["VAF1", "NRef1", "RAF1", "NR1", "NV1"], inplace=True)
