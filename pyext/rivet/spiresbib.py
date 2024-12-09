#! /usr/bin/env python3

# TODO: Separate this library from the executable script (currently via shebang above and `if __name__` below.

import logging, re
from urllib.request import urlopen

usage = """%prog <spiresid> [<spiresid2> ...]

Given Inspire and SPIRES paper IDs, fetch the corresponding BibTeX db entry from
the SPIRES Web interface and write it to stdout. Prefix the code with I or S
appropriately.

TODO:
 - remove SPIRES option and handling code
"""

def fetch_bibtex(iscode, refid):
    "Download and decode BibTeX"

    logging.debug("Fetching BibTeX")
    if iscode.upper() == "I":
        url = f"https://inspirehep.net/api/literature/{refid}?format=bibtex"
        logging.debug("Downloading Inspire BibTeX from %s" % url)
    else:
        raise Exception("SPIRES lookup no longer supported")

    hresp = urlopen(url)
    bibtex_bytes = hresp.read()
    encoding = hresp.headers.get_content_charset('utf-8')
    bibtex = bibtex_bytes.decode(encoding)
    hresp.close()

    return bibtex


def extract_bibtex(bibtxt):
    "Extract BibTeX key from block"

    logging.debug("Extracting BibTeX detail")
    bib = bibtxt.strip()

    ## Get BibTeX key
    re_bibtexkey = re.compile(r'^@.+?{(.+?),$', re.MULTILINE)
    m = re_bibtexkey.search(bib)
    if m is None:
        return None, bib
    key = m.group(1)

    ## Return key and BibTeX
    return key, bib


def get_bibtex_from_repo(iscode, refid):
    html = fetch_bibtex(iscode, refid)
    key, bibtex = extract_bibtex(html)
    return key, bibtex


def get_bibtexs_from_repos(iscodes_refids):
    bibdb = {}
    for iscode, refid in iscodes_refids:
        key, bibtex = get_bibtex_from_repo(iscode, refid)
        if key and bibtex:
            bibdb[refid] = (key, bibtex)
    return bibdb


if __name__ == '__main__':
    ## Parse command line options
    import argparse
    ap = argparse.ArgumentParser(usage=usage)
    ap.add_argument("CODES", nargs="*")
    args = parser.parse_args()

    ## Make individual bibinfo files
    for code in args.CODES:
        iscode = code[0]
        refid = code[1:]
        key, bibtex = get_bibtex_from_repo(iscode, refid)
        import sys
        f = sys.stdout
        f.write("BibKey: %s\n" % key)
        f.write("BibTeX: '%s'\n" % bibtex)

    # ## Build ref db
    # bibdb = get_bibtexs_from_spires(args.CODES)
    # for sid, (key, bibtex) in bibdb.iteritems():
    #     print key, "=>\n", bibtex

    # ## Pickle ref db9151176
    # import cPickle as pickle
    # fpkl = open("spiresbib.pkl", "w")
    # pickle.dump(bibdb)
    # fpkl.close()
