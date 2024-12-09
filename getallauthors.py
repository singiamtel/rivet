#! /usr/bin/env python3

import argparse
ap = argparse.ArgumentParser()
ap.add_argument("--blacklist", dest="BLACKLIST", help="file of analyses to exclude from address extraction")
args = ap.parse_args()

ana_blacklist = []
if args.BLACKLIST:
    with open(args.BLACKLIST) as blf:
        ana_blacklist = [l.strip() for l in blf.readlines()]
    #print(ana_blacklist)

import rivet, re
authors_emails = {}
for aname in rivet.AnalysisLoader.analysisNames():
    if aname in ana_blacklist:
        continue
    ana = rivet.AnalysisLoader.getAnalysis(aname)
    for au_em in ana.authors():
        au, em = None, None
        if "<" not in au_em:
            if "@" in au_em:
                au = au_em
                em = au_em
            else:
                au = au_em
        else:
            m = re.search("(.*)<(.*)>.*", au_em)
            if m:
                au, em = m.group(1).strip(), m.group(2).strip()
        if au or em:
            if em or au not in authors_emails:
                authors_emails[au] = em

for au, em in sorted(authors_emails.items()):
    print(r"{} <{}>".format(au, em))
    #.encode("utf-8"))
