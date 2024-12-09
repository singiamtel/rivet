#! /usr/bin/env python3

import argparse
ap = argparse.ArgumentParser()
ap.add_argument("AUTHORS", help="file of authors and email addresses from the other script")
args = ap.parse_args()

def surname(namestr):
    return namestr.split("~")[-1]

authstrs = []
import re
with open(args.AUTHORS) as af:
    for line in af.readlines():
        l = line.strip()
        m = re.match(r"^(.+) <(.*)>", l)
        if not m:
            print("ERROR: " + l)
            continue
        name = m.group(1)
        nameparts = name.split()
        for i in range(len(nameparts)-1):
            nameparts[i] = (nameparts[i][0] + ".")
        namestr = "~".join(nameparts)
        authstrs.append(namestr)

print(",\n".join(sorted(authstrs, key=surname)))
