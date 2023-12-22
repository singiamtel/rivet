#! /usr/bin/env bash

if [[ -z "$1" ]]; then
    echo "Usage: `basename $0` <path/to/some_analysis_S123456.cc>"
    exit 1
fi

CCPATH=$1
echo "Processing $CCPATH"

plugindir=$(dirname $CCPATH)
aname=$(basename $CCPATH .cc)
scode=$(echo $aname | sed -e 's/.*_S\([^_]\+\).*/\1/')
if [[ -z "$scode" ]]; then
    echo "  No SPIRES code identified in $aname"
    exit 1
fi
echo "  SPIRES code = $scode"

cpath=$plugindir/$aname.cc
ipath=$plugindir/$aname.info
ppath=$plugindir/$aname.plot
ypath=$plugindir/$aname.yoda.gz

icode=$(grep "InspireID" $ipath | sed -e 's/InspireID:\s\(.\+\)/\1/')
echo "  Inspire code = $icode"

newaname=$(echo $aname | sed -e s/S$scode/I$icode/)
newcpath=$(echo $cpath | sed -e s/S$scode/I$icode/)
newipath=$(echo $ipath | sed -e s/S$scode/I$icode/)
newypath=$(echo $ypath | sed -e s/S$scode/I$icode/)
newppath=$(echo $ppath | sed -e s/S$scode/I$icode/)
echo "  Making $plugindir/$newaname.{.cc,.info,.plot,.yoda}"

cat $cpath | sed -e s/"RIVET_DECLARE_ALIASED_PLUGIN(\($aname\), \($newaname\));/RIVET_DECLARE_ALIASED_PLUGIN(XXXXX, YYYYY);/" > $newcpath.tmp
cat $newcpath.tmp | sed -e s/S$scode/I$icode/g -e s/XXXXX/$newaname/g -e s/YYYYY/$aname/g > $newcpath
rm -f $newcpath.tmp

cat $ipath | sed -e s/$aname/$newaname/g > $newipath

cat $ppath | sed -e s/$aname/$newaname/g > $newppath

cp $ypath $newypath
gunzip $newypath
cat ${newypath%.gz} | sed -e s/$aname/$newaname/g > ${newypath%.gz}.tmp
mv ${newypath%.gz}.tmp ${newypath%.gz}
(cd $(dirname $newypath) && gzip $(basename ${newypath%.gz} | sed -e s/yoda.gz/yoda/))
rm -f ${newypath%.gz}

git rm $plugindir/$aname.*
