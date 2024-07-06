#! /usr/bin/env bash

echo PYTHONPATH=$PYTHONPATH
echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH
echo PATH=$PATH
echo

./testMerging

echo "Testing serialisation diff ..."
yodadiff -q ref-merge-serial.yoda out-merge-serial.yoda
echo "OK"

echo "Testing equivalent merging diff ..."
yodadiff -q ref-merge-equiv.yoda out-merge-equiv.yoda
echo "OK"

echo "Testing non-equivalent merging diff ..."
yodadiff -q ref-merge-nonequiv.yoda out-merge-nonequiv.yoda
echo "OK"

echo "Testing AH-merging merging diff ..."
yodadiff -q ref-merge-ahmerge.yoda out-merge-ahmerge.yoda
echo "OK"

