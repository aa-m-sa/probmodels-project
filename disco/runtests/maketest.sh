#!/bin/bash

name="$1"

testspath=$(pwd)
mainpath=$(dirname $testspath)

testpath=$testspath/"$name"
outputdir=$testpath/output
rm -rf $outputdir
mkdir $outputdir

cd $outputdir
while read cmd || [ -n "$cmd" ]
do
	cmd="${cmd//%INPUT%/$testpath/input}"
	$mainpath/$cmd >stdout.txt
	#status=$?
	#if [ $status -ne 0 ]; then
	#	echo "invalid return status" >&2
	#	exit 1
	#fi
done < $testspath/"$name"/commands

cd $testspath
