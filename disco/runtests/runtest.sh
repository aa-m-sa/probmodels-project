#!/bin/bash

name="$1"

#echo -n "Running '${name}'... "
printf "%-50s" "Running '${name}'... "

# cmds=$(cat "$name"/commands)
# for cmd in "$cmds"; do
# 	echo "komento: " $cmd
# done

testspath=$(pwd)
mainpath=$(dirname $testspath)
testpath=$testspath/"$name"
outputdir=$testpath/output

tmpdir=$testspath/$(mktemp $name.tmp.XXXX -d)

cd $tmpdir
while read cmd || [ -n "$cmd" ]
do
	cmd="${cmd//%INPUT%/$testpath/input}"
	$mainpath/$cmd >stdout.txt 2>/dev/null
	#status=$?
	#if [ $status -ne 0 ]; then
	#	echo "invalid return status" >&2
	#	exit 1
	#fi
done < $testspath/"$name"/commands

cmp -s <(ls $outputdir) <(ls $tmpdir) || {
	echo "Error: different output filenames";
	#diff <(ls $outputdir) <(ls $tmpdir);
	echo "**** ORIGINAL: *****"
	cat <(ls $outputdir)
	echo "**** CURRENT: ******"
	cat <(ls $tmpdir);
	echo "********************"
	exit 1;
}


for file in *; do
	cmp -s $outputdir/$file $tmpdir/$file || {
		echo "Error: different content in file $file";
		exit 1;
	}
done

cd $testspath

rm -r $tmpdir

echo "OK"
