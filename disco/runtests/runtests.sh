#!/bin/bash


for testname in test_*; do
	./runtest.sh $testname || {
		echo "Error: test $testname failed";
		exit 1;
	}
done

echo "All tests OK"
