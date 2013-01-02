#!/bin/bash
# 2012.3.1 
if test $# -lt 2
then
	echo
	echo "Usage: $0 output input1 input2 ..."
	echo
	echo "2nd Usage: $0 -o output input1 input2 ..."
	echo
	echo "Note:"
	echo "	#. It invokes unix command cat to do the actual job."
	echo "	#. To maintain interface compatibility with other python reducers, 2nd-usage is also supported. "
	echo
	echo "Example:"
	echo "	$0 Contig0.tsv Contig0_1_2000000.tsv Contig0_2000000_3000000.tsv"
exit
fi

catPath=/bin/cat

output=$1
if test "$output" = "-o"
then
	shift
	output=$1
fi
shift
no_of_inputs=$#
echo -n "Merge $no_of_inputs input files into $output ... "
$catPath $* >$output
echo ""
