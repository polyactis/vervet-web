#!/bin/bash
# 2011-10-05 run stuff after nucmer
if test $# -lt 6
then
	echo "Usage: $0 deltaFname coordsFname filterFname refFname qryFname plotPrefix"
	echo
	echo "Note:"
	echo "	A script does show-coords, delta-filter, mummerplot after nucmer/mummer is run."
	echo
	echo "Example:"
	echo "	$0 ref_qry.delta ref_qry.coords ref_qry.filter ref.fasta qry.fasta ref_qry"
exit
fi

deltaFname=$1
coordsFname=$2
filterFname=$3
refFname=$4
qryFname=$5
plotPrefix=$6

mummerBinaryPath=~/bin/MUMmer/

#-c	Include percent coverage columns in the output
#-r	Sort output lines by reference
#-l	Include sequence length columns in the output
$mummerBinaryPath/show-coords -r -l -c $deltaFname > $coordsFname
exitCode=$?
if test "$exitCode" != "0"
then
	exit $exitCode
fi

#-q	Query alignment using length*identity weighted LIS. For each query, leave only the alignments which form the longest consistent set for the query
#-r	Reference alignment using length*identity weighted LIS. For each reference, leave only the alignments which form the longest consistent set for the reference.
#-l int	Set the minimum alignment length (default 0)
$mummerBinaryPath/delta-filter -l 1000 -r $deltaFname > $filterFname
exitCode=$?
if test "$exitCode" != "0"
then
	exit $exitCode
fi

#-t string	Set the output terminal to x11, postscript or png
#--terminal	--x11 --postscript --png

numberOfLines=`wc $filterFname -l|awk -F ' ' '{print $1}'`
if test $numberOfLines -gt 3
then
	$mummerBinaryPath/mummerplot $filterFname -R $refFname -Q $qryFname --prefix $plotPrefix -t png
	gnuplotScript=$plotPrefix.gp
	#2012.8.15 correct gnuplot script bugs
	sed 's/set terminal png tiny/set terminal png /' $gnuplotScript > $gnuplotScript.tmp
	mv $gnuplotScript.tmp $gnuplotScript
	sed 's/set ticscale 0 0/set tics scale 0, 0/' $gnuplotScript > $gnuplotScript.tmp
	mv $gnuplotScript.tmp $gnuplotScript
	gnuplot $gnuplotScript
else
	echo "no synteny after filter."
fi
