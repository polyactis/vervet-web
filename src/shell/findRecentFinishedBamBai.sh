#!/bin/bash
#2011-11-6 Yu Huang
#Two criteria to find recently-finished non-corrupt bam.
# 1. the bai was last modified >60 minutes. (some MarkDuplicates crashed and had a partial bam and no bai.)
# 2. bai's bam is present and readable.



for i in `find  /u/home/eeskin/polyacti/NetworkData/vervet/db/individual_alignment/ -name *_vs_524_by_2_dupMark*bam.bai -mmin +60|sort`; do
	prefix=`echo $i|awk -F . '{print $1}'`;
	bamFname=$prefix.bam;
	baiFname=$i;
	if test -r $bamFname; then
		echo $bamFname;
	fi
done
