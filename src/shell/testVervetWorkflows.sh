#!/bin/bash

source ~/.bash_profile

if test $# -lt 2
then
        echo "Usage: $0 computingSiteID dbPassword"
        echo
        echo "This script would gemerate dag files for a few workflows and submit them. Used to test whether everything work together. Need to setup Pegasus workflow email notification."
        echo
        echo "Example:"
        echo "  $0 hcondor dbPassword"
exit
fi

siteID=$1
dbPassword=$2
#changing the directory
pushd ~/NetworkData/vervet/workflow

echo "Executing $0 $*."

date
echo "Testing TestMapReduceGenomeFileWorkflow.py (calculating LiftOver locus mapping probability) ... "
dagFname=dags/LiftPolymorphismCoordinates/Method109_To_225_LiftOverProbability.test.xml
~/script/pymodule/pegasus/TestMapReduceGenomeFileWorkflow.py --inputDir ./LiftPolymorphismCoordinates/FindNewRefCoordinates_Method109_vs_3488_BWA_F49.2013.Jul.19T141746/folderReduceGzip/ -H -C 40 -j $siteID -l $siteID -D ~/NetworkData/vervet/db/ -t ~/NetworkData/vervet/db/ -o $dagFname --inputSuffixList .tsv --db_user yh -z localhost --ref_ind_seq_id 3280 --ref_genome_tax_id 60711 --ref_genome_sequence_type_id 9 --ref_genome_version 2 --db_passwd $dbPassword
exitCode=$?
if test $exitCode -eq 0; then
	./submit $dagFname $siteID
fi

date
echo "Testing PlinkOnVCFWorkflow.py run-type 5, marking mendel errors missing ... "
#otherwise the dag-generating script would pause to avoid file overwriting
rm ./aux/Method225_MarkMendelErrorCallMissing_merge_list.test.txt
dagFname=dags/MarkMendelErrorCallMissing/MarkMendelErrorCallMissing_method225.test.xml
~/script/vervet//src/PlinkOnVCFWorkflow.py -I ~/NetworkData/vervet/db/genotype_file/method_225/ -o $dagFname -y5 --clusters_size 1 --needSSHDBTunnel --site_handler $siteID --input_site_handler $siteID --db_user yh --hostname localhost --local_data_dir ~/NetworkData/vervet/db/ --data_dir ~/NetworkData/vervet/db/ --maxContigID 100 --ref_ind_seq_id 3488 --mergeListFname ./aux/Method225_MarkMendelErrorCallMissing_merge_list.test.txt --db_passwd $dbPassword
exitCode=$?
if test $exitCode -eq 0; then
	./submit $dagFname $siteID
fi

date
echo "Testing FindNewRefCoordinatesGivenGenotypeMethodWorkflow.py (LiftOver from method 109 to ref sequence 3488 ..."
dagFname=dags/LiftPolymorphismCoordinates/FindNewRefCoordinates_Method109_vs_3488_BWA_F49.test.xml 
~/script/vervet/src/polymorphism/FindNewRefCoordinatesGivenGenotypeMethodWorkflow.py -I ~/NetworkData/vervet/db/genotype_file/method_109 --oldRefFastaFname ~/NetworkData/vervet/db/individual_sequence/3280_vervet_ref_6.0.3.fasta --newRefFastaFname ~/NetworkData/vervet/db/individual_sequence/3488_indID1_codeVRC_ref_sequencer3_seqType1_filtered0_version3.fasta --maxNoOfMismatches 2 -H -C 5 --no_of_aln_threads 1 -j $siteID -l $siteID -D ~/NetworkData/vervet/db/ -t ~/NetworkData/vervet/db/ -u yh -z localhost -o $dagFname --ref_ind_seq_id 3488 --alignmentMethodType 2 --intervalSize 40000 --flankingLength 49 --db_passwd $dbPassword
exitCode=$?
if test $exitCode -eq 0; then
	./submit $dagFname $siteID
fi


date
echo "Testing CalculateVCFStatPipeline.py (method id 247) ..."

refID=3499; country=135; pop=VRC; m=247; maxCID=2000; 
dagFname=dags/VCFStat/VCFStat_BeagleTrioCallerOnMethod$m.xml
~/script/vervet/src/popgen/CalculateVCFStatPipeline.py --ref_ind_seq_id $refID -I ~/NetworkData/vervet/db/genotype_file/method_$m/ -o $dagFname -j $siteID -l $siteID -u yh -z localhost -D ~/NetworkData/vervet/db/ -t ~/NetworkData/vervet/db/ --LDWindowSize 0 --minChrLengthForPlot 4000000 --minChrSize 1000000 -H --pop_country_id_ls $country --popHeader $pop --intervalSize 5000 --intervalOverlapSize 500 -x $maxCID -C 7 --samplingRate 0.01 --db_passwd $dbPassword
exitCode=$?
if test $exitCode -eq 0; then
	./submit $dagFname $siteID
fi


date
echo "Testing PlinkOnVCFWorkflow.py type 1, Mendel Error, (method id 247) ..."
m=247; refID=3499;
dagFname=dags/PlinkMendelError/PlinkMendelError_Method$m.xml
~/script/vervet//src/PlinkOnVCFWorkflow.py -I ~/NetworkData/vervet/db/genotype_file/method_$m/ -o $dagFname --ref_ind_seq_id $refID --checkEmptyVCFByReading --clusters_size 4 --needSSHDBTunnel --site_handler $siteID --input_site_handler $siteID --db_user yh --hostname localhost --local_data_dir ~/NetworkData/vervet/db/ --data_dir ~/NetworkData/vervet/db/ -y1 --maxContigID 3000 --db_passwd $dbPassword
exitCode=$?
if test $exitCode -eq 0; then
	./submit $dagFname $siteID
fi


#back to whatever initial folder
popd
