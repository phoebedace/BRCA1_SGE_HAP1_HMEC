#[CHANGE PATH NAMES BEFORE ATTEMPTING ANY COMMANDS]

#AFTER CREATING $out_dir, download the ASF sample sheet from the email they send when run is over, save it with a name specific to the run (see example below under asf_sample_sheet), and move it onto camp in the out_dir with the your file name, as example below

#Set up a project folder for data to be moved into in ~/home/shared/projects/SGE/BRCA1
#want files in seq_dir/fastq and seq_dir/collapsed (for Amplicon Sequencing service) to be copied over 
out_dir="/camp/home/dacep/home/shared/projects/SGE/BRCA1/241007_BRCA1i11_HMEC"

seq_dir="/camp/lab/findlayg/data/STPs/babs/inputs/phoebe.dace/asf/PM23140/20241004_LH00442_0057_A22NWWMLT3"

asf_sample_sheet="/camp/home/dacep/home/shared/projects/SGE/BRCA1/241007_BRCA1i11_HMEC/241007_PM23140_samplesheet.xlsx"


adapter1="CTGTCTCTTATACACATCTCCGAGCCCACGAGAC" #this is the nextera reverse adapter, as seen in R1
adapter2="CTGTCTCTTATACACATCTGACGCTGCCGACGA"  #this is the nextera forward adapter, as seen in R2
seq_type="L" # N for nextseq, M for Miseq, L for novaseq
ref_folder="/camp/home/dacep/home/shared/projects/SGE/BRCA1/fasta"

amplicon_list="BRCA1i11" #SGE experiment, splits on commas

#6 entries per experimental grouping must be specified, order must be:  experiment, pre, post, HDRLlib, neg, rna, as in here:
#splits on commas, then + : NOTE for BRCA1, using D5 in place of RNA to prevent errors

cigar_comparisons="BRCA1i11rHMEC1_D7+BRCA1i11rHMEC1_D14,BRCA1i11rHMEC2_D7+BRCA1i11rHMEC2_D14,BRCA1i11rHMEC1_D7+BRCA1i11rHMEC1_D14o,BRCA1i11rHMEC2_D7+BRCA1i11rHMEC2_D14o,BRCA1i11rHMEC1_D7+BRCA1i11rHMEC1_D21,BRCA1i11rHMEC2_D7+BRCA1i11rHMEC2_D21,BRCA1i11rHMEC1_D7+BRCA1i11rHMEC1_D21o,BRCA1i11rHMEC2_D7+BRCA1i11rHMEC2_D21o"

exp_groupings="BRCA1i11rHMEC1D7D14+BRCA1i11rHMEC1_D7+BRCA1i11rHMEC1_D14+BRCA1i11_HDRL+BRCA1i11_neg+BRCA1i11rHMEC1_D7,BRCA1i11rHMEC2D7D14+BRCA1i11rHMEC2_D7+BRCA1i11rHMEC2_D14+BRCA1i11_HDRL+BRCA1i11_neg+BRCA1i11rHMEC2_D7,BRCA1i11rHMEC1D7D14o+BRCA1i11rHMEC1_D7+BRCA1i11rHMEC1_D14o+BRCA1i11_HDRL+BRCA1i11_neg+BRCA1i11rHMEC1_D7,BRCA1i11rHMEC2D7D14o+BRCA1i11rHMEC2_D7+BRCA1i11rHMEC2_D14o+BRCA1i11_HDRL+BRCA1i11_neg+BRCA1i11rHMEC2_D7,BRCA1i11rHMEC1D7D21+BRCA1i11rHMEC1_D7+BRCA1i11rHMEC1_D21+BRCA1i11_HDRL+BRCA1i11_neg+BRCA1i11rHMEC1_D7,BRCA1i11rHMEC2D7D21+BRCA1i11rHMEC2_D7+BRCA1i11rHMEC2_D21+BRCA1i11_HDRL+BRCA1i11_neg+BRCA1i11rHMEC2_D7,BRCA1i11rHMEC1D7D21o+BRCA1i11rHMEC1_D7+BRCA1i11rHMEC1_D21o+BRCA1i11_HDRL+BRCA1i11_neg+BRCA1i11rHMEC1_D7,BRCA1i11rHMEC2D7D21o+BRCA1i11rHMEC2_D7+BRCA1i11rHMEC2_D21o+BRCA1i11_HDRL+BRCA1i11_neg+BRCA1i11rHMEC2_D7"


#A file listing the expected pam edits, the sge region in the amplicon, and which gRNA was used, etc.
editing_info_file="/camp/home/dacep/home/shared/projects/SGE/BRCA1/BRCA1_editing_data.txt"
#for determining whether or not to include a read in output of sam_to_edits
reads_threshold=".000002"
alignment_score_threshold="300"
#Is gene on "+" or "-" strand of genome (e.g. "-" for BRCA1)
orientation="-"
#Consensus coding sequence for gene, should match what's used in ClinVar
CCDS="CCDS11456.2"
#A tab-delimited text file listing the amplicons and their genomic coordinates (must match CADD human genome version, e.g. hg19 or hg38)
amplicon_coords="/camp/home/dacep/home/shared/projects/SGE/BRCA1/amplicon_coords_hg19.txt"
#directory with cadd files for each amplicon, named as {amplicon}.cadd e.g. BRCA1x11.cadd
cadd_dir="/camp/home/dacep/home/shared/projects/SGE/BRCA1/cadd"
#Must download a file from clinvar with all variants selected (e.g. at least 1 star on this date, and then convert it's format to tab-delimited txt file and reference below)
cadd_version="cadd.1.6" #this is the extension used on the cadd files -- currently 'cadd.1.6' e.g. BRCA1x11.cadd.1.6
#clinvar_file generation notes:  Access Clinvar on 7/15/2021. Searched "BRCA1". Filtered on "Single nucleotide" (8894 results). Downloaded as .txt file, imported to excel, deleted columns: 'Protein change', 'Accession', 'dbSNP ID', 'Canonical SPDI'. Created 4 new columns following previous BRCA1 example: cDNAVar, gDNAVar, CADD_Key, and Clinical_Significance. Ensure gDNAVar formula to give complementary base to cDNAVar (because - orientation of gene). Applied formulae to populate these columns to all variants. Saved as an excel file: 
clinvar_file="/camp/home/dacep/home/shared/projects/SGE/BRCA1/ClinVar_BRCA1_210715_Single_nucleotide.txt"
min_indel_freq=".00002"

ml Python/2.7.18-GCCcore-9.3.0

mkdir $out_dir
mkdir $out_dir/fastq
mkdir $out_dir/fastqc
#include collapsed folder to move over collapsed files from asf
cd $out_dir

#move asf sample sheet into folder

python ~/home/users/findlag/bin/run_move_NS_from_asf.py $seq_dir $asf_sample_sheet
sh run_move_from_asf.sh

cd $out_dir/fastq

sh ~/home/users/findlag/bin/merge_lanes.sh -o $out_dir/fastq_merged_lanes
rm $out_dir/fastq/*fastq.gz
mv $out_dir/fastq_merged_lanes/* $out_dir/fastq
rmdir $out_dir/fastq_merged_lanes


zgrep -c @$seq_type *R1_* | tee read_counts.txt
python ~/home/users/findlag/bin/read_count_fractions.py read_counts.txt read_counts2.txt
mv read_counts2.txt read_counts.txt


mkdir Seqprep
mkdir Seqprep/R1
mkdir Seqprep/R2
mkdir Seqprep/merged
python ~/home/users/findlag/bin/run_seqprep_210729.py $adapter1 $adapter2

echo "Running Seqprep on all samples in parallel."
sh run_seqprep.sh
cd Seqprep/merged 

echo "Seqprep done."
zgrep -c @$seq_type *.fastq.gz >> seqprep_read_counts.txt

mkdir no_Ns
python ~/home/users/findlag/bin/run_remove_n_bases.py #operates with getcwd() uses remove_n_bases.py script - this step seems only necessary for some NS runs
sh run_remove_n_bases.sh
echo "remove_n_bases done."
cd no_Ns


cDNA_info_file="/camp/home/dacep/home/shared/projects/SGE/BRCA1/BRCA1_cDNA_data.txt"
python ~/home/users/findlag/bin/run_cDNA_to_gDNA_SGE_pipeline.py $cDNA_info_file
sh run_cDNA_to_gDNA.sh #a script to work on all samples with RNA in the name, to check for perfect ends flanking exons and if they are present, replace with intronic sequence using a hard-coded reference file (cDNA_info_file) specifying each amplicon's 5'/3' cDNA ends, full length, and 5'/3' intronic flanks

mkdir sam
python ~/home/users/findlag/bin/run_needle_to_sam_VHL_pipeline.py $ref_folder

cd $out_dir/fastq/Seqprep/merged/no_Ns
ml purge
module load EMBOSS/6.6.0-foss-2016b

echo "Running needleall."

sh run_needle_to_sam.sh 
wait 
echo "Finished running needleall."

cd sam
mkdir cigar_counts
ml Python/2.7.18-GCCcore-9.3.0

python ~/home/users/findlag/bin/210217_SGE_pipeline_cigar_analyzer.py $cigar_comparisons #not informative for RNA based on full length processing
head -15 cigar_counts/*.txt >> cigar_counts/combined_cigar_counts.txt
mkdir variant_counts
#writes minimally thresholded output (set in hardcode (0.000001 in any) to variant_counts folder


python ~/home/users/findlag/bin/210317_sam_to_edits_pipeline.py $amplicon_list $exp_groupings $ref_folder $editing_info_file $reads_threshold $alignment_score_threshold

cd variant_counts
mkdir final

python ~/home/users/findlag/bin/210413_annotate_variants_pipeline.py $amplicon_list $ref_folder $editing_info_file $orientation $CCDS $amplicon_coords $cadd_dir $cadd_version $clinvar_file 
cd final
head -6 *summary.txt >> combined_editing_data.txt
