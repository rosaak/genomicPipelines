#!/usr/bin/bash
#
# VERSION 5
# Usage : bash uparse_pipeline.sh uparse_config.txt
# roshan.padmanabhan@inserm.fr
# 
#########################################################
##				
##					
##  UPARSE pipeline
##
##
#########################################################
#
# Requirements : usearch, vsearch, qiime scripts, fastqp
#
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Source the Config file >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# absolute path prefered
if [[ $1 == '' ]]; then
	echo "Please provide a config file"
	exit
else 
	source $1
fi
# checking if the dir_out is 
if [[ ! -d "$dir_out"  ]]; then
       #echo "creating an results out dir"
        mkdir $dir_out
fi

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Fucntions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

function log_msg () {
	echo -ne "[" `date | cut -d" " -f1,2,3,4 | sed s'/ /_/g' `" ]  "
	echo $@ 
}

function fasta_count (){
	grep "^>" -c $1
}

function fastq_count (){
	wc -l $1 | awk '!/total/{print $2"\t"$1/4}'
}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Dris >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# Baseflename 
base_fname=`basename $file_fastq  | cut -d. -f1`
#base_fname="pome1"

# creating subdirs
ar=( qualityCheck preprocess clusterOTU chimeraCheck OTUprocess biomFiles TaxAssigned )
#echo ${ar[@]}
for i in ${ar[@]}; do mkdir $dir_out/$i ; done;
#rmdir qualityCheck preprocess clusterOTU chimeraCheck OTUprocess biomFiles

# outfiles - qualityCheck
file_qc=$dir_out/qualityCheck/$base_fname

# outfiles - preprocess
file_trim250fa=$dir_out/preprocess/$base_fname"_trim250.fasta"
file_trim250fq=$dir_out/preprocess/$base_fname"_trim250.fastq"
file_derep=$dir_out/preprocess/$base_fname"_derep.fna"
file_filtered_singles=$dir_out/preprocess/$base_fname"_singletonRemoved.fna"

# outfiles - clusterOTU
file_otus=$dir_out/clusterOTU/$base_fname"_otu.fna"
file_otu_results=$dir_out/clusterOTU/$base_fname"_otu.up"

# outfiles - chimeraCheck
file_chimeras=$dir_out/chimeraCheck/$base_fname"_chimeric.fna"
file_nonchimeras=$dir_out/chimeraCheck/$base_fname"_chimeraRemoved.fna"
file_uchimeout=$dir_out/chimeraCheck/$base_fname"_results.uchime"

# outfiles - OTUprocess
file_labeledOTUs=$dir_out/OTUprocess/$base_fname"_labedledNonChimeric_otus.fasta"
file_map_uc=$dir_out/OTUprocess/$base_fname"_otuMap.uc"
file_otu_table=$dir_out/OTUprocess/$base_fname"_otuTable.txt"

# outfiles - biomFiles
file_otuTableBiom=$dir_out/biomFiles/$base_fname"_otu_table.biom"
file_otuTableBiom_wtaxa=$dir_out/biomFiles/$base_fname"_otu_table_withTaxa.biom"
file_otuTable_with_Tax_and_Metadata=$dir_out/biomFiles/$base_fname"_otu_table_withTaxa_and_MetaData.biom"
file_otuTable_with_Tax_and_Metadata_tsv=$dir_out/biomFiles/$base_fname"_otu_table_withTaxa_and_MetaData.tsv"

# outfiles - TaxAssigned
tf1=`basename  $file_labeledOTUs | sed s'/_otus.fasta/_otus_tax_assignments.txt/'`

# TOGGLE 

# by gg135
file_taxa_assigned=$dir_out/TaxAssigned_gg135/$tf1
file_utaxa_assigned=$dir_out/TaxAssigned_gg135/$base_fname"_utax_gg135_assigned.txt"

# by gg138
#file_taxa_assigned=$dir_out/TaxAssigned_gg138/$tf1
#file_utaxa_assigned=$dir_out/TaxAssigned_gg138/$base_fname"_utax_gg138_assigned.txt"

# by s119
#file_taxa_assigned=$dir_out/TaxAssigned_s119/$tf1
#file_utaxa_assigned=$dir_out/TaxAssigned_s119/$base_fname"_utax_s119_assigned.txt"


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Logging >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

# log file
logf=uparse_`date | cut -d" " -f1,2,3,4 | sed s'/ /_/g'`.log
logfile=$dir_out/$logf
#echo $logfile

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Basic Check >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
echo -ne "\n................UPARSE PIPELINE................\n" &>>$logfile
echo -ne "\n................Basic Check................\n" &>>$logfile
#Do a version check for the program
echo -ne "Usearch version : \t" &>>$logfile; echo `$u --version` | cut -d" " -f2 | cut -d. -f1 | sed s'/v//' &>>$logfile
echo -ne "Vsearch version : \t" &>>$logfile; echo `$v --version` | cut -d"," -f1,2 | sed s'/v//' &>>$logfile
echo "Absolute Paths to all the input params : " &>>$logfile
af1=( file_fastq dir_out file_db1 file_gg135_97fasta file_gg135_97tax py_scripts file_meta_data )
for i in ${af1[@]}; do echo $dir_out/$i &>>$logfile; done;
echo "Absolute Paths to all the dirs and files to be created : " &>>$logfile
af2=( $file_trim250fa $file_trim250fq $file_derep $file_filtered_singles $file_otus $file_chimeras $file_nonchimeras $file_uchimeout $file_labeledOTUs $file_map_uc $file_otu_table $file_otuTableBiom $file_otuTableBiom_wtaxa $file_otuTable_with_Tax_and_Metadata $file_otuTable_with_Tax_and_Metadata_tsv )
for i in ${af2[@]}; do echo $i &>>$logfile; done;

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Preprocess >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
echo -ne "\n................Preprocessing................\n" &>>$logfile
log_msg " " &>>$logfile
log_msg "STEP 1 : FASTA READ COUNT" &>>$logfile
log_msg `fasta_count $file_fasta` &>>$logfile
log_msg " " &>>$logfile
# If fastq is provided then I can enable this step to check the quality of the fastq file
#log_msg "STEP 1a : QUALITY CHECK FASTQ" &>>$logfile
#log_msg `fastqp -o $file_qc $file_fastq`  &>>$logfile

log_msg "STEP 2 : DEREPLICATION"  &>>$logfile
# uparse derep can't handle bcoz of 32 bit so this step has been replaced by vsearch dereplication command 
#$u -derep_fulllength $file_fasta -fastaout $file_derep -sizeout &>>$logfile
$v --derep_fulllength $file_fasta --output   $file_derep -sizeout &>>$logfile
log_msg " " &>>$logfile
log_msg "STEP 2a : FASTA READ COUNT - DEREPLICATED " &>>$logfile
log_msg `fasta_count $file_derep` &>>$logfile
log_msg " " &>>$logfile
log_msg "STEP 3 : DISCARD SINGLETONS"  &>>$logfile
$u -sortbysize $file_derep -fastaout $file_filtered_singles -minsize 2  &>>$logfile
log_msg " " &>>$logfile
log_msg "STEP 3a : FASTA READ COUNT - SINGLETON REMOVED" &>>$logfile
log_msg `fasta_count $file_filtered_singles` &>>$logfile
log_msg " " &>>$logfile

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> cluster OTU >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
echo -ne "\n................OTU clustering................\n" &>>$logfile
log_msg "STEP 4 : CLUSTER OTUs"  &>>$logfile
$u -
cluster_otus $file_filtered_singles -relabel OTU_ -sizeout -otus $file_otus -uparseout $file_otu_results  &>>$logfile
log_msg " " &>>$logfile
log_msg "STEP 4a : FASTA READ COUNT - CLUSTERED" &>>$logfile
log_msg `fasta_count $file_otus` &>>$logfile
log_msg " " &>>$logfile

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Chimera Check >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
echo -ne "\n................Reference based Chimera Check................\n" &>>$logfile
log_msg "STEP 5 : CHIMERA CHECKING USING REFERENCE DATABASE"  &>>$logfile
$u -uchime_ref $file_otus -db $file_db1 -strand plus -threads 30 -uchimeout $file_uchimeout  -chimeras $file_chimeras -nonchimeras $file_nonchimeras  &>>$logfile
log_msg " " &>>$logfile
log_msg "STEP 5a : FASTA READ COUNT of NonCHIMERIC SEQS" &>>$logfile
log_msg `fasta_count $file_nonchimeras` &>>$logfile
log_msg " " &>>$logfile
log_msg "STEP 5b : FASTA READ COUNT of CHIMERIC SEQS" &>>$logfile
log_msg `fasta_count $file_chimeras` &>>$logfile
log_msg " " &>>$logfile

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> OTU process >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
echo -ne "\n................Processing the OTUs................\n" &>>$logfile
echo -ne "\nOTU Labeling of nonchimeric fasta file\n" &>>$logfile
log_msg " " &>>$logfile
log_msg "STEP 6 : LABELING THE NonCHIMERIC FASTA WITH OTU_ "  &>>$logfile
python $py_scripts/fasta_number.py $file_nonchimeras OTU_ > $file_labeledOTUs  
log_msg " " &>>$logfile
echo -ne "\nMapping the OTUS back to the original trimmed fasta\n" &>>$logfile
log_msg " " &>>$logfile
log_msg "STEP 7 : MAPPING THE OTUS BACK TO THE ORIGINAL FASTA"  &>>$logfile
$u -usearch_global $file_fasta -db $file_labeledOTUs -strand plus -id 0.97 -uc $file_map_uc  &>>$logfile
log_msg " " &>>$logfile
echo -ne "\nCreate OTU Table\n" &>>$logfile
log_msg " " &>>$logfile
log_msg "STEP 8 : Create OTU Table"  &>>$logfile
python $py_scripts/uc2otutab_modified.py $file_map_uc > $file_otu_table   
log_msg " " &>>$logfile

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BiomFiles >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
echo -ne "\n................Create Biom File ...............\n" &>>$logfile
log_msg "STEP 9 : Create Biom file"  &>>$logfile
biom convert --table-type="OTU table" -i $file_otu_table -o $file_otuTableBiom --to-json  &>>$logfile
log_msg " " &>>$logfile

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Taxonomy Assignment >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
echo -ne "\n................Assigning Taxonomy to OTU table and Preparing Biom files...............\n" &>>$logfile

#
# require gg_97_tax and gg_97_fasta 
# not using usearch -utax option

# TOGGLE 

log_msg "STEP 10 : Assigning Taxonomy to OTU table with gg135"  &>>$logfile
$assigntaxonomy -t $file_gg135_97tax -r $file_gg135_97fasta -i $file_labeledOTUs -o $dir_out/TaxAssigned_gg135   &>>$logfile

#log_msg "STEP 10 : Assigning Taxonomy to OTU table with gg138"  &>>$logfile
#$assigntaxonomy -t $file_gg138_97tax -r $file_gg138_97fasta -i $file_labeledOTUs -o $dir_out/TaxAssigned_gg138   &>>$logfile

#log_msg "STEP 10 : Assigning Taxonomy to OTU table with s119"  &>>$logfile
#$assigntaxonomy -t $file_s119_97tax -r $file_s119_97fasta -i $file_labeledOTUs -o $dir_out/TaxAssigned_s119   &>>$logfile

# Alternate route 
#usearch -utax reads.fq -db 16s.udb -taxconfs 16s.tc -tt 16s.tt -utaxout results.txt
#$u -utax $file_labeledOTUs -db $u_db16s -tt $u_db16s_tt -utax_rawscore -utaxout $file_utaxa_assigned &>>$logfile
#$p_assigntaxonomy_rdp <-- check whether rdp classifier is installed 
#

log_msg " " &>>$logfile


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Biom Manipulation >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

#STEP 11-14 : Preparing the BIOM table with taxa and metadata information 

log_msg "STEP 11 : Adding Taxonomy to BIOM file"  &>>$logfile
biom add-metadata --sc-separated taxonomy --observation-header OTUID,taxonomy --observation-metadata-fp $file_taxa_assigned -i $file_otuTableBiom -o $file_otuTableBiom_wtaxa  &>>$logfile
log_msg " " &>>$logfile
log_msg "STEP 12:  Adding Meta Data to the BIOM file" &>>$logfile
biom add-metadata -i $file_otuTableBiom_wtaxa -m $file_meta_data -o $file_otuTable_with_Tax_and_Metadata
log_msg " " &>>$logfile
log_msg "STEP 13: Biom conversion into tsv file and changing taxonomy to Consensus Lineage" &>>$logfile
biom convert -i $file_otuTable_with_Tax_and_Metadata -o $file_otuTable_with_Tax_and_Metadata_tsv --to-tsv --header-key taxonomy &>>$logfile
sed -i s'/taxonomy/Consensus Lineage/' $file_otuTable_with_Tax_and_Metadata_tsv
log_msg " " &>>$logfile
log_msg "STEP 14: Summarize table of the final Biom file" &>>$logfile
biom summarize-table -i $file_otuTable_with_Tax_and_Metadata_tsv &>>$logfile
log_msg " " &>>$logfile
log_msg "Done"




