#!/usr/bin/bash
#
# VERSION 3
# Usage : bash PIPELINE_v3.sh uparse_config1.txt
#
#########################################################
	##												##
	##												##
	##			 UPARSE pipeline					##
	##												##
	##												##
#########################################################

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Source the Config file >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# absolute path prefered
#source /anas/roshan-current/i2mc/burcelin/frm/scripts/uparse_config1.txt
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
	echo -ne `date | cut -d" " -f1,2,3,4 | sed s'/ /_/g' `" : "
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
#base_fname="ssps"

# creating subdirs
ar=( qualityCheck preprocess clusterOTU chimeraCheck OTUprocess biomFiles TaxAssigned )
#echo ${ar[@]}
for i in ${ar[@]}; do mkdir $dir_out/$i ; done;
#rmdir qualityCheck preprocess clusterOTU chimeraCheck OTUprocess biomFiles

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

# outfiles - TaxAssigned
tf1=`basename  $file_labeledOTUs | sed s'/_otus.fasta/_otus_tax_assignments.txt/'`
file_taxa_assigned=$dir_out/TaxAssigned/$tf1


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Logging >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

# log file
logf=uparse_`date | cut -d" " -f1,2,3,4 | sed s'/ /_/g'`.log
logfile=$dir_out/$logf
#echo $logfile

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Basic Check >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

#Do a version check for the program
echo -ne "Usearch version : \t" &>>$logfile; echo `$u2 --version` | cut -d" " -f2 | cut -d. -f1 | sed s'/v//' &>>$logfile

echo 
echo "Absolute Paths to all the input params : " &>>$logfile
af1=( file_fastq dir_out file_db1 file_97fasta file_97tax py_scripts )
for i in ${af1[@]}; do echo $dir_out/$i &>>$logfile; done;
echo
echo "Absolute Paths to all the dirs and files to be created : " &>>$logfile
af2=( $file_trim250fa $file_trim250fq $file_derep $file_filtered_singles $file_otus $file_chimeras $file_nonchimeras $file_uchimeout $file_labeledOTUs $file_map_uc $file_otu_table $file_otuTableBiom $file_otuTableBiom_wtaxa )
for i in ${af2[@]}; do echo $i &>>$logfile; done;

exit
echo -ne "\n................UPARSE PIPELINE................\n"
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Preprocess >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
echo -ne "\n................Preprocessing................\n"

log_msg "" &>>$logfile
log_msg "UPARSE Pipeline" &>>$logfile
log_msg "" &>>$logfile
# STEP 1:  Check the number of reads in $file_fastq
log_msg "STEP 1 : FASTQ READ COUNT" &>>$logfile
log_msg `fastq_count $file_fastq` &>>$logfile
# log_msg "STEP 1a : QUALITY CHECK FASTQ" &>>$logfile
# log_msg `fastqp $file_fastq `

# STEP 2:  Length Truncation to 250 
log_msg "STEP 2 : TRIM 250"  &>>$logfile
$u2 -fastq_filter $file_fastq  -fastq_trunclen 250 -fastaout $file_trim250fa -fastqout $file_trim250fq  &>>$logfile

# get quality stats
#$u2 -fastq_eestats reads_trimmed.fastq -log reads_trimmed.log

# STEP 3: Depreplication
log_msg "STEP 3 : DEREPLICATION"  &>>$logfile
$u2 -derep_fulllength $file_trim250fa -fastaout $file_derep -sizeout &>>$logfile


log_msg "STEP 3a : FASTA READ COUNT - DEREPLICATED " &>>$logfile
log_msg `fasta_count $file_trim250fa` &>>$logfile

# STEP 4: Discard singletons 
log_msg "STEP 4 : DISCARD SINGLETONS"  &>>$logfile
$u2 -sortbysize $file_derep -fastaout $file_filtered_singles -minsize 2  &>>$logfile

log_msg "STEP 4a : FASTA READ COUNT - SINGLETON REMOVED" &>>$logfile
log_msg `fasta_count $file_filtered_singles` &>>$logfile

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> cluster OTU >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
echo -ne "\n................OTU clustering................\n"

# STEP 5 : Cluster-OTU
# This creates two files - a fasta file and a tab delimited file otu file
log_msg "STEP 5 : CLUSTER OTUs"  &>>$logfile
$u2 -cluster_otus $file_filtered_singles -relabel OTU_ -sizeout -otus $file_otus -uparseout $file_otu_results  &>>$logfile
log_msg "STEP 5a : FASTA READ COUNT - CLUSTERED" &>>$logfile
log_msg `fasta_count $file_otus` &>>$logfile


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Chimera Check >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
echo -ne "\n................Reference based Chimera Check................\n"
# STEP 6 :   Chimera filtering using reference database
log_msg "STEP 6 : CHIMERA CHECKING USING REFERENCE DATABASE"  &>>$logfile
$u2 -uchime_ref $file_otus -db $file_db1 -strand plus -threads 30 -uchimeout $file_uchimeout  -chimeras $file_chimeras -nonchimeras $file_nonchimeras  &>>$logfile
log_msg "STEP 6a : FASTA READ COUNT of NonCHIMERIC SEQS" &>>$logfile
log_msg `fasta_count $file_nonchimeras` &>>$logfile
log_msg "STEP 6b : FASTA READ COUNT of CHIMERIC SEQS" &>>$logfile
log_msg `fasta_count $file_chimeras` &>>$logfile

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> OTU process >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
echo -ne "\n................Processing the OTUs................\n"

# STEP 7 : label OTUs - using python script from uparese pipeline
# OTU_ can be repaced by any suitable name
echo -ne "\nOTU Labeling of nonchimeric fasta file ................\n"
log_msg "STEP 7 : LABELING THE NonCHIMERIC FASTA WITH OTU_ "  &>>$logfile
python $py_scripts/fasta_number.py $file_nonchimeras OTU_ > $file_labeledOTUs  

# STEP 8 : Map reads (including singletons) back to OTUs
# Because the names has been stripped so putting it back
echo -ne "\nMapping the OTUS back to the original trimmed fasta..................\n"
log_msg "STEP 8 : MAP THE TRIME250_READS TO LABELLED OTU FASTA"  &>>$logfile
$u2 -usearch_global $file_trim250fa -db $file_labeledOTUs -strand plus -id 0.97 -uc $file_map_uc  

# STEP 9 : Create OTU table
echo -ne "\nCreate OTU Table..This..May..Take..Some..Time...............\n"
log_msg "STEP 9 : Create OTU Table"  &>>$logfile
python $py_scripts/uc2otutab_modified.py $file_map_uc > $file_otu_table   


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> BiomFiles >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
echo -ne "\n................Create Biom File ...............\n"


# STEP 10 : Create a  biom file from file_otu_table
# converting into json format
log_msg "STEP 10 : Create Biom file"  &>>$logfile
biom convert --table-type="OTU table" -i $file_otu_table -o $file_otuTableBiom --to-json  &>>$logfile

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Taxonomy Assignment >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
echo -ne "\n................Assigning Taxonomy to OTU table ...............\n"

# STEP 11 : Assign taxonomy to the OTUs
# require gg_97_tax and gg_97_fasta 
# not using usearch -utax option as no taxconfs and taxtree is there

log_msg "STEP 11 : Assigning Taxonomy to OTU table"  &>>$logfile
assign_taxonomy.py -t $file_97tax -r $file_97fasta -i $file_labeledOTUs -o $dir_out/TaxAssigned   &>>$logfile


#STEP 12 : Adding the taxonomy to BIOM table
log_msg "STEP 12 : Adding Taxonomy to BIOM file"  &>>$logfile
biom add-metadata --sc-separated taxonomy --observation-header OTUID,taxonomy --observation-metadata-fp $file_taxa_assigned -i $file_otuTableBiom -o $file_otuTableBiom_wtaxa

# Convert the biom-json-with-tax file to hdf5 
# After this add-metadata from mapping file to the biom-wtax-hdf5 using --sample-metadata mappingfile.txt
# Later convert it back to tsv file before putting it into phyloseq 
# or try with the hdf5 file itself

# biom convert -i testV_otu_table_withTaxa.biom --table-type="OTU table" --header-key taxonomy --to-tsv -o T1.tsv
