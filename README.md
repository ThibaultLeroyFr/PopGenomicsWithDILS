# PopGenomicsWithDILS

This repository regroups all scripts used to generate the inputfiles in order to perform demographic analyses using DILS (https://github.com/popgenomics/DILS_web) from raw sequencing data to post-analysis R plots, including the mapping, variant calling, fasta sequence reconstruction and demographic inference steps.

Authors:<br>
Francesca Beclin (Master's project): a01346615_at_unet.univie.ac.at<br>
Thibault Leroy: thibault.leroy_at_univie.ac.at<br>

### 1/ From raw reads to read mapping (./PipelineMappingCalling)

<ins>1.1 Trimming</ins><br>
<em>java -jar /home/fs71105/francescab/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads [#CPUs] -phred33 [inputfile_rawreads_1] [inputfile_rawreads_2] [output_trimmed_paired_1] [output_trimmed_unpaired_1] [output_trimmed_paired_2] [output_trimmed_unpaired_2] [PARAMETERS_TRIMMING] </em><br>
e.g.: java -jar /home/fs71105/francescab/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 1 -phred33 /home/fs71105/francescab/pdav/pdav73R1.fastq.gz /home/fs71105/francescab/pdav/pdav73R2.fastq.gz /home/fs71105/francescab/trimming/pdav73_1_cleaned.fastq.gz /home/fs71105/francescab/trimming/pdav73_1_cleaned_unpaired.fastq.gz /home/fs71105/francescab/trimming/pdav73_2_cleaned.fastq.gz /home/fs71105/francescab/trimming/pdav73_2_cleaned_unpaired.fastq.gz  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50

<ins>1.2 Indexing references files</ins><br>
<p><em># bwa</em><br> 
/home/fs71105/francescab/software/bwa-0.7.17/bwa index /home/fs71105/francescab/Potra_genome2.2/Potra02_genome_softmasked.fasta <br>
<em># samtools</em><br> 
/home/fs71105/francescab/software/samtools/samtools faidx /home/fs71105/francescab/Potra_genome2.2/Potra02_genome_softmasked.fasta <br> 
<em># picard</em><br> 
java -jar /home/fs71105/francescab/software/picard/picard.jar CreateSequenceDictionary <br> 
      R=/home/fs71105/francescab/Potra_genome2.2/Potra02_genome_softmasked.fasta <br> 
      O=/home/fs71105/francescab/Potra_genome2.2/Potra02_genome_softmasked.dict <br> </p>

<ins>1.3 Mapping</ins><br>


### 2/ Generating gvcf files (./PipelineMappingCalling)

<p> listacc=$(echo "pdav73")  <em> #SampleID </em></br>
 refile=$(echo "/sandbox/users/tleroy/Francesca/Potra_genome2.2/Potra02_genome_softmasked.fasta" ) <em> #reference file (need to be indexed => script_index.sh) ! </em><br>
pathtodata=$(echo "/sandbox/users/tleroy/Francesca/mapping") <em> # the repertory containing all individus  </em><br>
 pathtoscripts=$(echo "/sandbox/users/tleroy/AfricanRice/scripts/PipelineMappingCalling/") <br>
<em># Please read the file 2_snpindel_callingGVCF.sh (& change file & programs paths)! </em><br>

module load java <em> # load java if needed for your computing cluster (# GATK requires java8) !</em><br>

cd /sandbox/users/tleroy/Francesca/gvcf/ <br>
<em># CMD: bash 2_snpindel_callingGVCF.sh [SampleIDn] [Reference_Genome] [output_directory] [Number_of_CPU_to_use] </em><br>
bash $pathtoscripts/2_snpindel_callingGVCF.sh $listacc $refile $pathtodata/ 4 <br><p>


### 3/ Joint Genotyping & SNP filtering (./PipelineMappingCalling)
<ins>3.1 Joint Genotyping (using intervals to perform this step on several CPUs)</ins><br>
<em>bash 3_intervals_jointgenotyping.sh [refernece_genome] [#CPUs]</em><br>
<em>e.g. bash 3_intervals_jointgenotyping.sh /sandbox/users/tleroy/Francesca/Potra_genome2.2/Potra02_genome_softmasked.fasta 10</em><br>
(see "3_intervals_jointgenotyping.sh" for details) <br>

<ins>3.2 Merging the outputs</ins><br>
intervals_used=$(echo "9")<br>
tmp_dir=$(echo "/sandbox/users/tleroy/Francesca/gvcf/tmp_vcf_Potra02_genome_softmasked.fasta")<br>
cd $tmp_dir<br>
for i in $(seq 1 $intervals_used); do file=$(echo "scatter""$i"".intervals.joint_bwa_mem_mdup_raw.vcf"); if [ $i == 1 ]; then cp $file merged_joint_bwa_mem_mdup_raw.vcf; else grep -v "#" $file >> merged_joint_bwa_mem_mdup_raw.vcf; fi; done<br>

<ins>3.3 Variant filtering </ins><br>
VariantFiltrationGVCF.py -q 2.0 -s 60.0 -m 40.0 -n -2.0 -r -2.0 -w 45000 -f [infile] > [outfile] <br>
e.g. /sandbox/users/tleroy/Francesca/scripts/VariantFiltrationGVCF.py -q 2.0 -s 60.0 -m 40.0 -n -2.0 -r -2.0 -w 45000 -f merged_joint_bwa_mem_mdup_raw.vcf > merged_joint_bwa_mem_mdup_raw.filtered.vcf <br>

### 4/ Reconstructing fasta sequences & extracting genomic blocks (or CDS) (./Generate_Sequence_Blocks)
Note: Depending of the number of individuals and the length of the longest scaffold (or chromosome), this step can be memory-intensive. 

<ins>4.1 Reconstruct fasta sequences</ins>

<p<em> File & directory names (full path)</em></br>
vcffile=$(echo "/bigvol/benoit/Analyses/Temp_Tibo/Francesca/joint_pdav/merged_joint_bwa_mem_mdup_raw.filtered.vcf")</br>
gfffile=$(echo "/bigvol/benoit/Analyses/Temp_Tibo/Francesca/Potra_genome2.2/Potra02_genes.gff.clean")</br>
outputdirscaffolds=$(echo "/bigvol/benoit/Analyses/Temp_Tibo/Francesca/joint_pdav/Pdavidiana_fasta_files_withoutquantiles_scaffold")</br>
cutoffqualitybases=$(echo "20") <em> # minimum illumina quality at the base position, here >= 20 </em></br>
cutoffcovmin=$(echo "3") <em> # minimum coverage per individual (cov < minimum => position will be hard masked = "N") </em></br>
cutoffcovmax=$(echo "50") <em> # maximum coverage per individual (cov > maximum => hard masked) </em></br>
outprefix=$(echo "Pdavidiana_withoutquantiles")</em></br>

<em> Main script (see script_VCF2Fasta_withcovqual.sh for all details)</em></br>
if [ -d "$outputdirscaffolds" ]; then</br>
___rm $outputdirscaffolds/*.fst</br>
else</br>
___mkdir "$outputdirscaffolds"</br>
fi</br>
bash /home/thibault/scripts/script_VCF2Fasta_withcovqual.sh $vcffile $outputdirscaffolds $cutoffqualitybases $cutoffcovmin $cutoffcovmax</br>

<ins>4.2 Extract genomic windows </ins>
<em>bash generate_slidwin.sh [assembly_scafflength] [input_dir_containing_scaffold.fst_files] [windows_size] [prefix] </em></br>
e.g bash generate_slidwin.sh /bigvol/benoit/Analyses/Temp_Tibo/Francesca/Potra_genome2.2/Potra02_genome_softmasked.fasta.scafflength Pdavidiana_fasta_files_withoutquantiles_scaffold 10000 Pdavidiana_pi <br>

The "scafflength" file is just a file containing the names of the scaffold (1st column) and the length of the corresponding scaffold (2nd column) and is expected to be generated by the python script during the step 4.1 (if not, execute python ./PipelineMappingCalling/script_scaff_length.py [input_fasta] > [output.scafflenth]). <br>
    
<ins>4.3 Extract some specific features, such as CDS or genes </ins>

<em>e.g. GET FASTA ON CDS</em><br>
gfffile=$(echo "/bigvol/benoit/Analyses/Temp_Tibo/Francesca/Potra_genome2.2/Potra02_genes.gff.clean")<br>
inputdirscaffolds=$(echo "/bigvol/benoit/Analyses/Temp_Tibo/Francesca/joint_pdav/Pdavidiana_fasta_files_withoutquantiles_scaffold")<br>
outprefix=$(echo "Pdavidiana_withoutquantiles")<br>

awk '$3 == "CDS" {print $0}' $gfffile | awk '{print $1}' | sort | uniq > $gfffile.scaffIDwithCDS<br>
outputdirCDS=$(echo "$inputdirscaffolds" | sed 's/scaffold/CDS/g' | sed 's/chromosome/CDS/g')<br>
if [ -d "$outputdirCDS" ]; then<br>
____rm $outputdirCDS/*.fst<br>
else<br>
____mkdir "$outputdirCDS"<br>
fi<br>
cd $outputdirCDS<br>
while read line; do python /media/bigvol/benoit/Scripts/cutSeqGff.py $outputdirscaffolds/$line.fst $gfffile $line CDS; done < $gfffile.scaffIDwithCDS<br>
cd ..<br>

<em>The gff files need to be as follows (the tag Name=XXX is particularly important; + same names for all CDS exons of the same gene)  </em></br>
<em>chr1    maker   gene    8865    11259   .       -       .       Name=Potra2n1c1;  </em></br>
<em>chr1    maker   mRNA    8865    10802   .       -       .       Name=Potra2n1c1.3; </em></br>
<em>chr1    maker   CDS     8865    9054    .       -       1       Name=Potra2n1c1.3; </em></br>
<em>chr1    maker   CDS     9487    9559    .       -       2       Name=Potra2n1c1.3; </em></br>
<em>chr1    maker   CDS     9669    9753    .       -       0       Name=Potra2n1c1.3; </em></br>

One file is expected to be generated per geneID (geneID.fst). Importantly, check the output of the different fasta files for CDS, most genes are expected to start by an "ATG" and to end by a stop codon. If it is not the case, try to use the "cutSeqGff_dec1bp.py" script rather than "cutSeqGff.py" to see if the issue can be fixed (e.g. this is typically needed for the Populus tremula genome v.2.2). <br>

### 5/ Computing summary statistics (./Compute_SumStats)

<ins>5.1 piN, piS and piN/piS ratiosn on CDS</ins>
<em>FILTER ALIGNMENTS & COMPUTE STATS</em><br>
gfffile=$(echo "/bigvol/benoit/Analyses/Temp_Tibo/Francesca/Potra_genome2.2/Potra02_genes.gff.clean")<br>
inputdirscaffolds=$(echo "/bigvol/benoit/Analyses/Temp_Tibo/Francesca/joint_pdav/Pdavidiana_fasta_files_withoutquantiles_scaffold")<br>
outprefix=$(echo "Pdavidiana_withoutquantiles")<br>

outputdirCDS=$(echo "$inputdirscaffolds" | sed 's/scaffold/CDS/g' | sed 's/chromosome/CDS/g')<br>
<em>\# list of all fasta CDS</em><br>
ls $outputdirCDS/ | grep ".fst" > $outprefix.list_CDS.txt<br>
cd $outputdirCDS<br>
<em>\# remove last codon</em><br>
/home/thibault/scripts/CDS_alignments_piNpiS/removeLastStopCodon -seq ../$outprefix.list_CDS.txt -f fasta -code univ<br>
cd ..<br>
<em>\# generate a new list with processed alignments</em><br>
ls $outputdirCDS/ | grep ".fst.clean.fst" > $outprefix.list_CDS.txt<br>
<em>\# clean alignments</em><br>
cd $outputdirCDS<br>
/home/thibault/Clean_Alignment -seq ../$outprefix.list_CDS.txt -f fasta -n 4<br>
cd ..<br>
<em>\# generate a new list with processed alignments</em><br>
ls $outputdirCDS/ | grep ".fst.clean.fst.clean.fst" > $outprefix.list_CDS.txt<br>
<em>\# compute summary statistics (-tstv = transition transervision ratio here fixed to 2 but can be set to another value)</em><br>
cd $outputdirCDS<br>
/home/thibault/scripts/CDS_alignments_piNpiS/seq_stat_coding -seq ../$outprefix.list_CDS.txt -f fasta -tstv 2 -code univ -o ../$outprefix.CDS.sumstats > ../$outprefix.CDS.sumstats.info<br>
cd ..<br>
<em>\# keep info of genes without premature stop codons</em><br>
grep "stop" $outprefix.CDS.sumstats.info | awk '{print $1}' > $outprefix.CDS.withprematurestopcodons<br>
<em>\# Compute piS, piN and piN/piS over the whole set of genes </em><br>
bash script_generate_piNpiS.sh $outprefix.CDS.sumstats<br>

### 6/ Generating input files for DILS

<em>The input files need to be formated as follows:</em></br>
<em>>GeneID|SpeciesName|sampleID|All1</em></br>
<em>ATGCGC...</em></br>
<em>>GeneID|SpeciesName|sampleID|All2</em></br>
<em>ATGCAC...</em></br></p>

<p><em>For example, to generate 10 sets of 1000 CDS randomly selected</em></br>
<em>Randomly sample 1000 genes among all the gene sequences (geneID.fst files). </em></br>
<em>List all sampled files in a "CDSsampling1.txt", do this step 10 times.</em></br></br>
for i in {1..10}; do</br>
____rm Sampling_1000CDSrandom_$i.fas</br>
____echo "generating file Sampling_1000CDSrandom_$i.fas"</br>
____while read line; do</br>
________locus=$(basename "$line" ".fst")</br>
________echo -n "." </br>
________less ./CDS_sampling$i/$line | sed "s/pdav/$locus|Populus_davidiana|pdav/g"  | sed 's/\.1$/|All1/g' | sed 's/\.2$/|All2/g' >> Sampling_1000CDSrandom_$i.fas</br>
____done < CDSsampling$i.txt</br>
done</br></p>
