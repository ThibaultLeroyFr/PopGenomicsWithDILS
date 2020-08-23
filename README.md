# PopGenomicsWithDILS

This repository regroups all scripts used to generate the inputfiles in order to perform demographic analyses using DILS (https://github.com/popgenomics/DILS_web) from raw sequencing data to post-analysis R plots, including the mapping, variant calling, fasta sequence reconstruction and demographic inference steps.

Authors:<br>
Francesca Beclin (Master's project): a01346615_at_unet.univie.ac.at<br>
Thibault Leroy: thibault.leroy_at_univie.ac.at<br>

### 1/ From raw reads to read mapping (./PipelineMappingCalling)

<ins>1.1 Trimming</ins><br>
java -jar /home/fs71105/francescab/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 1 -phred33 /home/fs71105/francescab/pdav/pdav73R1.fastq.gz /home/fs71105/francescab/pdav/pdav73R2.fastq.gz /home/fs71105/francescab/trimming/pdav73_1_cleaned.fastq.gz /home/fs71105/francescab/trimming/pdav73_1_cleaned_unpaired.fastq.gz /home/fs71105/francescab/trimming/pdav73_2_cleaned.fastq.gz /home/fs71105/francescab/trimming/pdav73_2_cleaned_unpaired.fastq.gz  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50

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


### 2/ Variant calling (./PipelineMappingCalling)

<p> listacc=$(echo "pdav73")  <em> #SampleID </em></br>
 refile=$(echo "/sandbox/users/tleroy/Francesca/Potra_genome2.2/Potra02_genome_softmasked.fasta" ) <em> #reference file (need to be indexed => script_index.sh) ! </em><br>
pathtodata=$(echo "/sandbox/users/tleroy/Francesca/mapping") <em> # the repertory containing all individus  </em><br>
 pathtoscripts=$(echo "/sandbox/users/tleroy/AfricanRice/scripts/PipelineMappingCalling/") <br>
<em># Please change file path in 1_mapping.sh and in 2_snpindel_callingGVCF.sh ! </em><br>

module load java <em> # load java if needed for your computing cluster (# GATK requires java8) !</em><br>

cd /sandbox/users/tleroy/Francesca/gvcf/ <br>
<em># CMD: bash 2_snpindel_callingGVCF.sh [SampleIDn] [Reference_Genome] [output_directory] [Number_of_CPU_to_use] </em><br>
bash $pathtoscripts/2_snpindel_callingGVCF.sh $listacc $refile $pathtodata/ 4 <br><p>


### 3/ Joint Genotyping & SNP filtering (./PipelineMappingCalling)




### 4/ Reconstructing fasta sequences & extracting genomic blocks (or CDS) (./)
Note: Depending of the number of individuals and the length of the longest scaffold (or chromosome), this step can be memory-intensive. 



### 5/ Computing summary statistics (./)


### 6/ Generating input files for DILS (./)
