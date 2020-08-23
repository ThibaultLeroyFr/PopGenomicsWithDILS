# PopGenomicsWithDILS

This repository regroups all scripts used to generate the inputfiles in order to perform demographic analyses using DILS (https://github.com/popgenomics/DILS_web) from raw sequencing data to post-analysis R plots, including the mapping, variant calling, fasta sequence reconstruction and demographic inference steps.

Authors:<br>
Francesca Beclin (Master's project): a01346615_at_unet.univie.ac.at<br>
Thibault Leroy: thibault.leroy_at_univie.ac.at<br>

### 1/ From raw reads to read mapping (./PipelineMappingCalling)


### 2/ Variant calling (./PipelineMappingCalling)

<code>
 <p> listacc=$(echo "pdav73")  <em> #SampleID </em></p>
<p> refile=$(echo "/sandbox/users/tleroy/Francesca/Potra_genome2.2/Potra02_genome_softmasked.fasta" ) <em> #reference file (need to be indexed => script_index.sh) ! </em></p>

<p> pathtodata=$(echo "/sandbox/users/tleroy/Francesca/mapping") <em> # the repertory containing all individus  </em></p>
<p> pathtoscripts=$(echo "/sandbox/users/tleroy/AfricanRice/scripts/PipelineMappingCalling/")  </p>
<p><em># Please change file path in 1_mapping.sh and in 2_snpindel_callingGVCF.sh ! </em></p>

<p>module load java <em> # load java if needed for your computing cluster (# GATK requires java8) !</em></p>

<p>cd /sandbox/users/tleroy/Francesca/gvcf/ </p>
<p><em># CMD: bash 2_snpindel_callingGVCF.sh [SampleIDn] [Reference_Genome] [output_directory] [Number_of_CPU_to_use] </em></p>
<p> bash $pathtoscripts/2_snpindel_callingGVCF.sh $listacc $refile $pathtodata/ 4 </p>
</code>

### 3/ Joint Genotyping (./PipelineMappingCalling)



