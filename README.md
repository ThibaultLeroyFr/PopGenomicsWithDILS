# PopGenomicsWithDILS

This repository regroups all scripts used to generate the inputfiles in order to perform demographic analyses using DILS (https://github.com/popgenomics/DILS_web) from raw sequencing data to post-analysis R plots, including the mapping, variant calling, fasta sequence reconstruction and demographic inference steps.

### 1/ From raw reads to read mapping (./PipelineMappingCalling)


### 2/ Variant calling (./PipelineMappingCalling)

<code>
<p> listacc=$(echo "pdav73") <p>
refile=$(echo "/sandbox/users/tleroy/Francesca/Potra_genome2.2/Potra02_genome_softmasked.fasta" )  #reference file (need to be indexed => script_index.sh) ! <br/>

pathtodata=$(echo "/sandbox/users/tleroy/Francesca/mapping") # the repertory containing all individus <br/>
pathtoscripts=$(echo "/sandbox/users/tleroy/AfricanRice/scripts/PipelineMappingCalling/") # set path for 1_mapping.sh and 2_snpindel_callingGVCF.sh <br><br>
\# Please change file path in 1_mapping.sh and in 2_snpindel_callingGVCF.sh ! <br/>

module load java # load java if needed for your cluster # GATK requires java8 ! <br/>

cd /sandbox/users/tleroy/Francesca/gvcf/ <br/>
\# CMD bash 2_snpindel_callingGVCF.sh [IDname] [Reference_Genome] [output_directory] [Number_of_CPU_to_use] <br/>
bash $pathtoscripts/2_snpindel_callingGVCF.sh $listacc $refile $pathtodata/ 4 <br/>
</code>

### 3/ Joint Genotyping (./PipelineMappingCalling)



Authors:<br>
Francesca Beclin (Master's project): a01346615_at_unet.univie.ac.at<br>
Thibault Leroy: thibault.leroy_at_univie.ac.at<br>
