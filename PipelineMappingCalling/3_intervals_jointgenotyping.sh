assembly=$1  # genome assembly
intervals=$2   # expected nb intervals (number of intervals created can be slightly lower or higher)

# tmp dir
assemblyname=$(basename $assembly)
mkdir /sandbox/users/tleroy/Francesca/gvcf/tmp_vcf_$assemblyname

# create intervals
python /sandbox/users/tleroy/Francesca/scripts/script_scaff_length.py $assembly > /sandbox/users/tleroy/Francesca/gvcf/tmp_vcf_$assemblyname/$assemblyname.scaffsize
python /sandbox/users/tleroy/Francesca/scripts/createintervalsfromscaffsize.py /sandbox/users/tleroy/Francesca/gvcf/tmp_vcf_$assemblyname/$assemblyname.scaffsize $intervals
mv scatter*.intervals /sandbox/users/tleroy/Francesca/gvcf/tmp_vcf_$assemblyname/
cd /sandbox/users/tleroy/Francesca/gvcf/tmp_vcf_$assemblyname

# loop over intervals
for i in scatter*.intervals; do
	#echo "time java -jar /media/bigvol/benoit/software/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -T GenotypeGVCFs -nt 1 -R $assembly -o $i.joint_bwa_mem_mdup_raw.vcf -V ../gcfFileList.list  -L $i" > lanceur_interval.sh
	echo "#$ -M a01346615@unet.univie.ac.at" > lanceur_interval.$i\.sh
	echo "#$ -m abe" >> lanceur_interval.$i\.sh
	echo "#$ -q max-1m.q" >> lanceur_interval.$i\.sh
	cd /sandbox/users/tleroy/Francesca/gvcf/tmp_vcf_$assemblyname >> lanceur_interval.$i\.sh
	echo "module load java" >> lanceur_interval.$i\.sh
	echo "time java -jar /sandbox/users/tleroy/Software/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -T GenotypeGVCFs -nt 1 -R $assembly -o /sandbox/users/tleroy/Francesca/gvcf/tmp_vcf_$assemblyname/$i.joint_bwa_mem_mdup_raw.vcf -V /sandbox/users/tleroy/Francesca/gvcf/gcfFileList.list  -L /sandbox/users/tleroy/Francesca/gvcf/tmp_vcf_$assemblyname/$i --includeNonVariantSites" >> lanceur_interval.$i\.sh
	qsub lanceur_interval.$i\.sh #> lanceur.$i & # qsub or sbatch for genotoul
done
# merging at the end
#for i in {1..$intervals}; do file=$(echo "scatter""$i"".intervals.joint_bwa_mem_mdup_raw.vcf"); if [ $i == 1 ]; then cp $file merged_joint_bwa_mem_mdup_raw.vcf; else grep -v "#" $file >> merged_joint_bwa_mem_mdup_raw.vcf; fi; done
