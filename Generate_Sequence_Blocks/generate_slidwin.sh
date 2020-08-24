# TL - 140318
# e.g. bash script_compute_pi_slidwin.sh /bigvol/Data/Taeniopygia_guttata/GCF_000151805.1_Taeniopygia_guttata-3.2.4_genomic_simplename.fna.scafflength Tguttata_fasta_files_scaffold_quantiles 100000 Tguttata_pi
scafflengthfile=$(echo "$1")
directory=$(echo "$2")
windowsize=$(echo "$3" | bc)
outprefix=$(echo "$4")
rm -r $outprefix.bed.tmp
mkdir $outprefix.bed.tmp
rm  $outprefix.$windowsize.bed
rm list_running_seq.txt

### generate a list of uniq individual
head -10000 merged_joint_bwa_mem_mdup_raw.filtered.vcf.qualsummary.ind | awk '{print $1}' | sort | uniq > merged_joint_bwa_mem_mdup_raw.filtered.vcf.IDs
while read line; do
    scaffold=$(echo "$line" | awk '{print $1}')
    lengthscaff=$(echo "$line" | awk '{print $2}')
    nbwindows=$(echo "($lengthscaff / $windowsize ) + 1" | bc)
    start=$(echo "1" )
    end=$(echo "$windowsize")
    for i in $(eval echo "{1..$nbwindows}"); do
       # generate a bed
        rm ./$outprefix.bed.tmp/$outprefix.$windowsize.tmp
        while read line; do # for each individual, print a line in a tmp bed file
            echo "$line.1 $start  $end" | sed 's/ \+/\t/g' >> ./$outprefix.bed.tmp/$outprefix.$windowsize.tmp
            echo "$line.2 $start  $end" | sed 's/ \+/\t/g' >> ./$outprefix.bed.tmp/$outprefix.$windowsize.tmp
        done < merged_joint_bwa_mem_mdup_raw.filtered.vcf.IDs
       # keep the information in a sumup bed file
        echo "$scaffold $start  $end" | sed 's/ \+/\t/g' >> $outprefix.$windowsize.bed
       # bedtools getfasta
        rm  ./$directory/$scaffold.fst.fai # rm to reinitiate the fai
        bedtools getfasta -fi ./$directory/$scaffold.fst -bed ./$outprefix.bed.tmp/$outprefix.$windowsize.tmp -fo ./$outprefix.bed.tmp/$outprefix.$scaffold.$start.$end.fasta
        echo "$outprefix.$scaffold.$start.$end.fasta" >> list_running_seq.txt
       # shift window
        previousend=$(echo "$end")
        start=$(echo "$end + 1" | bc)
        end=$(echo "$previousend + $windowsize" | bc)
    done
done < $scafflengthfile

