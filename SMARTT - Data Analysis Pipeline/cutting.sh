# usage: cutting.sh <name of paired file1 to be trimmed> 
#<name of paired file2 to be trimmed> <base name for output files>
sample_dir=$1
fastq_file1=$2
fastq_file2=$3
file_base=$4

#Adapter sequences listed below - input your own FWD_5 and REV3
#REV_5 and FWD_3 should be from the adapter

ADAPTER_FWD_5=^NNNNNNNNNNTATATGAGCATCGAAAGAA
ADAPTER_FWD_3=CTGTAGGCACCATCAATC
ADAPTER_REV_5=^NNNNNNNNNNGATTGATGGTGCCTACAG
ADAPTER_REV_3=TTCTTTCGATGCTCATATA

#Input your directory name here

output_dir=/

input1="$sample_dir"/"$fastq_file1"
input2="$sample_dir"/"$fastq_file2"
output1="$output_dir"/"$file_base"_trimmed1.fastq.gz
output2="$output_dir"/"$file_base"_trimmed2.fastq.gz
error_out="$output_dir"/"$file_base"_report

#this is for removing the 5' and 3' adaptors from the sample reads
# -q = quality score cutoff (trims 3' end sequence - 5 nt)
# -m = minimum length of the sequence AFTER removing adapters
# -O = min overlap cutoff (adaptor overlap)

cutadapt -m 45 -O 5 -a NNNNNNNNNNTATATGAGCATCGAAAGAA...CTGTAGGCACCATCAATC -A NNNNNNNNNNGATTGATGGTGCCTACAG...TTCTTTCGATGCTCATATA -o "$output1" -p "$output2" --pair-filter=any "$input1" "$input2"
