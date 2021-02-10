# Quality control and trimming with Fastp
for i in *.fastq; do
  fastp --thread 16 --length_required 40 --cut_right --cut_window_size 4 --cut_mean_quality 20 --overrepresentation_analysis \
  -i $i.fastq \
  -o ../2_trimmed/$i.fastq _R1_trimmed.fastq.gz;
done