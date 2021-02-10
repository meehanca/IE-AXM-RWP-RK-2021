for i in *.trimmed; do ~/wheat/Spikes/programs/kallisto/kallisto quant -i \
~/wheat/Spikes/reference/index \
-o ../3_aligned/$i.txt --single -l 49 -s 0.0758229 -t 64 $i; done