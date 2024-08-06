#depth
bamdst -p $Nipponbare_each_chrom_length.bed  $sample.bam -o $output_dir
#call rate
python3 call_rate_calc.py