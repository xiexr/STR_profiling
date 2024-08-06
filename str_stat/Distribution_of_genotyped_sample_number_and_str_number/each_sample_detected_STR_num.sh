for i in {1..12}
do
    nohup python3 each_sample_detected_STR_num.py Lastresultchr${i}.txt lobstrchrom${i}.vcf hipstrchrom${i}.vcf gatkchr$(printf "%02d" $i).csv &
done