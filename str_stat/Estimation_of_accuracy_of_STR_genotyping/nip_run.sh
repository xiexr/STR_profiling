
#!/bin/bash


#extract flanking sequence of STR sites
python3 nip_estimation

#extract repeat times of STRs sites
chromosomes=(1 2 3 4 5 6 7 8 9 10 11 12)

# Loop through each chromosome number
for chrom in "${chromosomes[@]}"
do
    python3 extract_rptimes.py "chrom${chrom}forTRF.txt" &
done