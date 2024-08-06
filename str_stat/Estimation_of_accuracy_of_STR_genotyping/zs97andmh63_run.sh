#extract sequence for blast
python3 extract_seq_for_blast.py

#blast the flanking sequence of nip to zs97 genome
blastn -out ${method}_${zs97ormh63}_100bp.txt -outfmt 6 -query seq_of_{$method}_100bp.txt -db MH63RS3_db -evalue 100 -strand both

#filter blast
python3 blast_filter.py

#mh63orzs97 extract flanking sequence of STR sites
python3 mh63_and_zs97_estimation.py

#extract repeat times
chromosomes=(1 2 3 4 5 6 7 8 9 10 11 12)

# Loop through each chromosome number
for chrom in "${chromosomes[@]}"
do
    python3 extract_rptimes.py "chrom${chrom}forTRF.txt" &
done