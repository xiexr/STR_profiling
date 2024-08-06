#filter snp
plink --bfile for127RNAseqout --exclude range geneposrange.txt --geno 0.1 --maf 0.1 --out filtersnp --recode
#transform format
plink --file filtersnp --out newfiltersnp --recode vcf-iid
#covariance
python3 covariance_between_samples.py