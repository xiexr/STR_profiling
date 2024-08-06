for SAMPLE in "${SAMPLES[@]}"; do
    echo "Generating GVCF for $SAMPLE..."
    $GATK --java-options "-Xmx4g" HaplotypeCaller \
        -R $REF \
        -I ${OUTPUT_DIR}/${SAMPLE}_sorted.bam \
        -O ${OUTPUT_DIR}/${SAMPLE}.g.vcf.gz \
        --emitRefConfidence GVCF \
        --variant_index_type LINEAR \
        --variant_index_parameter 128000 \
        --min-base-quality-score 20
done

# Step 3: Combine GVCFs into a single VCF
echo "Combining GVCFs..."
$GATK --java-options "-Xmx4g" CombineGVCFs \
    -R $REF \
    -O ${OUTPUT_DIR}/combined.g.vcf.gz \
    ${OUTPUT_DIR}/*.g.vcf.gz

# Step 4: Genotype the combined GVCF
echo "Genotyping combined GVCF..."
$GATK --java-options "-Xmx4g" GenotypeGVCFs \
    -R $REF \
    -V ${OUTPUT_DIR}/combined.g.vcf.gz \
    -O ${OUTPUT_DIR}/final_variants.vcf

# Step 5: Filter the variants
echo "Filtering variants..."
$GATK --java-options "-Xmx4g" VariantFiltration \
    -R $REF \
    -V ${OUTPUT_DIR}/final_variants.vcf \
    -O ${OUTPUT_DIR}/filtered_variants.vcf \
    --filter-expression "AC < 10 || DP < 50" \
    --filter-name "LowQual"

echo "Pipeline complete!"