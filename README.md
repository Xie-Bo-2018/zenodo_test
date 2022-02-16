# SnpEff (Building a database from GFF files)
In this step, you can refer to the documentation of SnpEff software. (https://pcingola.github.io/SnpEff/se_buildingdb/#option-2-building-a-database-from-gff-files)

## Add GRCh38 ALT sequence to the configuration file

1. Get a GFF file (into path/to/snpEff/data/GRCh38_ALT_PRSS/genes.gff.gz):
```
mkdir path/to/snpEff/data/GRCh38_ALT_PRSS
mv genes.gff.gz path/to/snpEff/data/GRCh38_ALT_PRSS/
```

2. Add the sequence to the config file:
```
mv sequences.fa path/to/snpEff/data/GRCh38_ALT_PRSS/
```

3. Edit the config file to create the new genome:
```
vim snpEff.config
```
  Add the following lines (you are editing snpEff.config):
```
# GRCh38_ALT_PRSS, GENCODE_v39 GRCh38_ALT_PRSS
GRCh38_ALT_PRSS.genome : path/to/snpEff/data//GRCh38_ALT_PRSS/sequences.fa
```

5. Create database (note the "-gff3" flag):
```
cd /path/to/snpEff
java -jar snpEff.jar build -gff3 -v GRCh38_ALT_PRSS
```


##########EAGLE test##########
```
ref_HX1_fa=hx1f4s4_3rdfixedv2.fa
ref_HX1_fai=hx1f4s4_3rdfixedv2.fa.fai
ref_NH1_fa=GWHAAAS00000000.genome.fasta
ref_HX1_filtered_fa=hx1f4s4_3rdfixedv2.filtered.fa

python Filter_HX1_fasta_contig_length.py ${ref_HX1_fa} ${ref_HX1_fai} ${ref_HX1_filtered_fa}

configureEAGLE.pl \
  --run-info=RunInfo_PairedReads2x251Cycles2x64Tiles.xml  \
  --reference-genome=${ref_HX1_filtered_fa} \
  --coverage-depth=30 \
  --motif-quality-drop-table=MotifQualityDropTables/DefaultMotifQualityDropTable.tsv \
  --template-length-table=TemplateLengthTables/TemplateLengthTableFrom2x250Run.tsv \
  --quality-table=NewQualityTable.read1.length251.qtable2 \
  --quality-table=NewQualityTable.read2.length251.qtable2
cd EAGLE
make fastq -j 15
```

##########001-a-raw reads mapping##########
```
fastq_1=EAGLE_S1_L001_R1_001.fastq.gz
fastq_2=EAGLE_S1_L001_R2_001.fastq.gz
time bwa mem -M -t 10 -R "@RG\tID:HX1\tSM:HX1\tLB:HX1\tPU:HX1\tPL:ILLUMINA" \
   ${ref_NH1_fa} ${fastq_1} ${fastq_2}  | samtools view -bS - > HX1.pe.bam
time samtools sort -@ 5 -m 4G HX1.refNH1.pe.bam NH1; samtools index HX1.bam
time java -jar -Xmx4g -Djava.io.tmpdir=HX1/ MarkDuplicates.jar \
	INPUT=HX1.bam OUTPUT=HX1.dedup.bam \
	VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true METRICS_FILE=NH1.txt ASSUME_SORTED=true CREATE_INDEX=true
rm HX1.bam;rm HX1.pe.bam
```

##########001-b-Use GATK software to do variants calling##########
```
ref_fai_path=${ref_NH1_fai}
ref_chr_list=`cat $ref_fai_path | awk '{print $1}'`
for chr in $ref_chr_list
    do
        time gatk --java-options \"-Xmx3G -XX:ParallelGCThreads=2 -Dsamjdk.compression_level=5 \" HaplotypeCaller \
        -R ${ref_NH1_fa} -ploidy 1 -L $chr -I HX1.dedup.bam -O HX1.$chr.g.vcf.gz -ERC GVCF -G StandardAnnotation \
        -G AS_StandardAnnotation -G StandardHCAnnotation --seconds-between-progress-updates 30
    done
sh Combine_list.sh;sh 170.JointCalling.sh
for chr in $ref_chr_list
    do
        sh 170.${chr}.sh
    done
sh 170.combine.sh
```

##########001-c-GATK hard-filtering##########
```
time gatk SelectVariants -select-type SNP -V HX1.genomewide.hc.vcf.gz -O HX1.genomewide.hc.snp.vcf.gz
time gatk VariantFiltration -V HX1.genomewide.hc.snp.vcf.gz \
    --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0" \
    --filter-name "Filter" -O HX1.genomewide.hc.snp.filter.vcf.gz
time gatk SelectVariants -select-type INDEL -V HX1.genomewide.hc.vcf.gz -O HX1.genomewide.hc.indel.vcf.gz
time gatk VariantFiltration -V HX1.genomewide.hc.indel.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0" \
    --filter-name "Filter" -O HX1.genomewide.hc.indel.filter.vcf.gz
time gatk MergeVcfs -I HX1.genomewide.hc.snp.filter.vcf.gz -I HX1.genomewide.hc.indel.filter.vcf.gz \
    -O HX1.genomewide.hc.filter.vcf.gz
vcftools --gzvcf HX1.genomewide.hc.filter.vcf.gz --remove-filtered-all \
--recode -c |bgzip -c > HX1.genomewide.hc.filtered.vcf.gz
```

##########rtg-tools##########
```
paftools_wgs_call.sh ${ref_NH1_fa} ${ref_HX1_fa};
bgzip NH1.HX1.vcf;tabix NH1.HX1.vcf.gz
base_vcf=NH1.HX1.vcf.gz
query_vcf=HX1.genomewide.hc.filtered.vcf.gz
rtg format -o NH1.sdf ${ref_NH1_fa}
rtg vcfeval -b base_vcf -c query_vcf -o output -t NH1.sdf
```


##########liftoff and snpeff##########
```
gff3_db=gencode.v34.annotation.gff3_db
ref_NH1_fa=GWHAAAS00000000.genome.fasta
ref_GRCh38_fa=GRCh38_full_analysis_set_plus_decoy_hla.fa
HX1_vcf=HX1.genomewide.hc.filtered.vcf.gz
HX1_medically_gene_vcf=HX1.genomewide.hc.medically_gene.filtered.vcf.gz
HX1_medically_gene_snpEff_vcf=HX1.genomewide.hc.medically_gene.snpEff_ann.filtered.vcf.gz
liftoff -db ${gff3_db} -a 0.9 -s 0.9 -exclude_partial -p 10 \
-o NH1.gencode.v34.gff -u NH1_unmapped.txt ${ref_NH1_fa} ${ref_GRCh38_fa}
python get_NH1_medically_genes.py NH1.gencode.v34.gff
bcftools view -R NH1.medically_gene.bed ${HX1_vcf} |bgzip > ${HX1_medically_gene_vcf} &
python snpEff_config.py /path/to/snpEff/ /path/to/data/
java -Xmx4g -jar snpEff.jar -v NH1 ${HX1_medically_gene_vcf} -c /path/to/snpEff.config \
-datadir /path/to/data/ | bgzip 1> ${HX1_medically_gene_snpEff_vcf} 
```
