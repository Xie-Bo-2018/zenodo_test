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
${ref_NH1_fa} ${fastq_1} ${fastq_2} | samtools view -bS - > NH1.pe.bam
time samtools sort -@ 5 -m 4G NH1.pe.bam NH1; samtools index NH1.bam
time java -jar -Xmx4g -Djava.io.tmpdir=NH1/ MarkDuplicates.jar \
	INPUT=NH1.bam OUTPUT=NH1.dedup.bam \
	VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true METRICS_FILE=NH1.txt ASSUME_SORTED=true CREATE_INDEX=truerm 
rm NH1.bam
rm NH1.pe.bam
```
