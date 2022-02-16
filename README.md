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
