#!/usr/bin/env bash

data_dir="../data"

ensembl_release="101"
dbsnp_relase="154"

cd $data_dir


### get reference genome - WGS | WES
curl ftp://ftp.ensembl.org/pub/release-$ensembl_release/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz


### download only coding sequences
curl ftp://ftp.ensembl.org/pub/release-$ensembl_release/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz > Homo_sapiens.GRCh38.cds.all.fa.gz

### transform fasta to bed and make coordinates 0-based 
gunzip -c Homo_sapiens.GRCh38.cds.all.fa.gz | grep ">" | awk '{print $3, $1, $2}' | sed 's/>//g' | sed -E 's/(chromosome|scaffold):GRCh38://g' | sed 's/:/\t/g' | awk '{$4=""; print $0}' | awk '{$2=($2 - 1); print $0}' | sed 's/\s/\t/g' > Homo_sapiens.GRCh38.cds.all.bed

#grep -v CHR Homo_sapiens.GRCh38.cds.all.bed > Homo_sapiens.GRCh38.cds.withoutSpecial.bed  ### vardict cant deal with weird CHR

gunzip -c Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | grep ">" | awk '{print $3, $1, $2}' | sed 's/>//g' | sed -E 's/(chromosome|scaffold):GRCh38://g' | sed 's/:/\t/g' | awk '{$4=""; print $0}' | awk '{$2=($2 - 1); print $0}' | sed 's/\s/\t/g' > Homo_sapiens.GRCh38.dna.primary_assembly.bed


### get common chromosomes between reference and coding seq bed and subset
comm -12 <(awk '{print $1}' Homo_sapiens.GRCh38.dna.primary_assembly.bed | sort | uniq) <(awk '{print $1}' Homo_sapiens.GRCh38.cds.all.bed | sort | uniq) > common_both.txt

### subset bed for only considered chromosomes - important for vardict + strelka
Rscript -e 'write.table(subset(read.table("Homo_sapiens.GRCh38.cds.all.bed"), V1 %in% read.table("common_both.txt")$V1), "Homo_sapiens.GRCh38.cds.all_filtered.bed", col.names=F, quote=F, row.names=F, sep="\t")'
rm common_both.txt





### short test genome
#curl ftp://ftp.ensembl.org/pub/release-$ensembl_release/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz > Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

