#!/bin/bash

# ==============================================================================
#                       PREP FOR ORTHOLOG ANALYSIS
# ==============================================================================
mkdir -p data/ref/aeeg data/ref/dmel

# Get the Aedes proteome
aedes_proteome_url=https://vectorbase.org/common/downloads/release-49/AaegyptiLVP_AGWG/fasta/data/VectorBase-49_AaegyptiLVP_AGWG_AnnotatedProteins.fasta
wget -P data/ref/aeeg "$aedes_proteome_url"

# Get the Drosophila proteome, GFF and GAF
#? https://wiki.flybase.org/wiki/FlyBase:Downloads_Overview#Gene_Association_File_-_GAF_.28gene_association.fb.gz.29
#? https://wiki.flybase.org/wiki/FlyBase:Gene_Ontology_(GO)_Annotation
wget -P data/ref/dmel ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.54_FB2023_05/fasta/dmel-all-translation-r6.54.fasta.gz
wget -P data/ref/dmel ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.54_FB2023_05/gff/dmel-all-r6.54.gff.gz
wget -P data/ref/dmel ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.54_FB2023_05/fasta/dmel-all-chromosome-r6.54.fasta.gz
wget -P data/ref/dmel ftp://ftp.flybase.net/releases/FB2023_05/precomputed_files/go/gene_association.fb.gz
gunzip data/ref/dmel/*gz

# Fix the Dmel GFF, which contains a weird first line, and has the DNA sequences at the end of the file
sed -i 's/.*##gff-version 3/##gff-version 3/' data/ref/dmel/dmel-all-r6.54.gff
sed '/##FASTA/Q' data/ref/dmel/dmel-all-r6.54.gff > data/ref/dmel/dmel-all-r6.54_nofasta.gff
sed -n '/##FASTA/,$p' data/ref/dmel/dmel-all-r6.54.gff | tail -n+2 > data/ref/dmel/dmel-all-r6.54_fromgff.fa

# Check the nr of genes in Aedes
grep -v "##" data/ref/aedes/VectorBase-49_AaegyptiLVP_AGWG.gff | grep -Pc "\tgene\t" # 19,804 genes
