#!/bin/bash

# Files
aedes_faa=data/ref/all_faa/aedes.faa
dmel_faa=data/ref/all_faa/dmel.faa
aedes_translens=results/ortho/annot/aedes_translens.tsv
dmel_translens=results/ortho/annot/dmel_translens.tsv
aedes_longtrans=results/ortho/annot/aedes_longtrans.faa
dmel_longtrans=results/ortho/annot/dmel_longtrans.faa
dmel_db=results/ortho/diamond/dmel_db/dmel_longtrans.dmnd
aedes_db=results/ortho/diamond/aedes_db/aedes_longtrans.dmnd

# Copy proteome files to a single dir
mkdir -p data/ref/all_faa
ln -sv "$PWD"/data/ref/dmel/dmel-all-translation-r6.54.fasta "$dmel_faa"
ln -sv "$PWD"/data/ref/aedes/VectorBase-49_AaegyptiLVP_AGWG_AnnotatedProteins.fasta "$aedes_faa"

# Get longest transcript per gene for aedes
grep ">" "$aedes_faa" |
    sed -E 's/>([^ ]+).*gene=([^ ]+).*protein_length=([^ ]+).*/\1\t\2\t\3/' |
    sort -k2,2 -k3nr > "$aedes_translens"
seqkit grep -f <(sort -u -k2,2 "$aedes_translens" | cut -f1) "$aedes_faa" > "$aedes_longtrans"

# Get longest transcript per gene for dmel
grep ">" "$dmel_faa" |
    sed -E 's/>([^ ]+).*parent=([^,]+).*length=([^;]+).*/\1\t\2\t\3/' |
    sort -k2,2 -k3nr > "$dmel_translens"
seqkit grep -f <(sort -u -k2,2 "$dmel_translens" | cut -f1) "$dmel_faa" > "$dmel_longtrans"

# Run DIAMOND
sbatch mcic-scripts/blast/diamond_db.sh -i "$dmel_longtrans" -o results/ortho/diamond/dmel_db
sbatch mcic-scripts/blast/diamond_db.sh -i "$aedes_longtrans" -o results/ortho/diamond/aedes_db

sbatch mcic-scripts/blast/diamond.sh -i "$aedes_longtrans" --db "$dmel_db" \
    -o results/ortho/diamond/to_aedes --pct_id 30 --max_target_seqs 1 --sens ultra-sensitive
sbatch mcic-scripts/blast/diamond.sh -i "$dmel_longtrans" --db "$aedes_db" \
    -o results/ortho/diamond/to_dmel --pct_id 30 --max_target_seqs 1 --sens ultra-sensitive

# Get reciprocal best blast hits
todmel=results/ortho/diamond/to_aedes/diamond_out.tsv
toaedes=results/ortho/diamond/to_dmel/diamond_out.tsv
join -t $'\t' -1 1 -2 2 \
    <(tail -n+4 "$todmel" | cut -f 1,2 | sort -k1,1) \
    <(tail -n+4 "$toaedes" | cut -f 1,2 | sort -k2,2) |
    awk '$2 == $3' | cut -f1,2  > results/ortho/diamond/rbbh.tsv
wc -l results/ortho/diamond/rbbh.tsv #> 7767


# ==============================================================================
#                                   TESTS
# ==============================================================================
# Run Orthofinder for comparison
sbatch mcic-scripts/compgenom/orthofinder.sh -i results/ortho/annot -o results/ortho/orthofinder
#> OrthoFinder assigned 24849 genes (86.6% of total) to 9002 orthogroups.
#> Fifty percent of all genes were in orthogroups with 2 or more genes (G50 was 2) and were contained in the largest 3754 orthogroups (O50 was 3754).
#> There were 7869 orthogroups with all species present and 6159 of these consisted entirely of single-copy genes.

#awk "/^>/ {n++} n>5{exit} {print}" "$aedes_longtrans" > aedes_test.faa
#sbatch mcic-scripts/blast/diamond.sh -i aedes_test.faa --db "$dmel_db" -o DIAMOND_TEST --pct_id 30 --sens fast

# Run JustOrthologs
# See scripts/justorthologs.sh

# Check OrthoDB - https://data.orthodb.org/download/README.txt
#aedes_code=$(zgrep -P "Aedes aegypti\t" odb11v0_species.tab.gz | cut -f2) #7159_0
#zgrep -P "\t${aedes_code}\t" odb11v0_genes.tab.gz | head
#7159_0:000000   7159_0  XP_021693426.1  5568916 A0A6I8TTE2              5568916 5568916
#https://www.ncbi.nlm.nih.gov/protein/XP_021693426.1/
