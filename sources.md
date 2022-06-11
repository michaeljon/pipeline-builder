# File and reference sources

## Reference assembly and feature files

These two file sets have to be pulled via a UI link right now. Navigate to the following and select the "Download Assembly" link in the upper-right corner. From there select RefSeq as the source database, then select "Genomic FASTA (.fna) for the FASTA file. When that starts downloading use the file type to select "Genomic GFF (.gff)" from under features and annotations.

Repeat this for both sets. You will get two `tar` files. Un-tar them, then navigate into the 'ncbi' directory to find the associated `.gz`. Those are the files you want.

- GRCh38.p13 - https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39/
- GRCh38.p14 - https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.40/

## Reference "Known sites" VCF

This replaces The Broad Institute's "dbsnp", "known_indels", and "Mills_and_1000G" indel sets with a single file that's updated on a regular basis. Match the reference assembly to the version. The patch is not going to be a direct match.

- https://ftp.ncbi.nih.gov/snp/latest_release/VCF/
