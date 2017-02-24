# Data in this directory

### Alt alignments (dir: alt_alignments)
Alt alignments are collected from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_structure/).
The alignmentfiles have been further processed into a file on the format
`main_chr_start_position, main_chr_end_position, alt_locus_start_position, alt_locus_end_position, cigar string`

All coordinates are 0-based with exlusive end position.

### Genes (dir: genes)
Genes are collected from the UCSC Table Browser (http://genome.ucsc.edu/cgi-bin/hgTables?command=start).

For refseq genes the track _Refseq genes_ was used. For Gencode genes, the track _ALL GENCODE V24_  was used.
This tracks were processed into one file per chromosome.

### grch38.chrom.sizes
Downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes.

### Position of alternative loci (grch38_alt_loci.txt)
Downloaded from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.7_GRCh37.p6/GCA_000001405.7_GRCh37.p6_assembly_regions.txt

grch38_alt_loci.txt is a preprocessed version of that file, containing the following fields:
alt_locus_id main_chr alt_locus_start alt_locus_end alt_locus_length
