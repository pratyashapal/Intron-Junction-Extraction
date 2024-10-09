This Python script processes RNA-seq alignments in SAM format to identify intron-exon junctions. It extracts split reads that span introns, counts how many times each junction is supported, and matches these junctions to gene boundaries provided in a separate gene summary file.

Features
Extract Split Reads: Detects split reads from SAM files using CIGAR strings.
Junction Counting: For each unique intron-exon junction, counts the supporting reads.
Gene Boundary Mapping: Matches junctions to genes based on given genomic coordinates.
Customizable Input: Accepts any SAM file and gene summary file as input.
Input Files
SAM File: Contains RNA-seq alignments. Key columns used:

Chromosome name (RNAME)
Alignment start position (POS)
CIGAR string for alignment description (CIGAR)
Read alignment count (NH:i
)
Gene Summary File: Tab-separated file with:

Gene ID
Transcript ID
Gene location in the format: Chromosome:start..end(strand)
Output
The script generates a tab-separated output file listing gene IDs, intron junctions (start/end positions), and the number of supporting reads for each junction.
