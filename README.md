# varcount

dependencies: htslib

to compile: `make`

to run: `./varcount <[vb]cf file> <[sbam file]>

output: For each id in the vcf file, outputs a count for the number of reads
overlapping the ref allele and the number of reads overlapping the alt allele

`snp id ref_count alt_count`

