# varcount

dependencies: htslib

to compile: `make`

to run: `./varcount <[vb]cf file> <[sb]am file>`

output: For each id in the vcf file, outputs a count for 

1) the number of reads overlapping the ref allele and 

2) the number of reads overlapping the alt allele.

Space separated.

`snp id ref_count alt_count`

