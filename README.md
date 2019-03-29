# varcount

## dependencies: 

* htslib
* c++11 compliant compiler.

## Compilation

`make`

## Description and Usage: 

```
varcount
========

description: given a V/BCF file and a S/BAM file, outputs a VCF with
    INFO/REFCOUNT and INFO/ALTCOUNT fields such that INFO/REFCOUNT is the count
    of alignments covering the ref allele and INFO/ALTCOUNT is the count of
    alignments covering the alt allele. 
    Optionally, output a genotype based on a given threshold for the difference
    between the alt and ref count

usage: ./varcount <vcf> <sam> <sample_name> <threshold>
        <vcf>: [bv]cf file name
        <sam>: [bs]sam file name
        <sample_name>: sample name in VCF output
        <threshold>: threshold of difference between ref and alt count to 'predict' genotype
            threshold>=0: ref/alt chosen if coverage exceeds the other allele
            by more than threshold. a tie defaults to ref
            threshold<0 : all variants and their counts are output, and genotypes are undefined

output: a VCF (to stdout) with INFO/REFCOUNT, INFO/ALTCOUNT, and\n\
    FORMAT/GT fields for the sample 'sample_name'\n");
```
