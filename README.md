# varcount

## dependencies: 

* htslib
* c++11 compliant compiler.

## Compilation

`make`

## Description and Usage: 

```
Description:

Given a VCF and SAM file, calculate the alignment coverage over each ALT and
REF allele in the VCF. Outputs in VCF format to stdout.n\
Usage:

./varcount [options] <vcf> <sam>

<vcf>=STR               [bv]cf file name (required)
<sam>=STR               [bs]sam file name (required)
-s/--sample-name=STR    sample name in VCF output
                        (default: sample)
-g/--gt                 'predict' a genotype in the GT field based on threshold (-c/--threshold)
-c/--threshold=NUM      Requires -g/--gt. when NUM < 0: GT=1 if ALT >= 1.
                        when NUM >= 0: GT=1 if ALT-REF >= NUM; GT=0 if REF-ALT >= NUM; GT='.' otherwise.
                        (note: if NUM == 0 and REF == ALT, GT=0).
                        (default: 0)
-h/--help               print this help message
```
