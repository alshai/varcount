#include <cstdio>
#include <vector>
#include <string>
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <getopt.h>
#include "hts_util.hpp"

struct VcntArgs {
    std::string vcf_fname = "";
    std::string sam_fname = "";
    std::string sample_name = "sample";
    int thres = -1;
    int gt = 0;
    int keep = 0;
    int verbose = 0;
    int diploid = 0;
    int alt_default = 0;
    int alt_sensitive = 0;
    int homs_only = 0;
    int min_ac = 0;
    int min_rc = 0;
    int min_c = 0;
    int gl = 0;
    float e = 0.01;
};

static inline std::tuple<std::string,std::string> truncate_str_pair(const std::string& s1, const std::string& s2) {
    return std::forward_as_tuple(s1.substr(s1.size()-1), s2.substr(s1.size()-1));
}


static inline bool var_match(const hts_util::Var& lv, const hts_util::Var& rv) {
    if (lv.type != rv.type) return false;
    switch (lv.type) {
        case hts_util::VTYPE::V_INS:
            return lv.pos+static_cast<int32_t>(lv.ref.size())-1 == rv.pos && truncate_str_pair(lv.ref, lv.alt) == truncate_str_pair(rv.ref, rv.alt) ;
        case hts_util::VTYPE::V_DEL:
            return lv.pos == rv.pos && truncate_str_pair(lv.alt, lv.ref) == truncate_str_pair(rv.alt, rv.ref);
        case hts_util::VTYPE::V_SNP:
            return lv.pos == rv.pos && lv.alt == rv.alt;
        default:
            fprintf(stderr, "no support for non SNPs & INDELs yet\n");
            return false;
    }
}


std::array<int32_t, 2> gt_by_threshold(const VcntArgs& args, const hts_util::Var& v) {
    std::array<int32_t, 2> gts;
    gts[0] = bcf_gt_missing;
    if (args.thres < 0 && v.ad[1]) {
        gts[0] = bcf_gt_phased(v.ad[1] != 0);
    } else if ((v.ad[1] || v.ad[0]) && std::abs(v.ad[0] - v.ad[1]) >= args.thres) {
        if (v.ad[0] < v.ad[1]) {
            gts[0] = bcf_gt_phased(1);
        } else if (v.ad[0] > v.ad[1])  {
            gts[0] = bcf_gt_phased(0);
        } else { // v.ad[0] == v.ad[1]
            gts[0] = bcf_gt_phased(args.alt_default);
        }
    } 
    gts[1] = gts[0];
    return gts;
}

std::array<int32_t, 2> gt_by_alt_evidence(const VcntArgs& args, const hts_util::Var& v) {
    std::array<int32_t, 2> gts;
    if (v.ad[1]) {
        gts[0] = bcf_gt_phased(1); gts[1] = bcf_gt_phased(1);
    } else if (v.ad[0]) {
        gts[0] = bcf_gt_phased(0); gts[1] = bcf_gt_phased(0);
    } else {
        gts[0] = bcf_gt_missing; gts[1] = bcf_gt_missing;
    }
    return gts;
}

std::array<int32_t, 2> gt_by_likelihood(const VcntArgs& args, const std::array<int32_t, 3>& pls) {
    std::array<int32_t, 2> gts;
    int min = pls[0], minidx = 0;
    if (pls[1] < min) { min = pls[1]; minidx = 1; }
    if (pls[2] < min) { min = pls[2]; minidx = 2; }
    switch (minidx) {
        case 0:
            gts[0] = bcf_gt_phased(0); gts[1] = bcf_gt_phased(0);
            break;
        case 1:
            gts[0] = bcf_gt_unphased(1); gts[1] = bcf_gt_unphased(0);
            break;
        case 2:
            gts[0] = bcf_gt_phased(1); gts[1] = bcf_gt_phased(1);
            break;
        default:
            fprintf(stderr, "something went wrong\n");
            exit(1);
            break;
    }
    return gts;
}

void varcount(const VcntArgs& args) {
    samFile* sam_fp = sam_open(args.sam_fname.data(), "r");
    bam_hdr_t* sam_hdr = sam_hdr_read(sam_fp);
    bam1_t* aln = bam_init1();

    vcfFile* vcf_fp = bcf_open(args.vcf_fname.data(), "r");
    bcf_hdr_t* vcf_hdr = bcf_hdr_read(vcf_fp);
    if (bcf_hdr_set_samples(vcf_hdr, NULL, 0)) {
        fprintf(stderr, "error setting samples\n");
        exit(1);
    } 

    hts_util::contig2map_map<hts_util::Var> contig2vars(hts_util::bcf_to_map(vcf_fp, vcf_hdr));

    int32_t pid = -1;
    bam1_core_t* c = nullptr;
    hts_util::pos2var_map<hts_util::Var>* vmap = nullptr;
    while (sam_read1(sam_fp, sam_hdr, aln) >= 0) {
        c = &aln->core;
        // we only care about uniquely mapped reads
        if((c->flag & BAM_FUNMAP) == 0 && (c->flag & BAM_FSECONDARY) == 0) {
            if (c->tid != pid) {
                // load new hash table
                // TODO: maybe we should do the VCF reading step here (subsetting by region) instead of reading the whole thing beforehand
                vmap = &(contig2vars[sam_hdr->target_name[c->tid]]);
            }
            pid = c->tid;

            std::vector<hts_util::Var> aln_vars(hts_util::Var::from_bam(aln));
            if (args.verbose) {
                fprintf(stderr, "a %s ", bam_get_qname(aln));
                // hts_util::print_varlist(aln_vars, stderr);
            }
            // can we vectorize or at least parallelize this? would it be worth it?
            // as of now, the performance bottleneck is the VCF reading
            for (int32_t i = c->pos; i <= bam_endpos(aln); ++i) {
                auto found = vmap->find(i);
                if (found != vmap->end()) {
                    // hts_util::Var* v_cached = &found->second[0];
                    for (auto&& vv: found->second) {
                        /*
                        if (vv.rec_start) // keep pointer to this record
                            v_cached = &vv;
                            */
                        bool x = 0;
                        for (const auto& av: aln_vars) {
                            x |= var_match(vv, av);
                        }
                        // x ? ++(v_cached->ad[1]) : ++(v_cached->ad[0]);
                        x ? ++(vv.ad[1]) : ++(vv.ad[0]);
                    }
                    if (args.verbose) {
                        fprintf(stderr, "v %s ", bam_get_qname(aln));
                        hts_util::print_varlist(found->second, stderr);
                    }
                }
            }
        }
    }

    bam_destroy1(aln);
    bam_hdr_destroy(sam_hdr);
    sam_close(sam_fp);


    // output in VCF format!!
    vcfFile* out_vcf_fp = bcf_open("-", "w");
    bcf_hdr_t* out_vcf_hdr = bcf_hdr_dup(vcf_hdr);

    bcf_hdr_add_sample(out_vcf_hdr, args.sample_name.data());
    bcf_hdr_remove(out_vcf_hdr, BCF_HL_INFO, NULL);
    bcf_hdr_remove(out_vcf_hdr, BCF_HL_FMT, NULL);
    // cleanup from original VCF
    bcf_hdr_remove(out_vcf_hdr, BCF_HL_GEN, "fileDate");
    bcf_hdr_remove(out_vcf_hdr, BCF_HL_GEN, "bcftools_viewVersion");
    bcf_hdr_remove(out_vcf_hdr, BCF_HL_GEN, "bcftools_viewCommand");
    bcf_hdr_remove(out_vcf_hdr, BCF_HL_STR, NULL);
    bcf_hdr_set_version(out_vcf_hdr, "VCFv4.3");
    // bcf_hdr_append(out_vcf_hdr, "##INFO=<ID=ALTCNT,Number=1,Type=Integer,Description=\"Count of reads covering alt allele\">");
    // bcf_hdr_append(out_vcf_hdr, "##INFO=<ID=REFCNT,Number=1,Type=Integer,Description=\"Count of reads covering ref allele\">");
    bcf_hdr_append(out_vcf_hdr, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic Depth\">");
    if (args.gt)
        bcf_hdr_append(out_vcf_hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    if (args.gl) 
        bcf_hdr_append(out_vcf_hdr, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">");
    if (bcf_hdr_write(out_vcf_fp, out_vcf_hdr)) {
        fprintf(stderr, "bcf_hdr_write error\n");
        exit(1);
    }

    bcf1_t* out_vcf_rec = bcf_init();
    bcf_clear(out_vcf_rec);

    // TODO: sort the output?
    for (const auto& vs: contig2vars) { // vs: {seq, map}
        for (const auto& v: vs.second)  { // v: {pos, Varlist}
            auto it = v.second.begin();
            while (it != v.second.end()) {
                /* TODO: print each allele as its own thing */
                std::string allele_str(it->ref + "," + it->alt);
                // auto nit = std::next(it);
                // while (nit != v.second.end() && !nit->rec_start) {
                //     allele_str += "," + nit->alt;
                //     ++nit;
                // }
                out_vcf_rec->rid = bcf_hdr_name2id(out_vcf_hdr, vs.first.data());
                out_vcf_rec->pos = it->pos;
                out_vcf_rec->rlen = it->ref.size();
                bcf_update_alleles_str(out_vcf_hdr, out_vcf_rec, allele_str.data());
                bcf_update_id(out_vcf_hdr, out_vcf_rec, it->id.data());
                bcf_update_info_int32(out_vcf_hdr, out_vcf_rec, "ALTCNT", &(it->ad[1]), 1);
                bcf_update_info_int32(out_vcf_hdr, out_vcf_rec, "REFCNT", &(it->ad[0]), 1);
                std::array<int32_t, 3> pls;
                std::array<int32_t, 2> gts;
                if (args.gl || (args.gt && args.thres <= 0)) {
                    /* TODO: handle indels smarter. Use individual base qualities */
                    pls = hts_util::get_pls_naive_normalized(it->ad[0], it->ad[1], args.e);
                }
                int gt_pass = 1;
                if (args.gt) { // naive genotyping
                    if (args.alt_sensitive) {
                        gts = gt_by_alt_evidence(args, *it);
                    } else if (args.thres >= 0) {
                        gts = gt_by_threshold(args, *it);
                    } else { // default: genotype by likelihood
                        gts = gt_by_likelihood(args, pls);
                    }
                    bcf_update_genotypes(out_vcf_hdr, out_vcf_rec, gts.data(), 2);
                    gt_pass = !bcf_gt_is_missing(gts[0]);
                    gt_pass &= (!args.homs_only || (args.homs_only && bcf_gt_is_phased(gts[0])));
                }
                if (args.keep || ( (it->ad[0] + it->ad[1] >= args.min_c && it->ad[1] >= args.min_ac && it->ad[0] >= args.min_rc) && (!args.gt || (args.gt && gt_pass)))) {
                    if (args.gl) bcf_update_format_int32(out_vcf_hdr, out_vcf_rec, "PL", pls.data(), 3);
                    bcf_update_format_int32(out_vcf_hdr, out_vcf_rec, "AD", it->ad.data(), 2);
                    if (bcf_write(out_vcf_fp, out_vcf_hdr, out_vcf_rec)) {
                        fprintf(stderr, "bcf_write error\n");
                        exit(1);
                    }
                }
                bcf_clear(out_vcf_rec);

                // it = nit; // DELETED
                ++it;
            }
        }
    }
    bcf_hdr_destroy(out_vcf_hdr);
    bcf_close(out_vcf_fp);

    bcf_hdr_destroy(vcf_hdr);
    bcf_close(vcf_fp);
}

void print_help() {
fprintf(stderr,
"Description: \n\
\n\
Given a VCF and SAM file, calculate the alignment coverage over each ALT and\n\
REF allele in the VCF. Outputs in VCF format to stdout.\n\n\
Usage:\n\
\n\
./varcount [options] <vcf> <sam>\n\
\n\
<vcf>=STR               [bv]cf file name (required)\n\
<sam>=STR               [bs]sam file name (required)\n\
-s/--sample-name=STR    sample name in VCF output\n\
                        (default: sample)\n\
-g/--gt                 'predict' a genotype in the GT field based on threshold (-c/--threshold)\n\
-c/--threshold=NUM      Requires -g/--gt. when NUM < 0: GT=1 if ALT >= 1.\n\
                        when NUM >= 0: GT=1 if ALT-REF >= NUM; GT=0 if REF-ALT >= NUM; GT='.' otherwise.\n\
                        (note: if NUM == 0 and REF == ALT, GT=0).\n\
                        (default: 0)\n\
-k/--keep               if threshold not met, print record anyway with undefined genotype\n\
-v/--verbose            prints detailed logging information to stderr\n\
-h/--help               print this help message\n\
\n");
}

int main(int argc, char** argv) {
    VcntArgs args;
    static struct option long_options[] {
        {"sample-name", required_argument, 0, 's'},
        {"threshold", required_argument, 0, 'c'},
        {"genotype", no_argument, &args.gt, 1},
        {"keep", no_argument, &args.keep, 1},
        {"diploid", no_argument, &args.diploid, 'd'},
        {"alt-default", no_argument, &args.alt_default, 1},
        {"alt-sensitive", no_argument, &args.alt_sensitive, 1},
        {"homs-only", no_argument, &args.homs_only, 1},
        {"min-alt-count", required_argument, 0, 'a'},
        {"min-ref-count", required_argument, 0, 'r'},
        {"min-total-count", required_argument, 0, 'm'},
        {"help", no_argument, 0, 'h'},
        {"verbose", no_argument, &args.verbose, 1},
        {0,0,0,0}
    };

    int ch;
    int argpos = 0;
    while ( (ch = getopt_long(argc, argv, "-:s:m:a:r:gvkhl", long_options, NULL)) != -1 ) {
        switch(ch) {
            case 0:
                break;
            case 1:
                if (argpos == 0) args.vcf_fname = std::string(optarg);
                else if (argpos == 1) args.sam_fname = std::string(optarg);
                else fprintf(stderr, "ignoring argument %s\n", optarg);
                ++argpos;
                break;
            case 2:
                break;
            case 's':
                args.sample_name = std::string(optarg);
                break;
            case 'c':
                args.thres = std::atoi(optarg);
                break;
            case 'a':
                args.min_ac = std::atoi(optarg);
                break;
            case 'r':
                args.min_rc = std::atoi(optarg);
                break;
            case 'm':
                args.min_c = std::atoi(optarg);
                break;
            case 'd':
                args.diploid = true;
                break;
            case 'g':
                args.gt = 1;
                break;
            case 'k':
                args.keep = 1;
                break;
            case 'v':
                args.verbose = 1;
                break;
            case 'l':
                args.gl = 1;
                break;
            case 'e':
                args.e = std::atof(optarg);
                break;
            case 'h':
                print_help();
                exit(0);
                break;
            case '?':
                print_help();
                fprintf(stderr, "error: unknown option -%c\n", optopt);
                exit(1);
                break;
            default:
                print_help();
                exit(1);
        }
    }

    if (args.vcf_fname == "" || args.sam_fname == "") {
        print_help();
        fprintf(stderr, "error: vcf and sam are mandatory arguments\n");
        exit(1);
    }

    varcount(args);
}
