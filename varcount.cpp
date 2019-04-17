#include <cstdio>
#include <vector>
#include <list>
#include <string>
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <getopt.h>
#include "mdparse.hpp"
#include "flat_hash_map.hpp"
#include "varcount.hpp"


vcnt::VarList vcnt::bam_to_vars(bam1_t* aln) {
    VarList vs;
    // look at md string and gather dels, snps
    char* md;
    std::list<std::pair<std::string, int32_t>> snps;
    std::list<std::pair<std::string, int32_t>> dels;
    if ((md = bam_aux2Z(bam_aux_get(aln, "MD")))) {
        std::vector<MDPos> mds = md_parse(md);
        for (auto m: mds) {
            if (m.st == MD_DEL) {
                std::string ref = "";
                ref += m.str;
                dels.push_back(std::pair<std::string, int32_t>(ref, aln->core.pos + m.p - 1));
            }
            else if (m.st == MD_SNP) {
                snps.push_back(std::pair<std::string, int32_t>(m.str, aln->core.pos + m.p));
            }
        }
    }

    // walk through the CIGAR string
    uint32_t* cs = bam_get_cigar(aln);
    int32_t qpos = 0;
    int32_t rpos = aln->core.pos;
    int32_t rlen, qlen;
    for (uint32_t i = 0; i < aln->core.n_cigar; ++i) {
        rlen = (bam_cigar_type(cs[i]) & 2) ? bam_cigar_oplen(cs[i]) : 0;
        qlen = (bam_cigar_type(cs[i]) & 1) ? bam_cigar_oplen(cs[i]) : 0;
        if (bam_cigar_op(cs[i]) == BAM_CINS) {
            std::string ins_seq = "";
            ins_seq += seq_nt16_str[bam_seqi(bam_get_seq(aln), qpos-1)]; // adding the previous character here for compatibility with VCF format
            for (uint32_t j = 0; j < bam_cigar_oplen(cs[i]); ++j) {
                ins_seq += seq_nt16_str[bam_seqi(bam_get_seq(aln), qpos+j)];
            }
            vs.push_back(Var(rpos-1, VTYPE::V_INS, ins_seq, ins_seq.substr(0,1)));
        } else if (bam_cigar_op(cs[i]) == BAM_CMATCH) {
            // check potential snps here
            for (auto it = snps.begin(); it != snps.end(); ) {
                int32_t s = it->second;
                if (s >= rpos && s < rpos + rlen) {
                    std::string snp = "";
                    snp += seq_nt16_str[bam_seqi(bam_get_seq(aln), qpos + (s - rpos))];
                    vs.push_back(Var(s, VTYPE::V_SNP, snp, it->first));
                    it = snps.erase(it); // we do this so we don't recheck snps that have already been determined, but... maybe deleting it would actually be more expensive?
                } else ++it;
            }
        } else if (bam_cigar_op(cs[i]) == BAM_CDEL) {
            for (auto it = dels.begin(); it != dels.end(); ) {
                int32_t d = it->second;
                if (d == rpos - 1) {
                    std::string alt = "";
                    alt += seq_nt16_str[bam_seqi(bam_get_seq(aln), qpos - 1)];
                    vs.push_back(Var(d, VTYPE::V_DEL, alt, alt + it->first));
                    it = dels.erase(it);
                } else ++it;
            }
        }
        rpos += rlen;
        qpos += qlen;
    }
    return vs;
}

// might encounter bug if STRLEN(REF) > 1 && STRLEN(ALT) > 1
vcnt::VarList vcnt::bcf_to_vars(bcf1_t* b) {
    VarList vs;
    char* ref = b->d.allele[0];
    for (uint32_t i = 1; i < b->n_allele; ++i) {
        char* alt = b->d.allele[i];
        if (alt[0] == '.') continue;
        if (strlen(alt)  < strlen(ref)) { // DEL
            vs.push_back(Var(b->pos, VTYPE::V_DEL, alt, ref, b->d.id)); // don't need alt here
        } else if (strlen(alt) > strlen(ref)) { // INS
            vs.push_back(Var(b->pos, vcnt::VTYPE::V_INS, alt, ref, b->d.id)); // don't need ref here
        } else { // SNP
            vs.push_back(Var(b->pos, VTYPE::V_SNP, alt, ref, b->d.id)); // don't need ref here
        }
    }
    vs[0].rec_start = 1;
    return vs;
}

void print_varlist(vcnt::VarList vs, FILE* out) {
    for (const auto& v: vs) {
        fprintf(out, "(%d %d %s %s", v.pos, static_cast<int>(v.type), v.alt.data(), v.ref.data());
        if (v.id.size()) fprintf(out, " %s", v.id.data());
        if (v.rc + v.ac) fprintf(out, " %d %d", v.rc, v.ac);
        fprintf(out, ") ");
    } fprintf(out, "\n");
}

void vcnt::varcount(const VcntArgs& args) {
    samFile* sam_fp = sam_open(args.sam_fname.data(), "r");
    bam_hdr_t* sam_hdr = sam_hdr_read(sam_fp);
    bam1_t* aln = bam_init1();

    vcfFile* vcf_fp = bcf_open(args.vcf_fname.data(), "r");
    bcf_hdr_t* vcf_hdr = bcf_hdr_read(vcf_fp);
    bcf_hdr_set_samples(vcf_hdr, NULL, 0); // no genotypes needed here
    bcf1_t* vcf_rec = bcf_init();

    contig2map_map contig2vars;
    for (int32_t i = 0; i < vcf_hdr->n[BCF_DT_CTG]; ++i) {
        const char* seqk = vcf_hdr->id[BCF_DT_CTG][i].key;
        contig2vars.insert_or_assign(seqk, pos2var_map());
    }

    int32_t pid = -1;
    int32_t ppos = -1;
    pos2var_map* vmap = nullptr;
    VarList vs;
    while (!bcf_read(vcf_fp, vcf_hdr, vcf_rec)) {
        bcf_unpack(vcf_rec, BCF_UN_STR);
        if (vcf_rec->pos != ppos) {
            if (vmap && vs.size()) {
                vmap->insert_or_assign(ppos, std::move(vs)); // vs is undefined after this
            }
            if (vcf_rec->rid != pid) { // change hash tables here
                vmap = &(contig2vars[bcf_hdr_id2name(vcf_hdr, vcf_rec->rid)]);
            }
            vs.clear();
        }
        VarList more_vs = bcf_to_vars(vcf_rec);
        vs.insert(vs.end(), more_vs.begin(), more_vs.end());
        ppos = vcf_rec->pos;
        pid = vcf_rec->rid;
    }

    if (vcf_rec->pos == ppos && vmap && vs.size()) {
        vmap->insert_or_assign(ppos, std::move(vs));
    }
    vs.clear();

    pid = -1;
    bam1_core_t* c = nullptr;
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

            VarList aln_vars(bam_to_vars(aln));
            if (args.verbose) {
                fprintf(stderr, "a %s ", bam_get_qname(aln));
                print_varlist(aln_vars, stderr);
            }
            // can we vectorize or at least parallelize this? would it be worth it?
            // as of now, the performance bottleneck is the VCF reading
            for (int32_t i = c->pos; i <= bam_endpos(aln); ++i) {
                auto found = vmap->find(i);
                if (found != vmap->end()) {
                    Var* v_cached = &found->second[0];
                    for (auto&& vv: found->second) {
                        if (vv.rec_start) // keep pointer to this record
                            v_cached = &vv;
                        if (aln_vars.size()) {
                            for (const auto& av: aln_vars) {
                                if (var_match(vv, av)) v_cached->ac += 1;
                                else v_cached->rc += 1;
                            }
                        } else v_cached->rc +=1;
                    }
                    if (args.verbose) {
                        fprintf(stderr, "v %s ", bam_get_qname(aln));
                        print_varlist(found->second, stderr);
                    }
                }
            }
        }
    }

    bam_destroy1(aln);
    bam_hdr_destroy(sam_hdr);
    sam_close(sam_fp);


    // output in VCF format!!
    // std::string out_fname = args.sample_name + ".vcf";
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
    bcf_hdr_append(out_vcf_hdr, "##INFO=<ID=ALTCNT,Number=1,Type=Integer,Description=\"Count of reads covering alt allele\">");
    bcf_hdr_append(out_vcf_hdr, "##INFO=<ID=REFCNT,Number=1,Type=Integer,Description=\"Count of reads covering ref allele\">");
    bcf_hdr_append(out_vcf_hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_write(out_vcf_fp, out_vcf_hdr);

    bcf1_t* out_vcf_rec = bcf_init();

    // TODO: sort the output!
    for (const auto& vs: contig2vars) { // vs: {seq, map}
        for (const auto& v: vs.second)  { // v: {pos, Varlist}
            auto it = v.second.begin();
            while (it != v.second.end()) {
                std::string allele_str(it->ref + "," + it->alt);
                auto nit = std::next(it);
                while (nit != v.second.end() && !nit->rec_start) {
                    allele_str += "," + nit->alt;
                    ++nit;
                }
                out_vcf_rec->rid = bcf_hdr_name2id(out_vcf_hdr, vs.first.data());
                out_vcf_rec->pos = it->pos;
                out_vcf_rec->rlen = it->ref.size();
                bcf_update_alleles_str(out_vcf_hdr, out_vcf_rec, allele_str.data());
                bcf_update_id(out_vcf_hdr, out_vcf_rec, it->id.data());
                bcf_update_info_int32(out_vcf_hdr, out_vcf_rec, "ALTCNT", &(it->ac), 1);
                bcf_update_info_int32(out_vcf_hdr, out_vcf_rec, "REFCNT", &(it->rc), 1);
                int32_t gts[1];
                if (args.gt) {
                    if (std::abs(it->rc - it->ac) >= args.thres) {
                        if ((it->ac && args.thres < 0) || it->rc < it->ac) {
                            gts[0] = bcf_gt_phased(1);
                        } else { // we default to 0 in the case of a tie and thres==0
                            gts[0] = bcf_gt_phased(0);
                        }
                    } else if (args.keep) {
                        gts[0] = bcf_gt_missing;
                    }
                } else {
                    gts[0] = bcf_gt_missing;
                }
                bcf_update_genotypes(out_vcf_hdr, out_vcf_rec, gts, 1);
                bcf_write(out_vcf_fp, out_vcf_hdr, out_vcf_rec);
                bcf_clear(out_vcf_rec);

                it = nit;
            }
        }
    }
    bcf_hdr_destroy(out_vcf_hdr);
    bcf_close(out_vcf_fp);

    bcf_destroy(vcf_rec);
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
    vcnt::VcntArgs args;
    static struct option long_options[] {
        {"sample-name", required_argument, 0, 's'},
        {"threshold", required_argument, 0, 'c'},
        {"genotype", no_argument, &args.gt, 1},
        {"keep", no_argument, &args.keep, 1},
        {"verbose", no_argument, &args.verbose, 1},
        {"help", no_argument, 0, 'h'},
        {0,0,0,0}
    };

    int ch;
    int argpos = 0;
    while ( (ch = getopt_long(argc, argv, "-:s:c:gvkh", long_options, NULL)) != -1 ) {
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
            case 'g':
                args.gt = 1;
                break;
            case 'k':
                args.keep = 1;
                break;
            case 'v':
                args.verbose = 1;
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

    vcnt::varcount(args);
}
