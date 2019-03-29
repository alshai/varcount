#include <cstdio>
#include <vector>
#include <list>
#include <string>
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include "mdparse.hpp"
#include "flat_hash_map.hpp"
#include "varcount.hpp"


vcnt::VarList vcnt::bam_to_vars(bam1_t* aln) {
    vcnt::VarList vs;
    // look at md string and gather dels, snps
    char* md;
    std::list<std::pair<std::string, int32_t>> snps;
    std::list<std::pair<std::string, int32_t>> dels;
    if ((md = bam_aux2Z(bam_aux_get(aln, "MD")))) {
        std::vector<MDPos> mds = md_parse(md);
        for (auto m: mds) {
            if (m.st == MD_DEL) {
                std::string ref = "X"; // place holder for that extra base we have to account for in a deletion
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
            ins_seq += seq_nt16_str[bam_seqi(bam_get_seq(aln), qpos-1)]; // make sure we add the previous seq here
            for (uint32_t j = 0; j < bam_cigar_oplen(cs[i]); ++j) {
                ins_seq += seq_nt16_str[bam_seqi(bam_get_seq(aln), qpos+j)];
            }
            vs.push_back(vcnt::Var(rpos-1, vcnt::VTYPE::V_INS, ins_seq, ins_seq.substr(0,1)));
        } else if (bam_cigar_op(cs[i]) == BAM_CMATCH) {
            // check potential snps here
            for (auto it = snps.begin(); it != snps.end(); ) {
                int32_t s = it->second;
                if (s >= rpos && s < rpos + rlen) {
                    std::string snp = "";
                    snp += seq_nt16_str[bam_seqi(bam_get_seq(aln), qpos + (s - rpos))];
                    vs.push_back(vcnt::Var(s, vcnt::VTYPE::V_SNP, snp, it->first));
                    it = snps.erase(it); // we do this so we don't recheck snps that don't need rechecking, but... maybe deleting it would actually be more expensive?
                } else ++it;
            }
        } else if (bam_cigar_op(cs[i]) == BAM_CDEL) {
            for (auto it = dels.begin(); it != dels.end(); ) {
                int32_t d = it->second;
                if (d == rpos - 1) {
                    std::string alt = "";
                    alt += seq_nt16_str[bam_seqi(bam_get_seq(aln), qpos - 1)];
                    vs.push_back(vcnt::Var(d, vcnt::VTYPE::V_DEL, alt, it->first));
                    it = dels.erase(it);
                } else ++it;
            }
        }
        rpos += rlen;
        qpos += qlen;
    }
    std::sort(vs.begin(), vs.end()); // we probably don't need to sort
    return vs;
}

// might encounter bug if STRLEN(REF) > 1 && STRLEN(ALT) > 1
vcnt::VarList vcnt::bcf_to_vars(bcf1_t* b) {
    vcnt::VarList vs;
    char* ref = b->d.allele[0];
    for (uint32_t i = 1; i < b->n_allele; ++i) {
        char* alt = b->d.allele[i];
        if (alt[0] == '.') continue;
        if (strlen(alt)  < strlen(ref)) { // DEL
            vs.push_back(vcnt::Var(b->pos, vcnt::VTYPE::V_DEL, alt, ref, b->d.id)); // don't need alt here
        } else if (strlen(alt) > strlen(ref)) { // INS
            vs.push_back(vcnt::Var(b->pos, vcnt::VTYPE::V_INS, alt, ref, b->d.id)); // don't need ref here
        } else { // SNP
            vs.push_back(vcnt::Var(b->pos, vcnt::VTYPE::V_SNP, alt, ref, b->d.id)); // don't need ref here
        }
    }
    return vs;
}

void print_varlist(vcnt::VarList vs, FILE* out) {
    for (auto v: vs) {
        fprintf(out, "(%d %d %s %s", v.pos, v.type, v.alt.data(), v.ref.data());
        if (v.id.size()) fprintf(out, " %s", v.id.data());
        fprintf(out, ") ");
    } fprintf(out, "\n");
}

void vcnt::varcount(const vcnt::VcntArgs& args) {
    samFile* sam_fp = sam_open(args.sam_fname.data(), "r");
    bam_hdr_t* sam_hdr = sam_hdr_read(sam_fp);
    bam1_t* aln = bam_init1(); 
    
    vcfFile* vcf_fp = bcf_open(args.vcf_fname.data(), "r");
    bcf_hdr_t* vcf_hdr = bcf_hdr_read(vcf_fp);
    bcf_hdr_set_samples(vcf_hdr, NULL, 0); // no genotypes needed here
    bcf1_t* vcf_rec = bcf_init();

    vcnt::contig2map_map contig2vars;
    for (int32_t i = 0; i < vcf_hdr->n[BCF_DT_CTG]; ++i) {
        const char* seqk = vcf_hdr->id[BCF_DT_CTG][i].key;
        contig2vars.insert_or_assign(seqk, vcnt::pos2var_map());
    }

    int32_t pid = -1;
    int32_t ppos = -1;
    vcnt::pos2var_map* vmap = nullptr;
    vcnt::VarList vs;
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
        vcnt::VarList more_vs = vcnt::bcf_to_vars(vcf_rec);
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
                vmap = &(contig2vars[sam_hdr->target_name[c->tid]]);
            }
            pid = c->tid;

            vcnt::VarList aln_vars(vcnt::bam_to_vars(aln));
            // can we vectorize or at least parallelize this? would it be worth it?
            for (int32_t i = c->pos; i <= bam_endpos(aln); ++i) {
                auto found = vmap->find(i);
                if (found != vmap->end()) {
                    for (auto&& vv: found->second) {
                        for (const auto& av: aln_vars) {
                            if (vv == av) vv.ac += 1;
                            else vv.rc += 1;
                        }
                    }
                }
            }
        }
    }

    bam_destroy1(aln);
    bam_hdr_destroy(sam_hdr);
    sam_close(sam_fp);


    // output in VCF format!!
    std::string out_fname = args.sample_name + ".vcf";
    vcfFile* out_vcf_fp = bcf_open(out_fname.data(), "w");
    bcf_hdr_t* out_vcf_hdr = bcf_hdr_dup(vcf_hdr);

    bcf_hdr_add_sample(out_vcf_hdr, args.sample_name.data());
    bcf_hdr_remove(out_vcf_hdr, BCF_HL_INFO, NULL);
    bcf_hdr_remove(out_vcf_hdr, BCF_HL_GEN, NULL);
    bcf_hdr_remove(out_vcf_hdr, BCF_HL_STR, NULL);
    bcf_hdr_remove(out_vcf_hdr, BCF_HL_FMT, NULL);
    bcf_hdr_set_version(out_vcf_hdr, "4.3"); // TODO: make sure this is printed at the top
    // bcf_hdr_append(out_vcf_hdr, "##fileformat=VCFv4.3");
    bcf_hdr_append(out_vcf_hdr, "##INFO=<ID=ALTCNT,Number=1,Type=Integer,Description=\"Count of reads covering alt allele\">");
    bcf_hdr_append(out_vcf_hdr, "##INFO=<ID=REFCNT,Number=1,Type=Integer,Description=\"Count of reads covering ref allele\">");
    bcf_hdr_append(out_vcf_hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_write(out_vcf_fp, out_vcf_hdr);

    bcf1_t* out_vcf_rec = bcf_init();

    // TODO: sort the output!
    for (const auto& vs: contig2vars) { // vs: {seq, map}
        for (const auto& v: vs.second)  { // v: {pos, Varlist}
            for (const auto& x: v.second) { // x: {Var}
                // only output genotype that we know from coverage evidence!
                if ((x.rc + x.ac) && std::abs(x.rc - x.ac) >= args.thres) {
                    out_vcf_rec->rid = bcf_hdr_name2id(out_vcf_hdr, vs.first.data());
                    out_vcf_rec->pos = v.first;
                    out_vcf_rec->rlen = x.ref.size();
                    std::string allele_str(x.ref + "," + x.alt);
                    bcf_update_alleles_str(out_vcf_hdr, out_vcf_rec, allele_str.data());
                    bcf_update_id(out_vcf_hdr, out_vcf_rec, x.id.data());
                    bcf_update_info_int32(out_vcf_hdr, out_vcf_rec, "ALTCNT", &x.ac, 1);
                    bcf_update_info_int32(out_vcf_hdr, out_vcf_rec, "REFCNT", &x.rc, 1);
                    int32_t gts[2];
                    if (x.rc >= x.ac) {
                        gts[0] = bcf_gt_phased(0);
                        gts[1] = bcf_gt_phased(0);
                    } else {
                        gts[0] = bcf_gt_phased(1);
                        gts[1] = bcf_gt_phased(1);
                    }
                    bcf_update_genotypes(out_vcf_hdr, out_vcf_rec, gts, 2);
                    bcf_write(out_vcf_fp, out_vcf_hdr, out_vcf_rec);
                    bcf_clear(out_vcf_rec);
                }
            }
        }
    }
    bcf_hdr_destroy(out_vcf_hdr);
    bcf_close(out_vcf_fp);

    bcf_destroy(vcf_rec);
    bcf_hdr_destroy(vcf_hdr);
    bcf_close(vcf_fp);

}

int main(int argc, char** argv) {
    if (argc < 5) {
        fprintf(stderr, "usage: ./program <vcf> <sam> output_prefix threshold\n");
        exit(1);
    }
    vcnt::VcntArgs args;
    args.vcf_fname = std::string(argv[1]);
    args.sam_fname = std::string(argv[2]);
    args.sample_name = std::string(argv[3]);
    args.thres = std::atoi(argv[4]);
    vcnt::varcount(args);
}
