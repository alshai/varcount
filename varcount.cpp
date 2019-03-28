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


VarList bam_to_vars(bam1_t* aln) {
    VarList vs;
    // look at md string and gather dels, snps
    char* md;
    std::list<int32_t> snps;
    if (md = bam_aux2Z(bam_aux_get(aln, "MD"))) {
        std::vector<MDPos> mds = md_parse(md);
        for (auto m: mds) {
            // we can just directly add deletions as Vars
            if (m.st == MD_DEL) {
                std::string ref = "X"; // place holder for that extra base we have to account for
                ref += m.str;
                vs.push_back(Var(aln->core.pos + m.p - 1, V_DEL, "", ref));
            } // snps need a little more processing to get the query base of the snp
            else if (m.st == MD_SNP) {
                snps.push_back(aln->core.pos + m.p);
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
            ins_seq += seq_nt16_str[bam_seqi(bam_get_seq(aln), qpos-1)];
            for (int j = 0; j < bam_cigar_oplen(cs[i]); ++j) {
                ins_seq += seq_nt16_str[bam_seqi(bam_get_seq(aln), qpos+j)];
            }
            vs.push_back(Var(rpos, V_INS, ins_seq, ""));
        } else if (bam_cigar_op(cs[i]) == BAM_CMATCH) {
            // check potential snps here
            for (auto it = snps.begin(); it != snps.end(); ) {
                int32_t s = *it;
                if (s >= rpos && s < rpos + rlen) {
                    std::string snp = "";
                    snp += seq_nt16_str[bam_seqi(bam_get_seq(aln), qpos + (s - rpos))];
                    vs.push_back(Var(s, V_SNP, snp, ""));
                    it = snps.erase(it); // we do this so we don't recheck snps that don't need rechecking, though deleting it might be more expensive?
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
VarList bcf_to_vars(bcf1_t* b) {
    VarList vs;
    char* ref = b->d.allele[0];
    for (uint32_t i = 1; i < b->n_allele; ++i) {
        char* alt = b->d.allele[i];
        if (alt[0] == '.') continue;
        if (strlen(alt)  < strlen(ref)) { // DEL
            vs.push_back(Var(b->pos, V_DEL, "", ref, b->d.id)); // don't need alt here
        } else if (strlen(alt) > strlen(ref)) { // INS
            vs.push_back(Var(b->pos, V_INS, alt, "", b->d.id)); // don't need ref here
        } else { // SNP
            vs.push_back(Var(b->pos, V_SNP, alt, "", b->d.id)); // don't need ref here
        }
    }
    return vs;
}

void print_varlist(VarList vs, FILE* out) {
    for (auto v: vs) {
        fprintf(out, "(%d %d %s %s", v.pos, v.type, v.alt.data(), v.ref.data());
        if (v.id.size()) fprintf(out, " %s", v.id.data());
        fprintf(out, ") ");
    } fprintf(out, "\n");
}

// invariant: a position on the chr corresponds to exactly one key/value pair in this map

void varcount(const char* vcf_fname, const char* sam_fname) {
    int ret;
    samFile* sam_fp = sam_open(sam_fname, "r");
    bam_hdr_t* sam_hdr = sam_hdr_read(sam_fp);
    bam1_t* aln = bam_init1(); 
    
    vcfFile* vcf_fp = bcf_open(vcf_fname, "r");
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
    bcf1_t* b = nullptr;
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
        // vcounts.insert_or_assign(std::string(vcf_rec->d.id), std::pair<uint32_t, uint32_t>(0,0));
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
                vmap = &(contig2vars[sam_hdr->target_name[c->tid]]);
            }
            pid = c->tid;

            VarList aln_vars(bam_to_vars(aln));
            for (uint32_t i = c->pos; i <= bam_endpos(aln); ++i) {
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

    for (const auto& vs: contig2vars) {
        for (const auto& v: vs.second)  {
            for (const auto& x: v.second) {
                fprintf(stdout, "%s %d %s %s %lu %lu\n", vs.first.data(), v.first, x.ref.data(), x.alt.data(), x.rc, x.ac);
            }
        }
    }
}

int main(int argc, char** argv) {
if (argc < 3) {
    fprintf(stderr, "usage: ./program <vcf> <sam>\n");
    exit(1);
    }
    varcount(argv[1], argv[2]);
}
