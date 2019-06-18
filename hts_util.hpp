#ifndef VARCOUNT_HPP
#define VARCOUNT_HPP

#include <cstdio>
#include <vector>
#include <list>
#include <string>
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <utility>
#include "mdparse.hpp"
#include "flat_hash_map.hpp"

namespace hts_util {
    enum class VTYPE {V_UNK, V_SNP, V_INS, V_DEL};


    struct Var {
        template<typename Ts>
        Var(int32_t p, VTYPE t, std::string a, std::string r, std::string i, Ts&&...) : pos(p), type(t), alt(a), ref(r), id(i) {}
        Var(int32_t p, VTYPE t, std::string a, std::string r, std::string i) : pos(p), type(t), alt(a), ref(r), id(i) {}
        Var(int32_t p, VTYPE t, std::string a, std::string r) : Var(p,t,a,r,"") {}
        // Var(Ts&&... args, std::string i) Var(std::forward args), id(i) {}
        int32_t pos = 0;
        VTYPE type = VTYPE::V_UNK;
        std::string alt = "";
        std::string ref = "";
        // the following WON'T be used in comparators!
        std::string id = ""; 
        int32_t rc = 0;
        int32_t ac = 0;
        bool rec_start = 0;
    };



    // trims s1.size()-1 characters from prefix of s1, s2 
    static inline std::tuple<std::string,std::string> truncate_str_pair(const std::string& s1, const std::string& s2) {
        return std::forward_as_tuple(s1.substr(s1.size()-1), s2.substr(s1.size()-1));
    }

    inline bool var_match(const Var& lv, const Var& rv) {
        if (lv.type != rv.type) return false;
        switch (lv.type) {
            case VTYPE::V_INS:
                return lv.pos+static_cast<int32_t>(lv.ref.size())-1 == rv.pos && truncate_str_pair(lv.ref, lv.alt) == truncate_str_pair(rv.ref, rv.alt) ;
            case VTYPE::V_DEL:
                return lv.pos == rv.pos && truncate_str_pair(lv.alt, lv.ref) == truncate_str_pair(rv.alt, rv.ref);
            case VTYPE::V_SNP:
                return lv.pos == rv.pos && lv.alt == rv.alt;
            default:
                fprintf(stderr, "no support for non SNPs & INDELs yet\n");
                return false;
        }
    }


    template <typename T=Var>
    using pos2var_map = ska::flat_hash_map< int32_t, std::vector<T> >;

    template <typename T=Var>
    using contig2map_map = ska::flat_hash_map<std::string, pos2var_map<T>>;

    // might encounter bug if STRLEN(REF) > 1 && STRLEN(ALT) > 1
    template <typename T>
    std::vector<T> bcf_to_vars(bcf_hdr_t* hdr, bcf1_t* b) {
        std::vector<T> vs;
        return vs;
    }
    
    template <>
    std::vector<Var> bcf_to_vars(bcf_hdr_t* hdr, bcf1_t* b) {
        (void) hdr;
        std::vector<Var> vs;
        char* ref = b->d.allele[0];
        for (uint32_t i = 1; i < b->n_allele; ++i) {
            char* alt = b->d.allele[i];
            if (alt[0] == '.') continue;
            if (strlen(alt)  < strlen(ref)) { // DEL
                vs.push_back(Var(b->pos, VTYPE::V_DEL, alt, ref, b->d.id)); // don't need alt here
            } else if (strlen(alt) > strlen(ref)) { // INS
                vs.push_back(Var(b->pos, VTYPE::V_INS, alt, ref, b->d.id)); // don't need ref here
            } else { // SNP
                vs.push_back(Var(b->pos, VTYPE::V_SNP, alt, ref, b->d.id)); // don't need ref here
            }
        }
        vs[0].rec_start = 1;
        return vs;
    }


    std::vector<Var> bam_to_vars(bam1_t* aln) {
        std::vector<Var> vs;
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


    template <typename T=Var>
    contig2map_map<T> bcf_to_map(vcfFile* vcf_fp, bcf_hdr_t* vcf_hdr) {
        bcf1_t* vcf_rec = bcf_init();
        contig2map_map<T> contig2vars;
        for (int32_t i = 0; i < vcf_hdr->n[BCF_DT_CTG]; ++i) {
            const char* seqk = vcf_hdr->id[BCF_DT_CTG][i].key;
            contig2vars.insert_or_assign(seqk, pos2var_map<T>());
        }

        int32_t pid = -1;
        int32_t ppos = -1;
        pos2var_map<T>* vmap = nullptr;
        std::vector<T> vs;
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
            std::vector<T> more_vs = bcf_to_vars<T>(vcf_hdr, vcf_rec);
            vs.insert(vs.end(), more_vs.begin(), more_vs.end());
            ppos = vcf_rec->pos;
            pid = vcf_rec->rid;
        }

        if (vcf_rec->pos == ppos && vmap && vs.size()) {
            vmap->insert_or_assign(ppos, std::move(vs));
        }
        vs.clear();
        bcf_destroy(vcf_rec);
        return contig2vars;
    }


    void print_varlist(std::vector<Var> vs, FILE* out) {
        for (const auto& v: vs) {
            fprintf(out, "(%d %d %s %s", v.pos, static_cast<int>(v.type), v.alt.data(), v.ref.data());
            if (v.id.size()) fprintf(out, " %s", v.id.data());
            if (v.rc + v.ac) fprintf(out, " %d %d", v.rc, v.ac);
            fprintf(out, ") ");
        } fprintf(out, "\n");
    }
};

#endif // VARCOUNT_HPP