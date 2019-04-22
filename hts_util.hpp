#ifndef VARCOUNT_HPP
#define VARCOUNT_HPP

#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <cstdio> 
#include <vector>
#include <string>
#include "flat_hash_map.hpp"

namespace hts_util {
    enum class VTYPE {V_UNK, V_SNP, V_INS, V_DEL};

    struct Var {
        Var(int32_t p, VTYPE t, std::string a, std::string r) : pos(p), type(t), alt(a), ref(r) {}
        Var(int32_t p, VTYPE t, std::string a, std::string r, std::string i) : pos(p), type(t), alt(a), ref(r), id(i) {}
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

    using VarList = std::vector<Var>;
    using pos2var_map = ska::flat_hash_map< int32_t, VarList >;
    using contig2map_map = ska::flat_hash_map<std::string, pos2var_map>;

    contig2map_map bcf_to_map(vcfFile* vcf_fp, bcf_hdr_t* vcf_hdr);

    VarList bam_to_vars(bam1_t* a);
    VarList bcf_to_vars(bcf1_t* b);

    void print_varlist(VarList vs, FILE* out);
};

#endif // VARCOUNT_HPP
