#ifndef VARCOUNT_HPP
#define VARCOUNT_HPP

#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <cstdio> 
#include <vector>
#include <string>

enum VTYPE {V_UNK, V_SNP, V_INS, V_DEL};

struct Var {
    Var(int32_t p, VTYPE t, std::string a, std::string r) : pos(p), type(t), ref(r), alt(a) {}
    Var(int32_t p, VTYPE t, std::string a, std::string r, std::string i) : pos(p), type(t), ref(r), alt(a), id(i) {}
    int32_t pos = 0;
    VTYPE type = V_UNK;
    std::string alt = "";
    std::string ref = "";
    // the following WON'T be used in comparators!
    std::string id = ""; 
    uint32_t rc = 0;
    uint32_t ac = 0;
};

// NOTE: we compare the ref SIZE rather than the string to account for the way VCF handles dels
// (VCF records del at pos before deletion, BAM/MD records del at pos at deletion)
// TODO: find a better way to handle this
inline bool operator==(const Var& l, const Var& r) {
    size_t s1 = l.ref.size(), s2 = r.ref.size();
    return std::tie(l.pos, l.type, l.alt, s1) == std::tie(r.pos, r.type, r.alt, s2);
}

// NOTE: we compare the ref SIZE rather than the string to account for the way VCF handles dels
// TODO: find a better way to handle this
inline bool operator<(const Var& l, const Var& r) {
    size_t s1 = l.ref.size(), s2 = r.ref.size();
    return std::tie(l.pos, l.type, l.alt, s1) < std::tie(r.pos, r.type, r.alt, s2);
}


using VarList = std::vector<Var>;
using pos2var_map = ska::flat_hash_map< int32_t, VarList >;
using contig2map_map = ska::flat_hash_map<std::string, pos2var_map>;

VarList bam_to_vars(bam1_t* a);
VarList bcf_to_vars(bcf1_t* b);
VarList intersect(const VarList& lhs, const VarList& rhs);

#endif // VARCOUNT_HPP
