#include <cstdio>
#include <getopt.h>
#include "hts_util.hpp"

struct CmpVcfArgs {
    std::string fname1;
    std::string fname2;
    std::string sample1;
    std::string sample2;
};

struct VarGT {
    VarGT(int32_t p,
          std::string a, std::string r, 
          int32_t g1, int32_t g2) : 
        pos(p),
        alt(a), ref(r), 
        gt1(g1), gt2(g2) {}
    int32_t pos = 0;
    std::string alt = "";
    std::string ref = "";
    int32_t gt1;
    int32_t gt2;
    // the following WON'T be used in comparators!
    bool rec_start = 0;

    static std::vector<VarGT> from_bcf(bcf_hdr_t* hdr, bcf1_t* b) {
        std::vector<VarGT> vs;
        char* ref = b->d.allele[0];
        int ngt, gt1, gt2;
        /* TODO: extract genotype from the hdr */
        int32_t* gts = hts_util::get_genotype(hdr, b, &ngt);
        gt1 = bcf_gt_allele(gts[0]);
        gt2 = bcf_gt_allele(gts[1]);
        if (ngt < 2) { fprintf(stderr, "error getting genotype!\n"); exit(1);}
        for (uint32_t i = 1; i < b->n_allele; ++i) {
            char* alt = b->d.allele[i];
            if (alt[0] == '.') continue;
            vs.push_back(VarGT(b->pos, alt, ref, gt1, gt2)); 
        }
        vs[0].rec_start = 1;
        return vs;
    }
};

static inline bool var_match(const VarGT& lv, const VarGT& rv) {
    return lv.pos == rv.pos && lv.alt == rv.alt && lv.ref == rv.ref && lv.gt1 == rv.gt1 && lv.gt2 == rv.gt2;
}

void compare_vcfs(CmpVcfArgs args) {
    vcfFile* fp1 = bcf_open(args.fname1.data(), "r");
    bcf_hdr_t* hdr1 = bcf_hdr_read(fp1);
    bcf_hdr_set_samples(hdr1, args.sample1.data(), 0);

    hts_util::contig2map_map<VarGT> contig2vars(hts_util::bcf_to_map<VarGT>(fp1, hdr1));

    vcfFile* fp2 = bcf_open(args.fname2.data(), "r");
    bcf_hdr_t* hdr2 = bcf_hdr_read(fp2);
    bcf_hdr_set_samples(hdr2, args.sample2.data(), 0); // no genotypes needed here
    bcf1_t* rec = bcf_init();

    vcfFile* out_fp = bcf_open("-", "w");
    bcf_hdr_t* out_hdr = bcf_hdr_dup(hdr2);
    // print the header of fp2
    bcf_hdr_write(out_fp, out_hdr);

    int32_t pid = -1;
    hts_util::pos2var_map<VarGT>* vmap = nullptr;
    // We're going to assume one variant per line for both VCF files for now.
    while (!bcf_read(fp2, hdr2, rec)) {
        bcf_unpack(rec, BCF_UN_STR || BCF_UN_FMT);
        auto vars2 = VarGT::from_bcf(hdr2, rec);
        if (rec->rid != pid) {
            vmap = &(contig2vars[bcf_hdr_id2name(hdr2, rec->rid)]);
        }
        std::string id(rec->d.id);
        auto vars1_iter = vmap->find(rec->pos);
        bool match = 0;
        if (vars1_iter != vmap->end()) {
            for (const auto& v2: vars2) {
                for (const auto& v1: vars1_iter->second) {
                    if (var_match(v1, v2)) {
                        match = 1;
                        break;
                    }
                }
                if (match) { // write record here
                    bcf_write(out_fp, out_hdr, rec);
                    break;
                } 
            } match = 0;
        }  // not found
        pid = rec->rid;
    }
}

void print_help() {
    fprintf(stderr, "usage: ./compare_vcfs -x sample1 -y sample2 <1.vcf> <2.vcf>...\n\
for all genotyped variants for sample2 in 2.vcf, checks if genotype matches\n\
that for sample 1 in 1.vcf. Ungenotyped variants in both files are ignored.\n\
\n\
output:\n\
matches(int) diffs(int) misses(int)\n");
}

int main(int argc, char** argv) {
    CmpVcfArgs args;
    static struct option long_options[] {
        {"sample1", required_argument, 0, 'x'},
        {"sample2", required_argument, 0, 'y'},
        {"help", no_argument, 0, 'h'},
        {0,0,0,0}
    };

    int ch;
    int argpos = 0;
    while ( (ch = getopt_long(argc, argv, "-:x:y:h", long_options, NULL)) != -1 ) {
        switch(ch) {
            case 0:
                break;
            case 1:
                if (argpos == 0) args.fname1 = std::string(optarg);
                else if (argpos == 1) args.fname2 = std::string(optarg);
                else fprintf(stderr, "ignoring argument %s\n", optarg);
                ++argpos;
                break;
            case 2:
                break;
            case 'x':
                args.sample1 = std::string(optarg); 
                break;
            case 'y':
                args.sample2 = std::string(optarg); 
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

    if (args.fname1 == "" || args.fname2 == "") {
        print_help();
        fprintf(stderr, "error: vcf files are mandatory arguments\n");
        exit(1);
    }

    if (args.sample1 == "" || args.sample2 == "") {
        print_help();
        fprintf(stderr, "error: sample names are mandatory arguments\n");
        exit(1);
    }

    compare_vcfs(args);
}

