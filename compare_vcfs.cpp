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
    VarGT(int32_t p, std::string i, int g) : pos(p), id(i), gt(g) {}
    int32_t pos = 0;
    std::string id = ""; 
    int gt = 0;
};


int get_first_genotype(bcf_hdr_t* hdr, bcf1_t* b) {
    int ngt;
    int32_t *gt_arr = NULL, ngt_arr = 0;
    ngt = bcf_get_genotypes(hdr, b, &gt_arr, &ngt_arr);
    if (ngt <= 0) { return -1; }
    return bcf_gt_allele(gt_arr[0]);
}

namespace hts_util {
template <>
std::vector<VarGT> bcf_to_vars(bcf_hdr_t* hdr, bcf1_t* b) {
    std::vector<VarGT> vs;
    // get the first genotype
    int gt = get_first_genotype(hdr, b);
    vs.push_back(VarGT(b->pos, b->d.id, gt));
    return vs;
}
};


void compare_vcfs(CmpVcfArgs args) {
    vcfFile* vcf_fp1 = bcf_open(args.fname1.data(), "r");
    bcf_hdr_t* vcf_hdr1 = bcf_hdr_read(vcf_fp1);
    bcf_hdr_set_samples(vcf_hdr1, args.sample1.data(), 0);

    hts_util::contig2map_map<VarGT> contig2vars(hts_util::bcf_to_map<VarGT>(vcf_fp1, vcf_hdr1));

    vcfFile* vcf_fp2 = bcf_open(args.fname2.data(), "r");
    bcf_hdr_t* vcf_hdr2 = bcf_hdr_read(vcf_fp2);
    bcf_hdr_set_samples(vcf_hdr2, args.sample2.data(), 0); // no genotypes needed here
    bcf1_t* brec = bcf_init();

    uint32_t fps = 0;
    uint32_t fns = 0;
    uint32_t tps = 0;
    uint32_t tns = 0;
    int32_t pid = -1;
    hts_util::pos2var_map<VarGT>* vmap = nullptr;
    while (!bcf_read(vcf_fp2, vcf_hdr2, brec)) {
        bcf_unpack(brec, BCF_UN_STR || BCF_UN_FMT);
        if (brec->rid != pid) {
            vmap = &(contig2vars[bcf_hdr_id2name(vcf_hdr2, brec->rid)]);
        }
        std::string id(brec->d.id);
        int gt = get_first_genotype(vcf_hdr2, brec);
        auto vars = vmap->find(brec->pos);
        bool found = vars != vmap->end();
        if (found) {
            for (const auto& v: vars->second) {
                // count missing genotypes as fps (if gt=0) or fns (if gt=1).
                if (v.id == id) {
                    if (gt < 0) {
                        v.gt ? ++fns : ++fps;
                    } else if (!gt) {
                        gt == v.gt ? ++tns : ++fns;
                    } else if (gt > 0) {
                        gt == v.gt ? ++tps : ++fps;
                    }
                    found = true;
                    break;
                } else found = false;
            }
        } 
        if (!found && gt >= 0) { // genotyped but not in reference
            gt ? ++fps : ++fns;
        }

        pid = brec->rid;
    }
    fprintf(stderr, "tps\ttns\tfps\tfns\n");
    fprintf(stdout, "%d\t%d\t%d\t%d\n", tps, tns, fps, fns);
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

