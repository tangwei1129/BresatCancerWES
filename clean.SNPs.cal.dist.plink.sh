### use plink to prune SNPs within 10000kb 10 SNP r2 0.1 ##
plink --bfile ../myvcf.snps.plink --indep-pairwise 10000 10 0.1

### extract the pruned SNPs and make a new file
plink --bfile ../myvcf.snps.plink --extract plink.prune.in --make-bed --out myvcf.snps.plink.ldprune

### caculate distance ####
plink --bfile myvcf.snps.plink.ldprune --distance square ibs 1-ibs allele-ct
