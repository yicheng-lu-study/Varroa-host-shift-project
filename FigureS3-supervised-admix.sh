
vcftools --vcf /data2/weirk/YSN/LYC/0095/ld/LDpruneR2-05/snps_0095_LDpruned.vcf  --plink --out plink_out
plink --vcf /data2/weirk/YSN/LYC/0095/ld/LDpruneR2-05/snps_0095_LDpruned.vcf  -recode12 --out plink_out
 

plink --file plink_out --allow-extra-chr --snps-only --make-bed --out QC

admixture --supervised QC.bed 2 -j20 | tee log_supervised_K2.out

paste <(cut -d ' ' -f1 QC.fam) <(tr ' ' '\t' < QC.2.Q) > QC.2.Q.tab