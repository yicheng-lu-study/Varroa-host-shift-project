vcftools --vcf /data2/weirk/YSN/LYC/0095/ld/LDpruneR2-05/snps_0095_LDpruned.vcf  --plink --out plink_out
plink --vcf /data2/weirk/YSN/LYC/0095/ld/LDpruneR2-05/snps_0095_LDpruned.vcf  -recode12 --out plink_out
 

plink --file plink_out --allow-extra-chr --snps-only --make-bed --out QC

for K in `seq 1 20`
do 
admixture --cv -B -s 1004 QC.bed $K  -j10 |tee log${K}.out
done
 

grep -h CV log*.out |sort -nk4  -t ' ' > cross-validation_error.txt


best_k=$(grep -h CV log*.out |sort -nk4  -t ' ' |head -n 1 |perl -p -e 's/.*K=(.+?)\):.*/$1/g')


for K in $(seq 1 20); do
    if [[ $K -lt 2 || $K -gt $best_k ]]; then
        rm -f QC.${K}.Q QC.${K}.P log${K}.out
    fi
done


for K in $(seq 2 $best_k); do
    if [[ -f "QC.${K}.Q" ]]; then
        paste <(cut -f 1 -d ' ' QC.fam) <(cat QC.${K}.Q | tr ' ' '\t') > QC.${K}.Q.tab
    fi
done