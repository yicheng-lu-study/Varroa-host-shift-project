#!/bin/bash

plink_bin="/root/anaconda3/envs/lupca/bin/plink"
run_id=$(date +%Y%m%d_%H%M%S)_$$

# 使用绝对路径，避免混淆不同目录下的 VCF 文件
vcf_files=(
    "/mnt/e/project_varroa_origin/pca_0624/LYC/0095/ld/LDpruneR2-05/snps_0095_LDpruned.vcf"
    "/mnt/e/project_varroa_origin/pca_0624/LYC/0095/ld/LDpruneR2-05/snps_0095_LDpruned111.vcf"
    "/mnt/e/project_varroa_origin/pca_0624/LYC/0095/ld/LDpruneR2-05/snps_0095_LDpruned39.vcf"
    "/mnt/e/project_varroa_origin/pca_0624/LYC/0095/ld/LDpruneR2-05/snps_0095_LDpruned80.vcf"
    "/mnt/e/project_varroa_origin/pca_0624/LYC/0095/ld/LDpruneR2-05/snps_0095_LDpruned86.vcf"
)

if [ ! -x "$plink_bin" ]; then
    echo "错误：未找到可执行的 PLINK：$plink_bin"
    echo "请先安装 lupca 环境，或检查 PLINK 路径。"
    exit 1
fi

# 循环处理每个文件
for vcf in "${vcf_files[@]}"; do
    echo "========================================"
    echo "正在处理: $vcf"

    # 检查文件是否存在
    if [ ! -f "$vcf" ]; then
        echo "错误：文件不存在，跳过！"
        continue
    fi

    # 获取VCF文件所在目录
    input_dir=$(dirname "$vcf")
    # 获取VCF文件名（不含路径和扩展名）
    vcf_basename=$(basename "$vcf" .vcf)
    # 输出目录：每次运行使用独立批次目录，避免覆盖旧结果
    output_dir="${input_dir}/pca_fullpc_${run_id}/${vcf_basename}"

    # 创建输出目录（如果不存在）
    mkdir -p "$output_dir"

    echo "输出目录: $output_dir"

    sample_count=$(awk -F'\t' '/^#CHROM/{print NF - 9; exit}' "$vcf")
    if [ -z "$sample_count" ] || [ "$sample_count" -lt 2 ]; then
        echo "错误：无法从 VCF 读取有效样本数，跳过本文件"
        continue
    fi
    pca_count=$((sample_count - 1))
    echo "样本数: $sample_count"
    echo "输出全部可用 PC 数: $pca_count"

    # 步骤1：将VCF转换为PLINK二进制格式
    "$plink_bin" --vcf "$vcf" \
                 --make-bed \
                 --allow-extra-chr \
                 --out "${output_dir}/plink_out"

    # 检查上一步是否成功（可选）
    if [ $? -ne 0 ]; then
        echo "错误：PLINK make-bed 失败，跳过本文件"
        continue
    fi

    # 步骤2：基于二进制文件计算全部可用主成分
    "$plink_bin" --bfile "${output_dir}/plink_out" \
                 --pca "$pca_count" \
                 --allow-extra-chr \
                 --out "${output_dir}/pca"

    if [ $? -eq 0 ]; then
        echo "完成: $vcf"
        awk 'BEGIN {
                 print "PC\teigenvalue\tnonnegative_eigenvalue\texplained_ratio\texplained_percent"
             }
             {
                 eigen[NR] = $1
                 nonnegative_eigen[NR] = ($1 < 0 ? 0 : $1)
                 total += nonnegative_eigen[NR]
             }
             END {
                 max = (NR < 10 ? NR : 10)
                 for (i = 1; i <= max; i++) {
                     ratio = (total > 0 ? nonnegative_eigen[i] / total : 0)
                     printf "PC%d\t%.10g\t%.10g\t%.10f\t%.6f\n", i, eigen[i], nonnegative_eigen[i], ratio, ratio * 100
                 }
             }' "${output_dir}/pca.eigenval" > "${output_dir}/pca.top10_explained_ratio.tsv"
        echo "前10个PC解释比例已写入: ${output_dir}/pca.top10_explained_ratio.tsv"
        awk 'BEGIN {
                 print "FID\tIID\tPC1\tPC2\tPC3\tPC4\tPC5"
             }
             {
                 print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7
             }' "${output_dir}/pca.eigenvec" > "${output_dir}/pca.pc1_5.tsv"
        echo "PC1-PC5 已写入: ${output_dir}/pca.pc1_5.tsv"
    else
        echo "错误：PLINK PCA 失败"
    fi
done

echo "所有文件处理完毕！"
