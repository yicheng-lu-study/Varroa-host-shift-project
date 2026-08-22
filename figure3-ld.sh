#!/bin/bash

INPUT_VCF="/data2/weirk/YSN/LYC/0095/snps_0095.recode.vcf"
OUTPUT_BASE="/data2/weirk/YSN/LYC/0095/ld"
ORIGINAL_CHROMS=(
    "NW_019211454.1"
    "NW_019211455.1"
    "NW_019211456.1"
    "NW_019211457.1"
    "NW_019211458.1"
    "NW_019211459.1"
    "NW_019211460.1"
)
LD_WINDOW=10
LD_STEP=5
LD_R2=0.5
GEN_BEAGLE="no"
BEAGLE_CHROMS=("${ORIGINAL_CHROMS[@]}")
THREADS=20
OUTDIR="${OUTPUT_BASE}/LDpruneR2-02"
mkdir -p "$OUTDIR" || { echo "无法创建输出目录: $OUTDIR"; exit 1; }


PREFIX="snps_0095"
READY_VCF="${OUTDIR}/${PREFIX}_ready4plink.vcf"
PLINK_PREFIX="${OUTDIR}/${PREFIX}_plink"
PRUNE_PREFIX="${OUTDIR}/${PREFIX}_pruned"
PRUNE_IN="${PRUNE_PREFIX}.prune.in"
PRUNE2KEEP="${OUTDIR}/prune2keep.txt"
LD_VCF="${OUTDIR}/${PREFIX}_LDpruned.vcf"
LD_VCF_RECODE="${OUTDIR}/${PREFIX}_LDpruned.recode.vcf"  

if [ ! -f "$INPUT_VCF" ]; then
    echo "错误: 输入VCF文件不存在: $INPUT_VCF"
    exit 1
fi

for cmd in vcftools plink; do
    if ! command -v $cmd &>/dev/null; then
        echo "错误: 未找到命令 $cmd，请安装或加入PATH"
        exit 1
    fi
done

cp "$INPUT_VCF" "$READY_VCF"
# 使用循环构建sed替换命令
for i in "${!ORIGINAL_CHROMS[@]}"; do
    chrom_orig="${ORIGINAL_CHROMS[$i]}"
    # 转义点号，因为sed中点为元字符
    chrom_orig_escaped=$(echo "$chrom_orig" | sed 's/\./\\./g')
    num=$((i+1))
    sed -i "s/${chrom_orig_escaped}/${num}/g" "$READY_VCF"
done

chr_params=""
for ((i=1; i<=${#ORIGINAL_CHROMS[@]}; i++)); do
    chr_params+=" --chr $i"
done
vcftools --vcf "$READY_VCF" --out "$PLINK_PREFIX" $chr_params --plink
if [ $? -ne 0 ]; then
    echo "错误: vcftools转换PLINK失败"
    exit 1
fi

plink --file "$PLINK_PREFIX" --indep-pairwise $LD_WINDOW $LD_STEP $LD_R2 --out "$PRUNE_PREFIX"
if [ $? -ne 0 ] || [ ! -f "$PRUNE_IN" ]; then
    echo "错误: PLINK LD修剪失败或未生成prune.in文件"
    exit 1
fi

cp "$PRUNE_IN" "$PRUNE2KEEP"
for i in "${!ORIGINAL_CHROMS[@]}"; do
    chrom_orig="${ORIGINAL_CHROMS[$i]}"
    num=$((i+1))

    sed -i "s/^${num}:/${chrom_orig}\t/" "$PRUNE2KEEP"
done

vcftools --vcf "$INPUT_VCF" --positions "$PRUNE2KEEP" --recode --recode-INFO-all --out "${OUTDIR}/${PREFIX}_LDpruned"
if [ $? -ne 0 ] || [ ! -f "$LD_VCF_RECODE" ]; then
    echo "错误: vcftools提取失败"
    exit 1
fi

mv "$LD_VCF_RECODE" "$LD_VCF"
echo "LD修剪完成，结果VCF: $LD_VCF"

