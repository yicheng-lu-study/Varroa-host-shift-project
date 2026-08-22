#!/bin/bash
# ==================================================
# 功能: 对VCF文件进行LD修剪，并可选输出BEAGLE格式,20260204使用版本
# 依赖: vcftools, plink
# ==================================================

# ----- 用户自定义参数 (请根据实际修改) -----
# 输入VCF文件路径
INPUT_VCF="/data2/weirk/YSN/LYC/0095/snps_0095.recode.vcf"
# 输出目录（将在此目录下创建LDprune子文件夹）
OUTPUT_BASE="/data2/weirk/YSN/LYC/0095/ld"
# 染色体原始名称列表（按顺序对应数字1..N）
ORIGINAL_CHROMS=(
    "NW_019211454.1"
    "NW_019211455.1"
    "NW_019211456.1"
    "NW_019211457.1"
    "NW_019211458.1"
    "NW_019211459.1"
    "NW_019211460.1"
)
# LD修剪参数 (plink --indep-pairwise)
LD_WINDOW=10
LD_STEP=5
LD_R2=0.2
#LD_R2=0.5
# 是否生成BEAGLE GL文件 (yes/no)
GEN_BEAGLE="no"
# BEAGLE染色体列表（若生成BEAGLE，需指定哪些染色体，通常与ORIGINAL_CHROMS一致）
BEAGLE_CHROMS=("${ORIGINAL_CHROMS[@]}")
# 线程数（用于并行处理BEAGLE转换）
THREADS=20
# ---------------------------------------------

# 设置输出子目录
OUTDIR="${OUTPUT_BASE}/LDpruneR2-02"
mkdir -p "$OUTDIR" || { echo "无法创建输出目录: $OUTDIR"; exit 1; }

# 定义临时文件前缀
PREFIX="snps_0095"
READY_VCF="${OUTDIR}/${PREFIX}_ready4plink.vcf"
PLINK_PREFIX="${OUTDIR}/${PREFIX}_plink"
PRUNE_PREFIX="${OUTDIR}/${PREFIX}_pruned"
PRUNE_IN="${PRUNE_PREFIX}.prune.in"
PRUNE2KEEP="${OUTDIR}/prune2keep.txt"
LD_VCF="${OUTDIR}/${PREFIX}_LDpruned.vcf"
LD_VCF_RECODE="${OUTDIR}/${PREFIX}_LDpruned.recode.vcf"  # vcftools输出会加上.recode

# 检查输入文件
if [ ! -f "$INPUT_VCF" ]; then
    echo "错误: 输入VCF文件不存在: $INPUT_VCF"
    exit 1
fi

# 检查依赖命令
for cmd in vcftools plink; do
    if ! command -v $cmd &>/dev/null; then
        echo "错误: 未找到命令 $cmd，请安装或加入PATH"
        exit 1
    fi
done

# ----- 步骤1: 复制VCF并替换染色体名称（原始名称 -> 数字）-----
echo "步骤1: 准备VCF，替换染色体名称..."
cp "$INPUT_VCF" "$READY_VCF"
# 使用循环构建sed替换命令
for i in "${!ORIGINAL_CHROMS[@]}"; do
    chrom_orig="${ORIGINAL_CHROMS[$i]}"
    # 转义点号，因为sed中点为元字符
    chrom_orig_escaped=$(echo "$chrom_orig" | sed 's/\./\\./g')
    num=$((i+1))
    sed -i "s/${chrom_orig_escaped}/${num}/g" "$READY_VCF"
done

# ----- 步骤2: 转换为PLINK格式（仅保留目标染色体）-----
echo "步骤2: 转换为PLINK格式..."
# 构建染色体参数：--chr 1 --chr 2 ... --chr N
chr_params=""
for ((i=1; i<=${#ORIGINAL_CHROMS[@]}; i++)); do
    chr_params+=" --chr $i"
done
vcftools --vcf "$READY_VCF" --out "$PLINK_PREFIX" $chr_params --plink
if [ $? -ne 0 ]; then
    echo "错误: vcftools转换PLINK失败"
    exit 1
fi

# ----- 步骤3: 使用PLINK进行LD修剪 -----
echo "步骤3: 执行LD修剪 (窗口=$LD_WINDOW, 步长=$LD_STEP, r²=$LD_R2)..."
plink --file "$PLINK_PREFIX" --indep-pairwise $LD_WINDOW $LD_STEP $LD_R2 --out "$PRUNE_PREFIX"
if [ $? -ne 0 ] || [ ! -f "$PRUNE_IN" ]; then
    echo "错误: PLINK LD修剪失败或未生成prune.in文件"
    exit 1
fi

# ----- 步骤4: 将保留SNP列表中的染色体数字还原为原始名称，并转换为两列格式-----
echo "步骤4: 还原染色体名称并准备位置列表..."
cp "$PRUNE_IN" "$PRUNE2KEEP"
for i in "${!ORIGINAL_CHROMS[@]}"; do
    chrom_orig="${ORIGINAL_CHROMS[$i]}"
    num=$((i+1))
    # 将 "数字:位置" 替换为 "原始名称\t位置"
    # 使用$'\t'插入制表符，确保跨平台兼容
    sed -i "s/^${num}:/${chrom_orig}\t/" "$PRUNE2KEEP"
done

# ----- 步骤5: 从原始VCF提取保留的SNP，生成LD修剪后的VCF -----
echo "步骤5: 提取保留SNP生成新VCF..."
vcftools --vcf "$INPUT_VCF" --positions "$PRUNE2KEEP" --recode --recode-INFO-all --out "${OUTDIR}/${PREFIX}_LDpruned"
if [ $? -ne 0 ] || [ ! -f "$LD_VCF_RECODE" ]; then
    echo "错误: vcftools提取失败"
    exit 1
fi
# 重命名或保持原样
mv "$LD_VCF_RECODE" "$LD_VCF"
echo "LD修剪完成，结果VCF: $LD_VCF"

