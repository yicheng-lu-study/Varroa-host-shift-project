#!/bin/bash

VCF_FILE="/Users/easonlu/Desktop/project_pi/snps_0095.recode.vcf"
WINDOW_SIZE=100000                       # 滑动窗口大小（100kb）
WINDOW_STEP=10000                        # 滑动窗口步长（10kb）
POPULATIONS=("lavender-AmK1-core" "purple-AmK1-USA&ZAF" "blue-AcK1-CHN" "mix-AmK1-JPN" "red-AmK1&AcK1-JPN" "purple-AmK1-USA_ZAF" "J1" "V1" "C")              

echo "开始进行多种群θπ分析..."
echo "VCF文件: $VCF_FILE"
echo "窗口大小: $WINDOW_SIZE bp"
echo "窗口步长: $WINDOW_STEP bp"
echo "分析种群: ${POPULATIONS[@]}"
echo "========================================="

for POP in "${POPULATIONS[@]}"
do
    echo ""
    echo "正在处理种群: $POP"
    echo "-----------------------------------------"
    
    POP_FILE="${POP}.pop" 
    
    # 检查种群文件是否存在
    if [ ! -f "$POP_FILE" ]; then
        echo "警告: 种群文件 $POP_FILE 不存在，跳过..."
        continue
    fi
    
    echo "1. 计算全基因组θπ（每位点）..."
    vcftools --vcf "$VCF_FILE" \
        --keep "$POP_FILE" \
        --site-pi \
        --out "${POP}_genomewide_pi"
    
    if [ $? -ne 0 ]; then
        echo "错误: vcftools执行失败，跳过该种群"
        continue
    fi
    
    echo "2. 计算滑动窗口θπ（窗口: ${WINDOW_SIZE}bp, 步长: ${WINDOW_STEP}bp）..."
    vcftools --vcf "$VCF_FILE" \
        --keep "$POP_FILE" \
        --window-pi "$WINDOW_SIZE" \
        --window-pi-step "$WINDOW_STEP" \
        --out "${POP}_window_pi_${WINDOW_SIZE}_${WINDOW_STEP}"
    

    echo "3. 计算全基因组θπ统计摘要..."
    
    if [ -f "${POP}_genomewide_pi.sites.pi" ]; then
        # 计算统计量
        awk -v pop="$POP" '
        NR>1 {  # 跳过标题行
            pi=$3
            if(pi != "NaN" && pi != "nan" && pi != "-nan" && pi != "") {
                pi_num=pi+0
                sum+=pi_num
                sum_sq+=pi_num*pi_num
                count++
                if(pi_num<min || count==1) min=pi_num
                if(pi_num>max || count==1) max=pi_num
                # 存储值用于计算中位数
                values[count]=pi_num
            }
        }
        END {
            if(count>0) {
                mean=sum/count
                if(count>1) {
                    variance=(sum_sq - sum*sum/count)/(count-1)
                    stdev=sqrt(variance)
                } else {
                    stdev=0
                }
                

                n=count

                for(i=1; i<=n; i++) {
                    for(j=i+1; j<=n; j++) {
                        if(values[i] > values[j]) {
                            temp=values[i]
                            values[i]=values[j]
                            values[j]=temp
                        }
                    }
                }
                
                if(n%2==1) {
                    median=values[int(n/2)+1]
                } else {
                    median=(values[n/2] + values[n/2+1])/2
                }
                

                print "========================================" > "'"${POP}"'_pi_summary.txt'"
                print "种群: " pop > "'"${POP}"'_pi_summary.txt'"
                print "========================================" > "'"${POP}"'_pi_summary.txt'"
                print "分析位点数量: " count > "'"${POP}"'_pi_summary.txt'"
                print "θπ平均值: " mean > "'"${POP}"'_pi_summary.txt'"
                print "θπ标准差: " stdev > "'"${POP}"'_pi_summary.txt'"
                print "θπ最小值: " min > "'"${POP}"'_pi_summary.txt'"
                print "θπ最大值: " max > "'"${POP}"'_pi_summary.txt'"
                print "θπ中位数: " median > "'"${POP}"'_pi_summary.txt'"
                print "========================================" > "'"${POP}"'_pi_summary.txt'"
                
                # 同时打印到屏幕
                print "全基因组θπ统计:"
                print "  位点数量: " count
                printf "  平均值: %.6f\n", mean
                printf "  标准差: %.6f\n", stdev
                printf "  范围: [%.6f, %.6f]\n", min, max
                printf "  中位数: %.6f\n", median
            } else {
                print "警告: 没有找到有效的位点数据" > "'"${POP}"'_pi_summary.txt'"
                print "警告: 没有找到有效的位点数据"
            }
        }' "${POP}_genomewide_pi.sites.pi"
    else
        echo "错误: 文件 ${POP}_genomewide_pi.sites.pi 不存在"
    fi
    
    echo "4. 计算滑动窗口θπ统计摘要..."
    WINDOW_FILE="${POP}_window_pi_${WINDOW_SIZE}_${WINDOW_STEP}.windowed.pi"
    
    if [ -f "$WINDOW_FILE" ]; then
        # 计算统计量
        awk -v pop="$POP" -v win_size="$WINDOW_SIZE" -v win_step="$WINDOW_STEP" '
        BEGIN {
            sum=0; sum_sq=0; count=0; min=999999; max=0
        }
        NR>1 {  # 跳过标题行
            pi=$6  # 第6列是PI值
            if(pi != "NaN" && pi != "nan" && pi != "-nan" && pi != "") {
                pi_num=pi+0  # 转换为数字
                sum+=pi_num
                sum_sq+=pi_num*pi_num
                count++
                if(pi_num<min || count==1) min=pi_num
                if(pi_num>max || count==1) max=pi_num
                # 存储值用于计算中位数
                values[count]=pi_num
            }
        }
        END {
            if(count>0) {
                mean=sum/count
                if(count>1) {
                    variance=(sum_sq - sum*sum/count)/(count-1)
                    stdev=sqrt(variance)
                } else {
                    stdev=0
                }
                
                # 计算中位数
                n=count
                # 排序
                for(i=1; i<=n; i++) {
                    for(j=i+1; j<=n; j++) {
                        if(values[i] > values[j]) {
                            temp=values[i]
                            values[i]=values[j]
                            values[j]=temp
                        }
                    }
                }
                
                if(n%2==1) {
                    median=values[int(n/2)+1]
                } else {
                    median=(values[n/2] + values[n/2+1])/2
                }
                
                # 输出到文件
                print "========================================" > "'"${POP}"'_window_summary.txt'"
                print "种群: " pop > "'"${POP}"'_window_summary.txt'"
                print "窗口大小: " win_size " bp" > "'"${POP}"'_window_summary.txt'"
                print "窗口步长: " win_step " bp" > "'"${POP}"'_window_summary.txt'"
                print "========================================" > "'"${POP}"'_window_summary.txt'"
                print "有效窗口数量: " count > "'"${POP}"'_window_summary.txt'"
                print "平均窗口θπ: " mean > "'"${POP}"'_window_summary.txt'"
                print "窗口θπ标准差: " stdev > "'"${POP}"'_window_summary.txt'"
                print "窗口θπ最小值: " min > "'"${POP}"'_window_summary.txt'"
                print "窗口θπ最大值: " max > "'"${POP}"'_window_summary.txt'"
                print "窗口θπ中位数: " median > "'"${POP}"'_window_summary.txt'"
                print "========================================" > "'"${POP}"'_window_summary.txt'"
                
                # 打印到屏幕
                print "滑动窗口θπ统计:"
                print "  有效窗口数量: " count
                printf "  平均窗口θπ: %.6f\n", mean
                printf "  范围: [%.6f, %.6f]\n", min, max
                printf "  中位数: %.6f\n", median
            } else {
                print "警告: 没有找到有效的窗口数据" > "'"${POP}"'_window_summary.txt'"
                print "警告: 没有找到有效的窗口数据"
            }
        }' "$WINDOW_FILE"
    else
        echo "警告: 窗口结果文件 $WINDOW_FILE 不存在"
    fi
    
    echo "种群 $POP 分析完成！"
    echo "-----------------------------------------"
done


echo ""
echo "创建分析汇总报告..."

> analysis_summary.txt


for POP in "${POPULATIONS[@]}"
do
    if [ -f "${POP}_pi_summary.txt" ]; then
        echo "【种群 $POP - site π】" >> analysis_summary.txt
        cat "${POP}_pi_summary.txt" >> analysis_summary.txt
        echo "" >> analysis_summary.txt
    fi
    if [ -f "${POP}_window_summary.txt" ]; then
        echo "【种群 $POP - window π】" >> analysis_summary.txt
        cat "${POP}_window_summary.txt" >> analysis_summary.txt
        echo "" >> analysis_summary.txt
    fi
done

echo "分析汇总报告已保存到: analysis_summary.txt"