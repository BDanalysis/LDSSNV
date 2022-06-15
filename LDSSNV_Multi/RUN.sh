###
 # @Author: cwx
 # @Date: 2022-05-16 17:32:28
 # @LastEditTime: 2022-05-18 10:26:50
 # @FilePath: /paper_snv/paper/LDSSNV/LDSSNV-Multi.sh
 # @Description: 在没有配对样本的情况下，使用多个肿瘤样本检测体细胞SNV
### 
# 1. build index and comparison to reference
bwa index -a bwtsw $reference
bwa mem -t 4 $reference $tumor1 $tumor2 > tumor.sam
# 2. sam to bam
samtools view -bS tumor.sam > tumor.bam
# 3. sort bam
samtools sort tumor.bam -o tumor.sort.bam
# 4. generate pileup file
samtools mpileup -f $reference tumor.sort.bam --output-QNAME > tumor.pileup

# 5. Detect CNVs using FREEC or other software
# 6. Get the CNV result file, please see README for csv file format

# 7. detect all SNVs
python detect_snv.py

# 8. get feature matrix
python get_feature_matrix.py

# 9. detect somatic SNVs
python detect_somatic_snv.py