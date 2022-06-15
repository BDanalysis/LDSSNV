###
 # @Author: cwx
 # @Date: 2022-05-16 17:32:28
 # @LastEditTime: 2022-05-18 10:24:54
 # @FilePath: /paper_snv/paper/LDSSNV/LDSSNV-Single.sh
 # @Description: 在没有配对样本的情况下，使用单肿瘤样本检测体细胞SNV
### 

reference="chr21.fa"
tumor1="tumor1.fastq"
tumor2="tumor2.fastq"

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

# 8. distinguish germline and somatic SNV
python detect_somatic_snv.py

# result

