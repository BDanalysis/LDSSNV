 # @Description:detect SNV with a single tumor sample
### 

# 1. build index and comparison to reference
bwa index -a bwtsw reference.fa
# pair-end alignment
bwa mem -t 4 reference.fa tumor_1.fastq tumor_2.fastq > tumor.sam
# 2. sam to bam
samtools view -bS tumor.sam > tumor.bam
# 3. sort bam
samtools sort tumor.bam -o tumor.sort.bam
# 4. generate pileup file
samtools mpileup -f $reference tumor.sort.bam --output-QNAME > tumor.pileup

# 5. Detect CNVs using FREEC or other software
# 6. Get the CNV result file, please see ReaDMe for csv file format

# 7. detect all SNVs
python detect_snv.py
