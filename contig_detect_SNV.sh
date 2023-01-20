
###description: detect SNV with a single tumor sample

# 1. transfer contig.fa to fastq file
wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/fasta-to-fastq/fasta_to_fastq.pl
perl fasta_to_fastq.pl contig.fasta > tumor.fastq

# 2. build index and comparison to reference
bwa index -a bwtsw $reference
bwa mem -t 4 reference.fa tumor.fastq > tumor.sam
# 3. sam to bam
samtools view -bS tumor.sam > tumor.bam
# 4. sort bam
samtools sort tumor.bam -o tumor.sort.bam
# 5. generate pileup file
samtools mpileup -f $reference tumor.sort.bam --output-QNAME > tumor.pileup

# 6. Detect CNVs using FREEC or other software
# 7. Get the CNV result file, please see ReaDMe for csv file format

# 8. detect all SNVs
python detect_snv.py
