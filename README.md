## 1. Environment Installation（Linux）
1. BWA
2. Samtools
3. XGBoost
4. pandas
5. csv
6. numpy

### 1.1 using conda 
```
sudo apt install bwa
sudo apt install samtools

conda install pandas csv numpy XGBoost
```
### 1.2 using pip
```
sudo apt install bwa
sudo apt install samtools

pip3 install xgboost
pip3 install pandas
pip3 install csv
pip3 install numpy
```

## 2. Run

### 2.1 SNV Detection

detect_SNV.sh (with short reads)

contig_detect_SNV.sh (with contigs)

### 2.2 Distinguish somatic-germline SNV (single-mode)
With the SNVs detected from a single sample in 2.1, run: 

python LDSSNV-single/detect_somatic_snv.py

### 2.3 Distinguish somatic-germline SNV (Multiple-mode)
With the SNVs detected from each of multiple samples in 2.1,  get feature matrix and distinguish somatic SNV: 

python LDSSNV-Multi/get_feature_matrix.py
python LDSSNV-Multi/detect_somatic_snv.py

### 2.4 SNV detection and distinguishing combined.

LDSSNV-Single/run.sh
LDSSNV-Multi/run.sh

## 4. generate multiple simulation samples

### 4.1 using the LD simulator
using germline variant-site files produced by LD_simulator/gen_LD_SNP.py, perform the simulation tool SInC to insert somatic copy number variation and indels and then generate sequencing reads.

### 4.2 using variant-site information from real data

using somatic and germline variant-site files from real data, perform the simulation tool SInC to insert copy number variation and indels and then generate sequencing reads.


## 5. Others
### 5.1 datasets download

Click on the following link to save, or open the "Alibaba Cloud Disk" APP. You don't need to download the fast online view, and the original video will be played at double speed.
「LDSSNV_DATA」https://www.aliyundrive.com/s/4yRKtX23mQo

### 5.2 meaning of numbers in files
> 0 - germline SNV
> 1 - somatic SNV
> 2 - Homozygous
> 3 - Heterozygous 
> 4 - real SNV
> 5 - false SNV
