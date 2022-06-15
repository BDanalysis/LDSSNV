<!--
 * @Author: cwx
 * @Date: 2022-05-16 17:32:28
 * @LastEditTime: 2022-05-16 19:11:06
 * @FilePath: /paper_snv/paper/LDSSNV/README-CH.md
 * @Description: LDSSNV
-->
## LDSSNV
DSSNV: A Linkage Disequilibrium-Based Method for the Detection of Somatic Single-nucleotide Variants


## 1. 安装环境（Linux）
如有以下安装包请忽略
1. BWA
2. Samtools
3. XGBoost
4. pandas
5. csv
6. numpy

### 1.1 conda 安装(推荐)
```
sudo apt install bwa
sudo apt install samtools

conda install pandas csv numpy XGBoost
```
### 1.2 pip 安装
```
sudo apt install bwa
sudo apt install samtools

pip3 install xgboost
pip3 install pandas
pip3 install csv
pip3 install numpy
```
## 2. 下载
```
git clone https://github.com/BDanalysis/LDSSNV
cd LDSSNV
```

## 3. 文件结构


## 4. 运行
检测体细胞SNV, 分为单样本模式和多样本模式。
### 4.1 单样本检测somatic SNV

### 4.2 多样本检测somatic SNV

## 5. 其它
### 5.1 测试数据下载及说明

#### 5.1.1 下载
「LDSSNV_DATA」https://www.aliyundrive.com/s/4yRKtX23mQo
点击链接保存，或者复制本段内容，打开「阿里云盘」APP ，无需下载极速在线查看，视频原画倍速播放。
#### 5.1.2 说明

### 5.2 文件中标签的定义
> 0 代表 germline SNV

> 1 代表 somatic SNV

> 2 代表 Homozygous 纯合子变异

> 3 代表 Heterozygous 杂合子变异

> 4 代表 真的SNV

> 5 代表 假的SNV

### 5.3 CNV检测结果文件格式（csv）
总共包含5列：染色体名称，开始位置，结束位置，扩增（I）或缺失（D），拷贝数的数量

例子如下：
```
chr21,666828,710377,I,10
chr21,4008845,4029066,D,1
chr21,5753247,5839641,I,8
...
```
### 5.4 检测的所有SNV的结果文件格式

### 5.5 检测的所有somatic SNV的结果文件格式



