

import re
import csv
import pandas as pd
import sys
import joblib
from sklearn.preprocessing import MinMaxScaler
import xgboost as xgb

VALID_BASE = set(['a', 'g', 'c', 't', 'A', 'G', 'C', 'T', '.', ','])
BASES = set(['a', 'g', 'c', 't', 'A', 'G', 'C', 'T'])

# 根据预测的SNV从Pileup文件提取每个碱基的信息
def get_base_info(input_pileup, predict_snv):
    outfile = open("result_pileup.txt", 'w')
    pileup_file = open(input_pileup, "r")
    result_file = pd.read_csv(predict_snv)
    num = result_file.index.stop
    for i in range(num):
        result_line = result_file.iloc[i].tolist()
        while True:
            pileup_line = pileup_file.readline()
            if not pileup_line:
                break
            pileup_line = pileup_line.strip().split('\t')
            # position 一样时输出匹配情况和参考碱基
            if result_line[1] == int(pileup_line[1]):
                # 匹配发生单位点变异的reads
                pattern = re.compile(r'[AGCTagct]')
                base_info = pattern.findall(pileup_line[4])

                # 删除indel
                pattern1 = re.compile(r'\+[0-9]+[ACGTNacgtn]+')
                pileup_line[4] = re.sub(pattern1, "", pileup_line[4])
                pattern2 = re.compile(r'\-[0-9]+[ACGTNacgtn]+')
                pileup_line[4] = re.sub(pattern2, "", pileup_line[4])
                # 删除^*
                pattern3 = re.compile(r'\^[A-Z0-9\!\[\]\<\>\:\/\\\?\;\$]')
                pileup_line[4] = re.sub(pattern3, "", pileup_line[4])
                # 删除$!^*, 方便后面提取特征
                baseinfo = list(pileup_line[4])
                new_baseinfo = []
                for i in range(len(baseinfo)):
                    if baseinfo[i] in VALID_BASE:
                        new_baseinfo.append(baseinfo[i])
                # 求变异等位基因频率
                AF = 0
                variant_base = pattern.findall(pileup_line[4])
                # print(variant_base)
                if len(variant_base) > 0:
                    AF = len(variant_base) / int(pileup_line[3])
                # chr, pos, AF,
                outfile.write(result_line[0] + '\t' + str(
                    int(result_line[1])) + '\t' + str(AF) + '\t')
                # RD ,baseinfo, quality, readID
                outfile.write(
                    pileup_line[3] + '\t' + ''.join(new_baseinfo) + '\t' + pileup_line[5] + '\t' + pileup_line[6] + '\n')
                break
            elif int(pileup_line[1]) > result_line[1]:
                break
    outfile.close()

# 计算r^2
def get_rsquare(pa, pb, pab):
    r2 = 0
    if (pa * (1 - pa) * pb * (1 - pb)) != 0:
        r2 = (pab - pa * pb) ** 2 / (pa * (1 - pa) * pb * (1 - pb))
        r2 = float(format(r2, '.5f'))
    return r2


# 对检测的SNV 提取LD和其它特征
def get_ld_feature(source_file, target_file):
    snv_infos = []
    result_file = open(source_file, "r")
    while True:
        line = result_file.readline()
        if not line:
            break
        line = line.strip().split("\t")
        snv_infos.append(line)
    result_file.close()

    out_file = open(target_file, "w", newline="")
    writer = csv.writer(out_file)
    num = len(snv_infos)
    print(num)
    for i in range(num):
        # 前两个
        if i < 2:
            snv1 = snv_infos[i]
            snv2 = snv_infos[i + 1]
            snv3 = snv_infos[i + 2]

            snv1_reads = snv1[6].split(',')
            snv2_reads = snv2[6].split(',')
            snv3_reads = snv3[6].split(',')

            snv1_variant_info = list(snv1[4])
            snv2_variant_info = list(snv2[4])
            snv3_variant_info = list(snv3[4])

            snv1_variant_reads = set()
            snv2_variant_reads = set()
            snv3_variant_reads = set()

            for j in range(len(snv1_variant_info)):
                if snv1_variant_info[j] in BASES:
                    snv1_variant_reads.add(snv1_reads[j])
            for j in range(len(snv2_variant_info)):
                if snv2_variant_info[j] in BASES:
                    snv2_variant_reads.add(snv2_reads[j])
            for j in range(len(snv3_variant_info)):
                if snv3_variant_info[j] in BASES:
                    snv3_variant_reads.add(snv3_reads[j])

            snv1_reads = set(snv1_reads)
            snv2_reads = set(snv2_reads)
            snv3_reads = set(snv3_reads)

            same_all_read1 = snv1_reads.intersection(snv2_reads)
            same_all_read1_num = len(same_all_read1)
            same_variant_read1_num = len(snv1_variant_reads.intersection(snv2_variant_reads))

            same_all_read2 = snv1_reads.intersection(snv3_reads)
            same_all_read2_num = len(same_all_read2)
            same_variant_read2_num = len(snv1_variant_reads.intersection(snv3_variant_reads))

            r_square_1 = 0
            r_square_2 = 0
            if same_all_read1_num != 0:
                pab = same_variant_read1_num / same_all_read1_num
                pa = len(snv1_variant_reads.intersection(same_all_read1)) / same_all_read1_num
                pb = len(snv2_variant_reads.intersection(same_all_read1)) / same_all_read1_num
                r_square_1 = get_rsquare(pa, pb, pab)
            if same_all_read2_num != 0:
                pab = same_variant_read2_num / same_all_read2_num
                pa = len(snv1_variant_reads.intersection(same_all_read2)) / same_all_read2_num
                pb = len(snv3_variant_reads.intersection(same_all_read2)) / same_all_read2_num
                r_square_1 = get_rsquare(pa, pb, pab)
            r_squrare = (r_square_1 + r_square_2) * 2

            # chr, pos, depth, AF, LD,quality
            row = [snv_infos[i][0], snv_infos[i][1], snv_infos[i][3], snv_infos[i][2], r_squrare]
            writer.writerow(row)

            # 倒数2个
        elif i >= (len(snv_infos) - 2):
            snv3 = snv_infos[i - 2]
            snv2 = snv_infos[i - 1]
            snv1 = snv_infos[i]

            snv1_reads = snv1[6].split(',')
            snv2_reads = snv2[6].split(',')
            snv3_reads = snv3[6].split(',')

            snv1_variant_info = list(snv1[4])
            snv2_variant_info = list(snv2[4])
            snv3_variant_info = list(snv3[4])

            snv1_variant_reads = set()
            snv2_variant_reads = set()
            snv3_variant_reads = set()

            for j in range(len(snv1_variant_info)):
                if snv1_variant_info[j] in BASES:
                    snv1_variant_reads.add(snv1_reads[j])
            for j in range(len(snv2_variant_info)):
                if snv2_variant_info[j] in BASES:
                    snv2_variant_reads.add(snv2_reads[j])
            for j in range(len(snv3_variant_info)):
                if snv3_variant_info[j] in BASES:
                    snv3_variant_reads.add(snv3_reads[j])

            snv1_reads = set(snv1_reads)
            snv2_reads = set(snv2_reads)
            snv3_reads = set(snv3_reads)

            same_all_read1 = snv1_reads.intersection(snv2_reads)
            same_all_read1_num = len(same_all_read1)
            same_variant_read1_num = len(snv1_variant_reads.intersection(snv2_variant_reads))

            same_all_read2 = snv1_reads.intersection(snv3_reads)
            same_all_read2_num = len(same_all_read2)
            same_variant_read2_num = len(snv1_variant_reads.intersection(snv3_variant_reads))

            r_square_1 = 0
            r_square_2 = 0
            if same_all_read1_num != 0:
                pab = same_variant_read1_num / same_all_read1_num
                pa = len(snv1_variant_reads.intersection(same_all_read1)) / same_all_read1_num
                pb = len(snv2_variant_reads.intersection(same_all_read1)) / same_all_read1_num
                r_square_1 = get_rsquare(pa, pb, pab)
            if same_all_read2_num != 0:
                pab = same_variant_read2_num / same_all_read2_num
                pa = len(snv1_variant_reads.intersection(same_all_read2)) / same_all_read2_num
                pb = len(snv3_variant_reads.intersection(same_all_read2)) / same_all_read2_num
                r_square_1 = get_rsquare(pa, pb, pab)
            r_squrare = (r_square_1 + r_square_2) * 2
            # chr, pos, depth, AF, LD,quality
            row = [snv_infos[i][0], snv_infos[i][1], snv_infos[i][3], snv_infos[i][2], r_squrare]
            writer.writerow(row)

        else:
            # 中间其余的snv位点
            snv_2 = snv_infos[i - 2]
            snv_1 = snv_infos[i - 1]
            snv = snv_infos[i]
            snv1 = snv_infos[i + 1]
            snv2 = snv_infos[i + 2]

            snv_2_reads = snv_2[6].split(',')
            snv_1_reads = snv_1[6].split(',')
            snv_reads = snv[6].split(',')
            snv1_reads = snv1[6].split(',')
            snv2_reads = snv2[6].split(',')

            snv_2_variant_info = list(snv_2[4])
            snv_1_variant_info = list(snv_1[4])
            snv_variant_info = list(snv[4])
            snv1_variant_info = list(snv1[4])
            snv2_variant_info = list(snv2[4])

            snv_2_variant_reads = set()
            snv_1_variant_reads = set()
            snv_variant_reads = set()
            snv1_variant_reads = set()
            snv2_variant_reads = set()

            for j in range(len(snv_2_variant_info)):
                if snv_2_variant_info[j] in BASES:
                    snv_2_variant_reads.add(snv_2_reads[j])
            for j in range(len(snv_1_variant_info)):
                if snv_1_variant_info[j] in BASES:
                    snv_1_variant_reads.add(snv_1_reads[j])
            for j in range(len(snv_variant_info)):
                if snv_variant_info[j] in BASES:
                    snv_variant_reads.add(snv_reads[j])
            for j in range(len(snv1_variant_info)):
                if snv1_variant_info[j] in BASES:
                    snv1_variant_reads.add(snv1_reads[j])
            for j in range(len(snv2_variant_info)):
                if snv2_variant_info[j] in BASES:
                    snv2_variant_reads.add(snv2_reads[j])

            snv_2_reads = set(snv_2_reads)
            snv_1_reads = set(snv_1_reads)
            snv_reads = set(snv_reads)
            snv1_reads = set(snv1_reads)
            snv2_reads = set(snv2_reads)
            # 计算r^2

            # 当前位点与前面倒数第 2个点
            same_all_reads_2 = snv_reads.intersection(snv_2_reads)
            same_all_reads_2_num = len(same_all_reads_2)
            same_variant_reads_2_num = len(snv_variant_reads.intersection(snv_2_variant_reads))  # 两个点同时变的num

            # 当前位点与前面倒数第 1个点
            same_all_reads_1 = snv_reads.intersection(snv_1_reads)
            same_all_reads_1_num = len(same_all_reads_1)
            same_variant_reads_1_num = len(snv_variant_reads.intersection(snv_1_variant_reads))

            # 当前位点与后面第 1 个点
            same_all_reads1 = snv_reads.intersection(snv1_reads)
            same_all_reads1_num = len(same_all_reads1)
            same_variant_reads1_num = len(snv_variant_reads.intersection(snv1_variant_reads))

            # 当前位点与后面第 2 个点
            same_all_reads2 = snv_reads.intersection(snv2_reads)
            same_all_reads2_num = len(same_all_reads2)
            same_variant_reads2_num = len(snv_variant_reads.intersection(snv2_variant_reads))

            r_square_1 = 0
            r_square_2 = 0
            r_square_3 = 0
            r_square_4 = 0
            if same_all_reads_2_num != 0:
                pab = same_variant_reads_2_num / same_all_reads_2_num
                pa = len(snv_2_variant_reads.intersection(same_all_reads_2)) / same_all_reads_2_num
                pb = len(snv_variant_reads.intersection(same_all_reads_2)) / same_all_reads_2_num
                r_square_1 = get_rsquare(pa, pb, pab)
            if same_all_reads_1_num != 0:
                pab = same_variant_reads_1_num / same_all_reads_1_num
                pa = len(snv_1_variant_reads.intersection(same_all_reads_1)) / same_all_reads_1_num
                pb = len(snv_variant_reads.intersection(same_all_reads_1)) / same_all_reads_1_num
                r_square_2 = get_rsquare(pa, pb, pab)
            if same_all_reads1_num != 0:
                pab = same_variant_reads1_num / same_all_reads1_num
                pa = len(snv1_variant_reads.intersection(same_all_reads1)) / same_all_reads1_num
                pb = len(snv_variant_reads.intersection(same_all_reads1)) / same_all_reads1_num
                r_square_3 = get_rsquare(pa, pb, pab)
            if same_all_reads2_num != 0:
                pab = same_variant_reads2_num / same_all_reads2_num
                pa = len(snv2_variant_reads.intersection(same_all_reads2)) / same_all_reads2_num
                pb = len(snv_variant_reads.intersection(same_all_reads2)) / same_all_reads2_num
                r_square_4 = get_rsquare(pa, pb, pab)

            r_squrare = r_square_1 + r_square_2 + r_square_3 + r_square_4
            # chr, pos, depth, AF, LD,quality
            row = [snv_infos[i][0], snv_infos[i][1], snv_infos[i][3], snv_infos[i][2], r_squrare]
            writer.writerow(row)
    out_file.close()

# snv predict result
def get_cn_feature(target_file, source_file, cnv_file):
    feature_file = pd.read_csv(source_file, header=None)
    cnv_file = pd.read_csv(cnv_file, header=None)
    final_file = open(target_file, "w", newline="")
    writer = csv.writer(final_file)
    row = ['Chromosome', 'Position','ReadDepth','AlleleFrequency','LinkageDisequilibrium','CopyNumber']
    writer.writerow(row)
    
    num = feature_file.index.stop
    cnv_num = cnv_file.index.stop
    cnv_index = 0
    cnv_info = cnv_file.iloc[cnv_index].tolist()
    cnv_start = cnv_info[1]
    cnv_end = cnv_info[2]
    result = []
    for i in range(num):
        snv = feature_file.iloc[i].tolist()
        if snv[1] < cnv_start:
            snv.append(2)
            result.append(snv)
            writer.writerow(snv)
        elif cnv_start <= snv[1] < cnv_end:
            snv.append(cnv_info[4])
            result.append(snv)
            writer.writerow(snv)
        elif snv[1] >= cnv_end:
            snv.append(2)
            result.append(snv)
            writer.writerow(snv)
            cnv_index += 1
            if cnv_index < cnv_num:
                cnv_info = cnv_file.iloc[cnv_index].tolist()
                cnv_start = cnv_info[1]
                cnv_end = cnv_info[2]
def get_features():
    get_ld_feature(source_file, target_file)
    get_cn_feature(target_file, source_file, cnv_file)

def predict_somatic_snv(feature_file):
    # input predict file
    all_snv_file = pd.read_csv(feature_file)
    all_snv_num = all_snv_file.index.stop

    y_test = all_snv_file.iloc[:, 2:6][0:all_snv_num]  # 2,3,4,5为RD, AF R^2 CN
  
    xgb_model = joblib.load("./model/LDSSNV_distinguish_Single.m")
    ypred = xgb_model.predict(y_test)

    # output predict result
    result_file = open("somatic_snv.csv", "w", newline="")
    writer = csv.writer(result_file)
    row = []
    for i in range(all_snv_num):
        line = all_snv_file.iloc[i].tolist()
        if ypred[i] == 1.0:
            row = [line[0], line[1], ypred[i]]  # position, (0, 1)germline or somatic
            writer.writerow(row)
    result_file.close()

# step1
get_bases_info("tumor.pileup", "out_part_snp.csv")

# step2
get_ld_param("result_pileup.txt", "predict_somatic_feature_part.csv")
# step3
get_cnv_feature()
# step4
xgbSomatic_predict("predict_somatic_feature.csv")



