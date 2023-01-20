import re
import pandas as pd
import joblib
import csv
import xgboost as xgb


def get_pre_selected_snv(pileup_file):
    feature_file = open("part.pileup",'w')
    pattern = re.compile(r'[AGCTNagctn]')
    tumor_pileup = open(pileup_file, 'r')
    while True:
        line = tumor_pileup.readline()
        if not line:
            break
        line_list = line.strip().split('\t')
        result = pattern.findall(line_list[4])
        if len(result) > 0:
            feature_file.write(line)
    feature_file.close()

def get_four_feature():

    feature_file = open("detect_snp_feature.csv",'w', newline='')
    writer = csv.writer(feature_file)
    row = ['chrName', 'position','Mismatched reads', "VAF", "Read depth", 'Ave of quality']
    writer.writerow(row) # output column head

    pattern1 = re.compile(r'\^[A-Z0-9\!\[\]\<\>\:\/\\\?\;\$]')
    pattern2 = re.compile(r'\$')
    pattern3 = re.compile(r'\+[0-9]+[ACGTNacgtn]+')
    pattern4 = re.compile(r'\-[0-9]+[ACGTNacgtn]+')
    pattern5 = re.compile(r'\.')
    pattern6 = re.compile(r'\,')
    pattern7 = re.compile(r'[AGCTagct]')
    pattern8 = re.compile(r'\*')

    tumor_pileup = open("part.pileup", 'r')
    while True:
        row = []
        line = tumor_pileup.readline()
        if not line:
            break
        line = line.strip().split('\t')
        row.append(line[0]) # chr
        row.append(line[1]) # position

        result = pattern1.findall(line[4])
        # row.append(len(result))
        line[4] = re.sub(r'\^[A-Z0-9\!\[\]\<\>\:\/\\\?\;\$]', "", line[4])

        result = pattern2.findall(line[4])
        # row.append(len(result))
        line[4] = re.sub(r'\$', "", line[4])

        result = pattern3.findall(line[4])
        # row.append(len(result))
        line[4] = re.sub(pattern3, "", line[4])

        result = pattern4.findall(line[4])
        # row.append(len(result))
        line[4] = re.sub(pattern4, "", line[4])

        result = pattern5.findall(line[4])
        # row.append(len(result))

        result = pattern6.findall(line[4])
        # row.append(len(result))
        
        # feature1: Mismatched reads
        result_mismatch = pattern7.findall(line[4])
        row.append(len(result_mismatch))

        # feature2: Variant Allele Frequency
        AF = 0
        if int(line[3]) != 0:
            AF = len(result_mismatch) / int(line[3])
        row.append(AF)

        # feature3: Read Depth
        RD = int(line[3])
        row.append(RD)

        quality = list(line[5])

        # 错配的质量分数
        # base_info = list(line[4])
        # mis_quality = [] # 记录错配read的质量值
        # for j in range(len(base_info)):
        #     if base_info[j] in BASES:
        #         mis_quality.append(quality[j])

        # mis_quality_value = 0
        # for j in range(len(mis_quality)):
        #     mis_quality_value += (ord(quality[j]) - 33)
        # quality_mean_value = 0
        # if len(mis_quality) != 0:
        #     quality_mean_value = mis_quality_value / len(mis_quality)
        # row.append(quality_mean_value)
        
        # feature4: Ave of the Quality
        quality_value = 0
        for j in range(len(quality)):
            quality_value += (ord(quality[j]) - 33)
        quality_mean_value = quality_value / len(quality)
        row.append(quality_mean_value)

        writer.writerow(row)
    feature_file.close()

# 
def get_cn_feature(target_file, source_file, cnv_result):
    feature_file = pd.read_csv(source_file)
    cnv_file = pd.read_csv(cnv_result, header=None)
    
    final_file = open(target_file, "w", newline="")
    writer = csv.writer(final_file)
    row = ['chrName', 'position','Mismatched reads', "VAF", "Read depth", 'Ave of quality', 'Copy number']
    writer.writerow(row)# output column head
    
    num = feature_file.index.stop
    cnv_num = cnv_file.index.stop
    cnv_index = 0
    cnv_info = cnv_file.iloc[cnv_index].tolist()
    cnv_start = cnv_info[1]
    cnv_end = cnv_info[2]
    for i in range(num):
        snv = feature_file.iloc[i].tolist()
        if snv[1] < cnv_start:
            snv.append(2)
            writer.writerow(snv)
        elif cnv_start <= snv[1] < cnv_end:
            snv.append(cnv_info[4])
           
            writer.writerow(snv)
        elif snv[1] >= cnv_end:
            snv.append(2)     
            writer.writerow(snv)
            cnv_index += 1
            if cnv_index < cnv_num:
                cnv_info = cnv_file.iloc[cnv_index].tolist()
                cnv_start = cnv_info[1]
                cnv_end = cnv_info[2]

def get_features():
    get_four_features()
    get_cnv_number(args[1], "feature.csv", CNV_TRUTH_FILE)

def predict_snv(PREDICT_FILE):
    # input predict file
    all_snv_file = pd.read_csv(PREDICT_FILE)
    all_snv_num = all_snv_file.index.stop
    print(str(all_snv_num))
    y_test = all_snv_file.iloc[:, 2:7][0:all_snv_num]  # [2,3,4,5,6] 5个特征
    y_test = xgb.DMatrix(y_test)
    # predicting
    xgb_model = joblib.load("./model/LDSSNV_detect.m")
    ypred = xgb_model.predict(y_test)

    print(str(len(ypred)))
    # output predict result
    ypred = (ypred >= 0.5)*1
    result_file = open("out_snp.csv", "w", newline="")
    writer = csv.writer(result_file)
    for i in range(all_snv_num):
        line = all_snv_file.iloc[i].tolist()
        if ypred[i] == 0.0: # is snv
            row = [line[0], line[1]]  # position, germline or somatic
            writer.writerow(row)
    result_file.close()


CNV_TRUTH_FILE = "/home/cwx/DATA/bio_exprement/data/cnv_truth.csv"

get_pre_selected_snv()
get_features()
predict_snv()
