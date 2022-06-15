'''
Author: cwx
Date: 2021-11-17 11:46:20
LastEditTime: 2022-05-16 17:20:32
Description: 多个样本间的种系SNP具有连锁不平衡性
FilePath: /paper_snv/gen_snp/genSnpSample_final.py
'''
#某个block所在的het位点全变在第 1 条染色体， 或者变在第 2 条，纯粹的样本生成
import numpy as np
import pandas as pd
import math
import random
from common import common_outfile


SAMPLE_NUM = 50  # 样本数
CHROM_NUM = 100  # 染色体条数
MAX_POSITION = 46709983  # 碱基数量
PA_PB_RANGE = 1 / SAMPLE_NUM
PA_PB_RANGE_2 = 1 / CHROM_NUM


'''

description: 获取0-20kb, 均值为10kb, 符合正态分布的随机数(LD块的长度)
param {*}
return {*} int
'''
def getBlockLen():
    block = np.random.normal(10000, 10000)
    block = round(math.fabs(block))
    return block
'''
description: 获取a-b之间的随机数，保留三位小数
param {*} a min = 0
param {*} b max = 1
return {*} 
'''
def getRandomByRange(a, b):
    rd = random.uniform(a, b)
    return float(format(rd, '.5f'))


def getHomRandomByRange(a, b):
    a_min = float(format((1 / SAMPLE_NUM), '.5f'))
    if a <= a_min:
        a = a_min
    if b > 0.5:
        b = 0.5
    return getRandomByRange(a, b)
def getHetRandomByRange(a, b):
    a_min = float(format((1 / CHROM_NUM), '.5f'))
    if a > 0.5 and b > 0.5:
        a = a / 2
        b = b / 2

    if a > a_min and a <= 0.5 and b > 0.5:
        b = 0.5
    elif a <= a_min and b > 0.5:
        a = a_min
        b = 0.5
    elif a <= a_min and b > a_min and b <= 0.5:
        a = a_min
    return getRandomByRange(a, b)
'''
description: 根据Pa, Pb 计算p(ab)，获取两个位点共同变异的样本个数
param {*} Pa
param {*} Pb
return {*}
'''
def getPab(Pa, Pb):
    a = Pa * Pb
    b = min(Pa, Pb)
    a = 0.1 * a + 0.9 * b
    return getRandomByRange(a, b)


def genSnpSample(all_snp_file):
    outLDBlock = open("LDBlock_len.txt", 'w')

    block = getBlockLen()
    outLDBlock.write(str(0) + '\t')
    outLDBlock.write(str(block) + '\t')
    outLDBlock.write(str(block) + '\n')

    results_1 = []  # 所有样本的第一条染色体，snp记录
    results_2 = []  # 所有样本的第二条染色体，snp记录
    for i in range(SAMPLE_NUM):
        results_1.append([])
        results_2.append([])

    all_sample = set()  # 所有样本编号集合
    for i in range(SAMPLE_NUM):
        all_sample.add(i)

    file_snp_2 = pd.read_csv(all_snp_file, header=None) 
    snp_num = file_snp_2.index.stop  # snp位点总数

    position = block  # 碱基序列当前位置
    block_first_pos = True
    prev_snp_prob = 0.5
    prev_snp_samples = []
    prev_snp_is_het = False
    # 保存前一个snp是het的时候变异所在的样本编号
    prev_snp_het_1 = []  # 第一条染色体
    prev_snp_het_2 = []  # 第二条染色体
    # ! 每到一个block，随机生成一次，大于0.5, 则该block所在的het位点全变在第 1 条染色体， 否则变在第 2 条
    block_single_het = random.random()

   

    for i in range(snp_num):
        line_snp = file_snp_2.iloc[i].tolist()
        chr_pos = line_snp[1]
        
        # 当前snp位点的位置小于当前块的右边界
        if chr_pos <= position:
            if block_first_pos:
                block_first_pos = False
                # 该位点是hom
                if line_snp[5] == 2:
                    cur_snp_prob = getHomRandomByRange(0, 0.5)
                    num = math.floor((cur_snp_prob * CHROM_NUM / 2))  # 样本数，成对的染色体
                    cur_snp_samples = random.sample(all_sample, num)
                    for j in range(num):
                        results_1[cur_snp_samples[j]].append(line_snp)
                        results_2[cur_snp_samples[j]].append(line_snp)
                    prev_snp_prob = cur_snp_prob
                    prev_snp_samples = cur_snp_samples
                   
                # 该位点是het
                else:
                    cur_snp_prob = getHetRandomByRange(0.25, 0.5)
                    num = math.floor((cur_snp_prob * CHROM_NUM))
                    cur_snp_samples = random.sample(all_sample, num);
                    print('first_pos, het: pa:{0}, samples:{1}'.format(cur_snp_prob, len(cur_snp_samples)))
                    if block_single_het > 0.5:
                        for j in range(num):
                            results_1[cur_snp_samples[j]].append(line_snp)
                            prev_snp_het_1.append(cur_snp_samples[j])
                    else:
                        for j in range(num):
                            results_2[cur_snp_samples[j]].append(line_snp)
                            prev_snp_het_2.append(cur_snp_samples[j])
                    prev_snp_is_het = True
                    prev_snp_prob = cur_snp_prob
                    prev_snp_samples = cur_snp_samples
            else:
                # 上一个位点是het
                if prev_snp_is_het:
                    # 该位点是hom
                    if line_snp[5] == 2:
                        cur_snp_samples = []
                        cur_snp_prob = getHomRandomByRange(2 * prev_snp_prob - PA_PB_RANGE,
                                                        2 * prev_snp_prob + PA_PB_RANGE)
                        left_value = 0.1 * cur_snp_prob * prev_snp_prob + 0.9 * min(0.5 * cur_snp_prob, prev_snp_prob)
                        Pab = getRandomByRange(left_value, min(0.5*cur_snp_prob, prev_snp_prob))
                        adjacentSnpCommonSampleNum = math.floor(Pab * CHROM_NUM)
                        curSnpAllSampleNum = math.floor(cur_snp_prob * SAMPLE_NUM)

                        # 与上一个snp所在样本的交集中选取adjacentSnpCommonSampleNum
                        tempList = random.sample(prev_snp_samples, adjacentSnpCommonSampleNum)
                        for j in range(adjacentSnpCommonSampleNum):
                            results_1[tempList[j]].append(line_snp)
                            results_2[tempList[j]].append(line_snp)
                           
                        cur_snp_samples = tempList
                        # 剩下的从上一个snp所在样本的差集中选取
                        remaing_samples = all_sample.symmetric_difference(set(prev_snp_samples))  # 取差集
                        curSnpSingleSampleNum = curSnpAllSampleNum - adjacentSnpCommonSampleNum
                      
                        # ! curSnpSingleSampleNum 大于 剩余样本数量
                        if curSnpSingleSampleNum > len(remaing_samples):
                            curSnpSingleSampleNum = len(remaing_samples)
                            cur_snp_prob = (curSnpSingleSampleNum + adjacentSnpCommonSampleNum) / CHROM_NUM
                       
                        tempList = random.sample(remaing_samples, curSnpSingleSampleNum);
                        for j in range(curSnpSingleSampleNum):
                            results_1[tempList[j]].append(line_snp)
                            results_2[tempList[j]].append(line_snp)
                            
                        cur_snp_samples = set(cur_snp_samples + tempList)
                        prev_snp_is_het = False
                        prev_snp_prob = cur_snp_prob;
                        prev_snp_samples = cur_snp_samples;
                       
                    # 该位点是het
                    else:
                        temp_snp_het_1 = []
                        temp_snp_het_2 = []
                        cur_snp_samples = []
                        cur_snp_prob = getHetRandomByRange(prev_snp_prob - PA_PB_RANGE_2, prev_snp_prob + PA_PB_RANGE_2)
                        Pab = getPab(prev_snp_prob, cur_snp_prob)
                        adjacentSnpCommonSampleNum = math.floor(Pab * CHROM_NUM)
                        curSnpAllSampleNum = math.floor(cur_snp_prob * CHROM_NUM)
                       
                        # 与上一个snp所在样本的交集中选取adjacentSnpCommonSampleNum
                        tempList = random.sample(prev_snp_samples, adjacentSnpCommonSampleNum)
                        if block_single_het > 0.5:
                            for j in range(adjacentSnpCommonSampleNum):
                                results_1[tempList[j]].append(line_snp)
                                temp_snp_het_1.append(tempList[j])
                        else:
                            for j in range(adjacentSnpCommonSampleNum):
                                results_2[tempList[j]].append(line_snp)
                                temp_snp_het_2.append(tempList[j])

                        cur_snp_samples = tempList
                        # 剩下的从上一个snp所在样本的差集中选取
                        remaing_samples = all_sample.symmetric_difference(set(prev_snp_samples))  # 取差集
                        curSnpSingleSampleNum = curSnpAllSampleNum - adjacentSnpCommonSampleNum;
                       
                        # ! curSnpSingleSampleNum 大于 剩余样本数量
                        if curSnpSingleSampleNum > len(remaing_samples):
                            curSnpSingleSampleNum = len(remaing_samples)
                            cur_snp_prob = (curSnpSingleSampleNum + adjacentSnpCommonSampleNum) / CHROM_NUM
                        
                        tempList = random.sample(remaing_samples, curSnpSingleSampleNum)
                        if block_single_het > 0.5:
                            for j in range(curSnpSingleSampleNum):
                                results_1[tempList[j]].append(line_snp)
                                temp_snp_het_1.append(tempList[j])
                        else:
                            for j in range(curSnpSingleSampleNum):
                                results_2[tempList[j]].append(line_snp)
                                temp_snp_het_2.append(tempList[j])
                        cur_snp_samples = set(cur_snp_samples + tempList)
                        prev_snp_het_1 = temp_snp_het_1
                        prev_snp_het_2 = temp_snp_het_2
                        prev_snp_is_het = True
                        prev_snp_prob = cur_snp_prob
                        prev_snp_samples = cur_snp_samples
                # 上一个位点hom
                else:
                    # 该位点是hom
                    if line_snp[5] == 2:

                        cur_snp_prob = getHomRandomByRange(prev_snp_prob - PA_PB_RANGE, prev_snp_prob + PA_PB_RANGE)
                        Pab = getPab(prev_snp_prob, cur_snp_prob)
                        adjacentSnpCommonSampleNum = math.floor(Pab * SAMPLE_NUM)
                        curSnpAllSampleNum = math.floor(cur_snp_prob * SAMPLE_NUM)
                      
                        # 与上一个snp所在样本的交集中选取adjacentSnpCommonSampleNum
                        tempList = random.sample(prev_snp_samples, adjacentSnpCommonSampleNum)
                        for j in range(adjacentSnpCommonSampleNum):
                            results_1[tempList[j]].append(line_snp)
                            results_2[tempList[j]].append(line_snp)
                        cur_snp_samples = tempList

                        # 剩下的从上一个snp所在样本的差集中选取
                        remaing_samples = all_sample.symmetric_difference(set(prev_snp_samples))  # 取差集
                        curSnpSingleSampleNum = curSnpAllSampleNum - adjacentSnpCommonSampleNum
                       
                        # ! curSnpSingleSampleNum 大于 剩余样本数量
                        if curSnpSingleSampleNum > len(remaing_samples):
                            curSnpSingleSampleNum = len(remaing_samples)
                            cur_snp_prob = (curSnpSingleSampleNum + adjacentSnpCommonSampleNum) / CHROM_NUM
                        
                        tempList = random.sample(remaing_samples, curSnpSingleSampleNum)
                        for j in range(curSnpSingleSampleNum):
                            results_1[tempList[j]].append(line_snp)
                            results_2[tempList[j]].append(line_snp)
                           

                        cur_snp_samples = set(cur_snp_samples + tempList)
                        prev_snp_is_het = False
                        prev_snp_prob = cur_snp_prob
                        prev_snp_samples = cur_snp_samples
                        
                    # 该位点是het
                    else:
                        temp_snp_het_1 = []
                        temp_snp_het_2 = []
                        cur_snp_samples = []
                        cur_snp_prob = getHetRandomByRange((prev_snp_prob - PA_PB_RANGE) * 0.5, 0.5 * (prev_snp_prob + PA_PB_RANGE))
                        left_value = 0.1 * cur_snp_prob * prev_snp_prob + 0.9 * min(cur_snp_prob, 0.5 * prev_snp_prob)
                        Pab = getRandomByRange(left_value, min(cur_snp_prob, 0.5 * prev_snp_prob))
                        adjacentSnpCommonSampleNum = math.floor(Pab * CHROM_NUM)
                        curSnpAllSampleNum = math.floor(cur_snp_prob * CHROM_NUM)
                        
                        # 与上一个snp所在样本的交集中选取adjacentSnpCommonSampleNum
                        tempList = random.sample(prev_snp_samples, adjacentSnpCommonSampleNum)
                        if block_single_het > 0.5:
                            for j in range(adjacentSnpCommonSampleNum):
                                results_1[tempList[j]].append(line_snp)
                                temp_snp_het_1.append(tempList[j])
                        else:
                            for j in range(adjacentSnpCommonSampleNum):
                                results_2[tempList[j]].append(line_snp)
                                temp_snp_het_2.append(tempList[j])

                        cur_snp_samples = tempList
                        # 剩下的从上一个snp所在样本的差集中选取
                        remaing_samples = all_sample.symmetric_difference(set(prev_snp_samples))  # 取差集
                        curSnpSingleSampleNum = curSnpAllSampleNum - adjacentSnpCommonSampleNum;

                        # ! curSnpSingleSampleNum 大于 剩余样本数量
                        if curSnpSingleSampleNum > len(remaing_samples):
                            curSnpSingleSampleNum = len(remaing_samples)
                            cur_snp_prob = (curSnpSingleSampleNum + adjacentSnpCommonSampleNum) / CHROM_NUM
                        
                        tempList = random.sample(remaing_samples, curSnpSingleSampleNum);
                        if block_single_het > 0.5:
                            for j in range(curSnpSingleSampleNum):
                                results_1[tempList[j]].append(line_snp)
                                temp_snp_het_1.append(tempList[j])
                        else:
                            for j in range(curSnpSingleSampleNum):
                                results_2[tempList[j]].append(line_snp)
                                temp_snp_het_2.append(tempList[j])
                                
                        cur_snp_samples = set(cur_snp_samples + tempList)
                        prev_snp_het_1 = temp_snp_het_1
                        prev_snp_het_2 = temp_snp_het_2
                        prev_snp_is_het = True
                        prev_snp_prob = cur_snp_prob
                        prev_snp_samples = cur_snp_samples
                       
        elif position < MAX_POSITION:
            block_single_het = random.random()
            # ! 生成该块内R^2的详细信息
            # 块第一个位点 该位点是hom
            if line_snp[5] == 2:
                cur_snp_prob = getHomRandomByRange(0, 0.5)
                num = math.floor((cur_snp_prob * CHROM_NUM / 2))  # 样本数，成对的染色体
                cur_snp_samples = random.sample(all_sample, num)
                for j in range(num):
                    results_1[cur_snp_samples[j]].append(line_snp)
                    results_2[cur_snp_samples[j]].append(line_snp)

                prev_snp_prob = cur_snp_prob
                prev_snp_samples = cur_snp_samples
                prev_snp_is_het = False

            # 该位点是het
            else:
                prev_snp_het_1 = []
                prev_snp_het_2 = []
                cur_snp_prob = getHetRandomByRange(0, 0.5)
                num = math.floor((cur_snp_prob * CHROM_NUM))
                cur_snp_samples = random.sample(all_sample, num)
                # print('first_pos, het: pa:{0}, samples:{1}'.format(cur_snp_prob, len(cur_snp_samples)))
                if block_single_het > 0.5:
                    for j in range(num):
                        results_1[cur_snp_samples[j]].append(line_snp)
                        prev_snp_het_1.append(cur_snp_samples[j])
                else:
                    for j in range(num):
                        results_2[cur_snp_samples[j]].append(line_snp)
                        prev_snp_het_2.append(cur_snp_samples[j])

                prev_snp_is_het = True
                prev_snp_prob = cur_snp_prob
                prev_snp_samples = cur_snp_samples

            block = getBlockLen()
            outLDBlock.write(str(position) + '\t')
            outLDBlock.write(str(block) + '\t')
            position = position + block  # 碱基序列当前位置
            outLDBlock.write(str(position) + '\n')
        else:
            print("区块划分完成！")
    outLDBlock.close()
    common_outfile.outAllSnpSample(results_1, '1', SAMPLE_NUM)
    common_outfile.outAllSnpSample(results_2, '2', SAMPLE_NUM)


genSnpSample("chr21_snp1.csv")

