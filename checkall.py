#!/usr/bin/env python3
#rsid	chromosome	position	genotype

#rivername = '111-1726-4025_202011202334'
#ryanname  = '111-1727-7094_202011232301'
#kellyname = '111-1727-9111_202011232302'

rivername = '111-1726-4025_202011231011'
ryanname  = '111-1727-7094_202011231611'
kellyname = '111-1727-9111_202011231612'

import os
import numpy as np
#from numba import njit

def readData(filename, forceReload = False):
    npyName = filename + '.npy'
    txtName = filename + '.txt'
    if os.path.isfile(npyName) and forceReload == False:
        arrayData = np.load(npyName)
        print('loaded from npy',npyName)
    else:
        f = open(txtName,'r')
        lines = f.readlines()
        f.close()

        listData = []
        for line in lines:
            if line.startswith('#'):
                continue
            info = line.split()
            if len(info) != 4:
                continue
            listData.append(info)
        arrayData = np.array(listData)
        np.save(npyName, arrayData)
        print('loaded from txt',txtName)
    return arrayData

def statChromosome(arrayData):
    dictChromosome = {}
    for i in range(len(arrayData)):
        ch = arrayData[i][1]
        if ch not in dictChromosome:
            dictChromosome[ch] = 1
        else:
            dictChromosome[ch] += 1
    for ch in dictChromosome:
        print(ch, dictChromosome[ch])
    return

def getXYM(arrayData):
    dictX = {}
    dictY = {}
    dictM = {}
    for i in range(len(arrayData)):
        if arrayData[i][1] == 'X':
            dictX[arrayData[i][0]] = arrayData[i][3]
            continue
        if arrayData[i][1] == 'Y':
            dictY[arrayData[i][0]] = arrayData[i][3]
            continue
        if arrayData[i][1] == 'MT':
            dictM[arrayData[i][0]] = arrayData[i][3]
            continue
    return dictX, dictY, dictM

def checkDiff(dict_1, dict_2, CH):
    listDiff = []
    listDiff_0 = []
    if CH == 'Y' or CH == 'MT':
        for rsid in dict_1:
            if dict_1[rsid] != dict_2[rsid] or (dict_1[rsid] == '--' and dict_2[rsid] == '--'):
                listDiff.append((rsid, CH, dict_1[rsid], dict_2[rsid]))
                if dict_1[rsid] != '--' and dict_2[rsid] != '--':
                    listDiff_0.append((rsid, CH, dict_1[rsid], dict_2[rsid]))
    if CH == 'X':
        for rsid in dict_1:
            if dict_2[rsid][0] != dict_2[rsid][1] or (dict_1[rsid] == '--' and dict_2[rsid] == '--'):
                listDiff.append((rsid, CH, dict_1[rsid], dict_2[rsid]))
                if dict_1[rsid] != '--' and dict_2[rsid] != '--' and dict_2[rsid][0] not in dict_1[rsid]:
                    listDiff_0.append((rsid, CH, dict_1[rsid], dict_2[rsid]))
    return listDiff, listDiff_0

def checkMiss(arrayData, female=False):
    listMiss = []
    numTotal = 0
    for i in range(len(arrayData)):
        gt = arrayData[i][3]
        ch = arrayData[i][1]
        if not (female == True and ch == 'Y'):
            numTotal += 1
            if gt == '--':
                listMiss.append(arrayData[i]) 
    return numTotal, listMiss

def checkX(dictData):
    listWrong = []
    for rsid in dictData:
        gt = dictData[rsid]
        if gt[0] != gt[1]:
            listWrong.append((rsid, 'X', gt))
    return len(dictData), listWrong

def transArray2Dict(arrayData):
    dictData = {}
    for i in range(len(arrayData)):
        dictData[arrayData[i][0]] = arrayData[i]
    return dictData

def checkMutation(arrayD1, arrayD2, arrayD3, dictD1, dictD2, dictD3):
    listMutation = []
    countMissing = 0
    countCheck = 0
    if not (len(dictD1) == len(dictD2) and len(dictD1) == len(dictD3)):
        print('Data Len Error!')
        return countMissing, countCheck, listMutation
    print('checkMutation len:', len(dictD1))
    for rsid in dictD1:
        if dictD1[rsid][1] == 'Y' or dictD1[rsid][1] == 'X' or dictD1[rsid][1] == 'MT' :
            continue
        countCheck += 1
        if dictD1[rsid][3] == '--' or dictD2[rsid][3] == '--' or dictD3[rsid][3] == '--':
            countMissing += 1
            continue
        if not ((dictD1[rsid][3][0] in dictD2[rsid][3] and dictD1[rsid][3][1] in dictD3[rsid][3]) or (dictD1[rsid][3][0] in dictD3[rsid][3] and dictD1[rsid][3][1] in dictD2[rsid][3])):
            listMutation.append((rsid, dictD1[rsid][1], dictD1[rsid][3], dictD2[rsid][3], dictD3[rsid][3]))
    return countMissing, countCheck, listMutation

def main():
    riverdata = readData(rivername)
    print('riverdata.shape', riverdata.shape)

    ryandata = readData(ryanname)
    print('ryandata.shape ', ryandata.shape)

    kellydata = readData(kellyname)
    print('kellydata.shape', kellydata.shape)

    dictDataRiver = transArray2Dict(riverdata)
    dictDataRyan  = transArray2Dict(ryandata)
    dictDataKelly = transArray2Dict(kellydata)

    #statChromosome(riverdata)

    dictX_River, dictY_River, dictM_River = getXYM(riverdata)
    print('len dictX_River',len(dictX_River))
    print('len dictY_River',len(dictY_River))
    print('len dictM_River',len(dictM_River))

    dictX_Ryan, dictY_Ryan, dictM_Ryan = getXYM(ryandata)
    print('len dictX_Ryan ',len(dictX_Ryan))
    print('len dictY_Ryan ',len(dictY_Ryan))
    print('len dictM_Ryan ',len(dictM_Ryan))

    dictX_Kelly, dictY_Kelly, dictM_Kelly = getXYM(kellydata)
    print('len dictX_Kelly',len(dictX_Kelly))
    print('len dictY_Kelly',len(dictY_Kelly))
    print('len dictM_Kelly',len(dictM_Kelly))

    listDiff_River_Ryan, listDiff_River_Ryan_0 = checkDiff(dictY_River, dictY_Ryan, 'Y')
    print('River_Ryan Y  Diff:    {}/{}={}%'.format(len(listDiff_River_Ryan), len(dictY_River), round(100*len(listDiff_River_Ryan)/len(dictY_River),2)))
    print('River_Ryan Y  Diff_0:  {}/{}={}%'.format(len(listDiff_River_Ryan_0), len(dictY_River), round(100*len(listDiff_River_Ryan_0)/len(dictY_River),2)))

    listDiff_Kelly_Ryan, listDiff_Kelly_Ryan_0 = checkDiff(dictM_Kelly, dictM_Ryan, 'MT')
    print('Kelly_Ryan MT Diff:    {}/{}={}%'.format(len(listDiff_Kelly_Ryan), len(dictM_Kelly), round(100*len(listDiff_Kelly_Ryan)/len(dictM_Kelly),2)))
    print('Kelly_Ryan MT Diff_0:  {}/{}={}%'.format(len(listDiff_Kelly_Ryan_0), len(dictM_Kelly), round(100*len(listDiff_Kelly_Ryan_0)/len(dictM_Kelly),2)))

    listCheck_Kelly_Ryan, listCheck_Kelly_Ryan_0 = checkDiff(dictX_Kelly, dictX_Ryan, 'X')
    print('Kelly_Ryan X  Check:   {}/{}={}%'.format(len(listCheck_Kelly_Ryan), len(dictX_Kelly), round(100*len(listCheck_Kelly_Ryan)/len(dictX_Kelly),2)))
    print('Kelly_Ryan X  Check_0: {}/{}={}%'.format(len(listCheck_Kelly_Ryan_0), len(dictX_Kelly), round(100*len(listCheck_Kelly_Ryan_0)/len(dictX_Kelly),2)))

    #print('listDiff_River_Ryan_0')
    #for info in listDiff_River_Ryan_0:
    #    print(info)

    #print('listDiff_Kelly_Ryan_0')
    #for info in listDiff_Kelly_Ryan_0:
    #    print(info)

    #print('listDiff_River_Ryan')
    #for info in listDiff_River_Ryan:
    #    print(info)

    #print('listDiff_Kelly_Ryan')
    #for info in listDiff_Kelly_Ryan:
    #    print(info)

    numberXRiver, listXWRiver = checkX(dictX_River)
    print('River Wrong X: {}/{}={}%'.format(len(listXWRiver), numberXRiver, round(100*len(listXWRiver)/numberXRiver,2)))

    numberXRyan, listXWRyan = checkX(dictX_Ryan)
    print('Ryan  Wrong X: {}/{}={}%'.format(len(listXWRyan), numberXRyan, round(100*len(listXWRyan)/numberXRyan,2)))

    numberTotalRiver, listMissRiver = checkMiss(riverdata)
    print('River Miss: {}/{}={}%'.format(len(listMissRiver), numberTotalRiver, round(100*len(listMissRiver)/numberTotalRiver,2)))

    numberTotalRyan, listMissRyan= checkMiss(ryandata)
    print('Ryan  Miss: {}/{}={}%'.format(len(listMissRyan), numberTotalRyan, round(100*len(listMissRyan)/numberTotalRyan,2)))

    numberTotalKelly, listMissKelly= checkMiss(kellydata, female=True)
    print('Kelly Miss: {}/{}={}%'.format(len(listMissKelly), numberTotalKelly, round(100*len(listMissKelly)/numberTotalKelly,2)))

    numberMissing, numberCheck, listMutation = checkMutation(ryandata, riverdata, kellydata, dictDataRyan, dictDataRiver, dictDataKelly)
    print('Missing: {} , Mutate: {}/{}={}%'.format(numberMissing, len(listMutation), numberCheck, round(100*len(listMutation)/numberCheck,2)))

    #for mutation in listMutation:
    #    print(mutation)

    #print('listCheck_Kelly_Ryan')
    #for i in range(min(10,len(listCheck_Kelly_Ryan))):
    #    print(listCheck_Kelly_Ryan[i])
    #print('listCheck_Kelly_Ryan_0')
    #for i in range(min(10,len(listCheck_Kelly_Ryan_0))):
    #    print(listCheck_Kelly_Ryan_0[i])

    #print('listXWRiver')
    #for i in range(min(10,len(listXWRiver))):
    #    print(listXWRiver[i])
    #print('listXWRyan')
    #for i in range(min(10,len(listXWRyan))):
    #    print(listXWRyan[i])

    return

if __name__ == '__main__':
    main()
