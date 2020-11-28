#!/usr/bin/env python3
#rsid	chromosome	position	genotype

#Name = 'Ryan'
Name = 'River'
#Name = 'Kelly'

rivername = '111-1726-4025_202011231011.txt'
ryanname  = '111-1727-7094_202011231611.txt'
kellyname = '111-1727-9111_202011231612.txt'
dictName = {
    'Ryan':ryanname,
    'River':rivername,
    'Kelly':kellyname,
}
filename = dictName[Name].rstrip('.txt')

snpname = 'snpinfo.inf'

import sqlite3
import os
import numpy as np

def loadSnpInfo(snpname):
    dictSNP = {}
    f = open(snpname, 'r', encoding='UTF-8')
    lines = f.readlines()
    f.close()
    head = 0
    key = ''
    rsid = ''
    for line in lines:
        line = line.strip()
        if line.startswith('#'):
            continue
        if line == '':
            head = 0
            key = ''
            rsid = ''
            continue
        if head == 0:
            key = line.strip()
            if key not in dictSNP:
                dictSNP[key] = {}
            head += 1
            continue
        if head == 1:
            ch = ' '
            if ':' in line: ch = ':'
            if '：' in line: ch = '：'
            info = line.split(ch)
            rsid = info[-1]
            if rsid in dictSNP[key]:
                print('warning reduplicated rsid', rsid)
            dictSNP[key][rsid] = []
            head += 1
            continue
        if head == 2:
            ch = ' '
            if ':' in line: ch = ':'
            if '：' in line: ch = '：'
            info = line.split(ch)
            dictSNP[key][rsid].append(info[-1])
            head += 1
            continue
        if head >= 3:
            ch = ' '
            if ':' in line: ch = ':'
            if '：' in line: ch = '：'
            info = line.split(ch, 1)
            if len(dictSNP[key][rsid]) == 1:
                dictSNP[key][rsid].append({})
            dictSNP[key][rsid][1][info[0]] = info[1]
            head += 1
            continue
    return dictSNP

def printSnpInfo(dictSnpInfo):
    for key in dictSnpInfo:
        print(key)
        data = dictSnpInfo[key]
        for rsid in data:
            print('rsid', rsid)
            print('基因', data[rsid][0])
            if len(data[rsid]) > 1:
                dictInfo = data[rsid][1]
                for gene in dictInfo:
                    print(gene, dictInfo[gene])
            else:
                print('No Answer')
            print('-'*10)
        print()
    return

def saveDB(filename, arrayData):
    dbName = filename + '.db'
    conn = sqlite3.connect(dbName)
    cursor = conn.cursor()
    cursor.execute('drop table if exists user')
    cursor.execute('create table user(rsid varchar(32) primary key, chromosome varchar(2), position int, genetype varchar(2))')
    #cursor.execute('insert into user (rsid, chromosome, position, genetype) values (?,?,?,?)', arrayData[0])
    cursor.executemany('insert into user (rsid, chromosome, position, genetype) values (?,?,?,?)', arrayData)
    print('save Sqlite3 DB', cursor.rowcount)
    cursor.close()
    conn.commit()
    conn.close()
    return

def readDB(filename):
    dbName = filename + '.db'
    conn = sqlite3.connect(dbName)
    cursor = conn.cursor()
    cursor.execute('select * from user')
    listData = cursor.fetchall()
    print('load from Sqlite3', len(listData))
    cursor.close()
    conn.commit()
    conn.close()
    return np.array(listData)

def loadData(filename, forceReload=False):
    npyName = filename + '.npy'
    txtName = filename + '.txt'
    if os.path.isfile(npyName) and forceReload == False:
        arrayData = np.load(npyName)
        print('loaded from npy', npyName)
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
        print('loaded from txt', txtName)
    return arrayData

def searchSNP(arrayData, dictSNP):
    dictResult = {}
    for key in dictSNP:
        dictResult[key] = []
    for i in range(len(arrayData)):
        rsid = arrayData[i][0]
        for key in dictSNP:
            for SNP in dictSNP[key] :
                if rsid == SNP:
                    gene = arrayData[i][3]
                    if len(dictSNP[key][rsid]) == 1:
                        dictResult[key].append(list(arrayData[i])+['NoAnswer'])
                    elif gene not in dictSNP[key][rsid][1]:
                        dictResult[key].append(list(arrayData[i])+['基因异常']+[dictSNP[key][rsid][1]])
                    else:
                        dictResult[key].append(list(arrayData[i])+[dictSNP[key][rsid][1][gene]])
    return dictResult

def printResult(dictResult, dictSNP):
    for key in dictResult:
        print(key, '{}/{}'.format(len(dictResult[key]), len(dictSNP[key])))
        listrsid = []
        for result in dictResult[key]:
            print(result)
            listrsid.append(result[0])
        for info in dictSNP[key]:
            if info not in listrsid:
                print([info])
    print('-'*10)
    return

def saveResult(Name, filename, dictResult, dictSNP):
    output = '# {} \n'.format(Name)
    output += '{} \n'.format(filename)
    for key in dictResult:
        output += '#### ' + key + ' ({}/{}) \n'.format(len(dictResult[key]), len(dictSNP[key]))
        output += '|{}|{}|{}|{}|{}| \n'.format('rsid','chromo','position','type','结论')
        output += '|:---:|:---:|:---:|:---:|:---:| \n'
        listrsid = []
        for result in dictResult[key]:
            output += '|{}|{}|{}|{}|**{}**| \n'.format(result[0],result[1],result[2],result[3],result[4])
            listrsid.append(result[0])
        for info in dictSNP[key]:
            if info not in listrsid:
                output += '|{}|{}|{}|{}|**{}**| \n'.format(info,'-', '-', '-', '-')
    output += '---'
    f = open(Name+'.md','w',encoding='utf-8')
    f.write(output)
    f.close()
    return

def statChromosome(arrayData):
    dictStat = {}
    for i in range(len(arrayData)):
        data = arrayData[i]
        if data[1] not in dictStat:
            dictStat[data[1]] = 1
        else:
            dictStat[data[1]] += 1
    listKey = list(dictStat.keys())
    #listKey.sort()
    for key in listKey:
        print(key, dictStat[key])
    print('-'*10)
    return

def main():
    global Name, filename, snpname

    print('Check for {}'.format(Name))
    print(filename)

    dictSnpInfo = loadSnpInfo(snpname)
    #printSnpInfo(dictSnpInfo)
    print('SNPInfo Loaded')

    #arrayData = readDB(filename)
    #print('arrayData.shape', arrayData.shape)

    arrayData = loadData(filename)
    print('arrayData.shape', arrayData.shape)

    #saveDB(filename, arrayData)

    #statChromosome(arrayData)

    dictResult = searchSNP(arrayData, dictSnpInfo)
    printResult(dictResult, dictSnpInfo)

    saveResult(Name, filename, dictResult, dictSnpInfo)

    return

if __name__ == '__main__':
    main()

