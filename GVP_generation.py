# -*- coding: utf-8 -*-


import re
import argparse

parser = argparse.ArgumentParser(description = '\n该脚本用于根据SNP的注释信息提取出GVP信息', add_help = False, usage = '\npython3 GVP_generation.py -e [SNP.exonic_variant_function] -f [SNP.fasta] -o [SNP_GVP.table]\npython3 GVP_generation.py --exonic [SNP.exonic_variant_function] --fasta [SNP.fasta] --output [SNP_GVP.table]')
required = parser.add_argument_group('必选项')
optional = parser.add_argument_group('可选项')
required.add_argument('-e', '--exonic', metavar = '[SNP.exonic_variant_function]', help = 'ANNOVAR注释的外显子区域的SNP信息的文件，以tab作为分隔', required = True)
required.add_argument('-f', '--fasta', metavar = '[input.fasta]', help = 'ANNOVAR注释的外显子区域的SNP突变前后对应的蛋白质序列文件，fasta格式', required = True)
required.add_argument('-o', '--output', metavar = '[output.fasta]', help = '输出文件，以tab作为分隔的表格格式', required = True)
optional.add_argument('-h', '--help', action = 'help', help = '帮助信息')
args = parser.parse_args()


def readExonicInfo(filename):
    lineExondic = {}
    with open(filename, 'r') as file:
        for line in file:
            if line.split('\t')[1] != 'synonymous SNV' and line.split('\t')[1] != 'unknown':
                for muInfo in line.split('\t')[2].split(',')[:-1]:
                    key = line.split('\t')[0] + '-' + muInfo.split(':')[1]
                    lineExondic[key] = {}
                    lineExondic[key]['rsID'] = line.split('\t')[13]
                    lineExondic[key]['chr'] = line.split('\t')[3]
                    lineExondic[key]['Location'] = line.split('\t')[4]
                    lineExondic[key]['gene'] = muInfo.split(':')[0]
                    lineExondic[key]['hgvs.c'] = muInfo.split(':')[3]
                    lineExondic[key]['hgvs.p'] = muInfo.split(':')[-1]
                    lineExondic[key]['mutation_info'] = line.split('\t')[2]
                    lineExondic[key]['rs_info'] = line.replace('\n','').split('\t')[-1]
    return lineExondic


def KR_digestion(seq):
    pepls = []
    KR = sorted([i.start() for i in re.finditer('K', seq)] + [i.start() for i in re.finditer('R', seq)])
    pre = 0
    for i in KR:
        pep = seq[pre:i+1]
        if pep not in pepls:
            pepls.append(pep)
        pre = i+1
    pep = seq[pre:]
    if pep not in pepls:
        pepls.append(pep)
    return pepls


def get_dif_peptide(rawSeq, snpSeq):
    rawpepls = KR_digestion(rawSeq)
    snppepls = KR_digestion(snpSeq)
    
    same = [x for x in rawpepls if x in snppepls]
    dif_rawpepls = [y for y in rawpepls if y not in same]
    dif_snppepls = [z for z in snppepls if z not in same]
    return dif_rawpepls, dif_snppepls


def readFastafile(filename):
    lineGVPdic = {}
    with open(filename, 'r') as file:
        tempdic = {}
        lineNumNM = '0'
        for line in file:
            if line.startswith('>') and 'WILDTYPE' in line:
                if len(tempdic) == 2:
                    dif_rawpepls, dif_snppepls = get_dif_peptide(tempdic['rawSeq'], tempdic['snpSeq'])
                    lineGVPdic[lineNumNM] = {}
                    lineGVPdic[lineNumNM]['raw_pepls'] = dif_rawpepls
                    lineGVPdic[lineNumNM]['snp_pepls'] = dif_snppepls
                    tempdic = {}

                key = 'rawSeq'
                lineNumNM = line.split(' ')[0][1:] + '-' + line.split(' ')[1]
                tempdic[key] = ''
            elif line.startswith('>'):
                key = 'snpSeq'
                tempdic[key] = ''
            else:
                tempdic[key] = tempdic[key] + line.replace('\n','').replace('*','')
            
        dif_rawpepls, dif_snppepls = get_dif_peptide(tempdic['rawSeq'], tempdic['snpSeq'])
        lineGVPdic[lineNumNM] = {}
        lineGVPdic[lineNumNM]['raw_pepls'] = dif_rawpepls
        lineGVPdic[lineNumNM]['snp_pepls'] = dif_snppepls
            
        return lineGVPdic
                

def mergeInfo2GeneLevel(pepLengthCutoff):
    global lineExondic, lineGVPdic
    
    geneRsGVPdic = {}
    for (key, dic) in lineGVPdic.items():
        if len(dic['snp_pepls']) == 1 and len(dic['snp_pepls'][0]) >= pepLengthCutoff[0] and len(dic['snp_pepls'][0]) <= pepLengthCutoff[1]:
            gene = lineExondic[key]['gene']
            rsID = lineExondic[key]['rsID']
            if gene in geneRsGVPdic.keys():
                if rsID not in geneRsGVPdic[gene].keys():
                    geneRsGVPdic[gene][rsID] = lineExondic[key]
                    geneRsGVPdic[gene][rsID]['GVP'] = dic['snp_pepls'][0]
                    geneRsGVPdic[gene][rsID]['WTpep'] = '/'.join(dic['raw_pepls'])
            else:
                geneRsGVPdic[gene] = {}
                geneRsGVPdic[gene][rsID] = lineExondic[key]
                geneRsGVPdic[gene][rsID]['GVP'] = dic['snp_pepls'][0]
                geneRsGVPdic[gene][rsID]['WTpep'] = '/'.join(dic['raw_pepls'])

    return geneRsGVPdic


def save2table(tablefile):
    global geneRsGVPdic
    
    head = '#gene\trsID\tchr\tlocation\thgvs.c\thgvc.p\tGVP\tWildtype peptide\tannotation_exon\tannotation_VCF\n'
    with open(tablefile, 'w') as file:
        file.write(head)
        for (gene, rsGVPdic) in geneRsGVPdic.items():
            file.write(gene)
            for (rs, GVPdic) in rsGVPdic.items():
                ls = [rs, GVPdic['chr'], GVPdic['Location'],GVPdic['hgvs.c'],GVPdic['hgvs.p'],
                      GVPdic['GVP'],GVPdic['WTpep'],GVPdic['mutation_info'],GVPdic['rs_info']]
                file.write('\t' + '\t'.join(ls) + '\n')


pepLengthCutoff = [4,50]
exonInputfile = args.exonic
fastaInputfile = args.fasta
tablefile = args.output

lineExondic = readExonicInfo(exonInputfile)
lineGVPdic = readFastafile(fastaInputfile)
geneRsGVPdic = mergeInfo2GeneLevel(pepLengthCutoff)
save2table(tablefile)

