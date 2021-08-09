# -*- coding: utf-8 -*-


import os
import os.path
import argparse

parser = argparse.ArgumentParser(description = '\n该脚本用于提取搜库结果PSM文件中的GVP', add_help = False, usage = '\npython3 GVP_pickResult.py -t [[GVP_ref.table] -db [reference_databse_GVP.fasta] -i [result_PSM.txt] -pcol [peptide_column] -fcol [fdr_column]\npython3 GVP_pickResult.py -table [[GVP_ref.table] -GVPdb [reference_databse_GVP.fasta] -input [result_PSM.txt] -pepColumn [peptide_column] -fdrColumn [fdr_column]')
required = parser.add_argument_group('必选项')
optional = parser.add_argument_group('可选项')
required.add_argument('-t', '--table', metavar = '[GVP_ref.table]', help = 'db_format.py生成的包含GVP信息的表格文件，以tab分隔', required = True)
required.add_argument('-db', '--GVPdb', metavar = '[reference_databse_GVP.fasta]', help = 'db_format.py生成的包含GVP的参考蛋白质序列数据库，fasta格式', required = True)
required.add_argument('-i', '--input', metavar = '[result_PSM.txt]', help = '搜库软件输出的PSM结果文件，tab分隔', required = True)
required.add_argument('-pcol', '--pepColumn', metavar = '[peptide_column]', help = 'PSM文件中肽段所在的列', required = True)
required.add_argument('-fcol', '--fdrColumn', metavar = '[fdr_column]', help = 'PSM文件中肽段所在的列', required = True)
optional.add_argument('-f', '--fdr', metavar = '[fdr_cutoff]', help = '筛选肽段的fdr值，默认为0.01', required = False)
optional.add_argument('-h', '--help', action = 'help', help = '帮助信息')
args = parser.parse_args()


def read_GVPinfo(GVPrefFile, fastaFile):
    GVPinfodic = {}
    with open(GVPrefFile, 'r') as file:
        for line in file:
            GVP = line.split('\t')[-2]
            GVPinfodic[GVP] = line.replace('\n','').split('\t')
    
    GVPseq = ''
    rawseq = ''
    GVPflag = 'T'
    with open(fastaFile, 'r') as file:
        for line in file:
            if line.startswith('>'):
                if '.m' not in line:
                    GVPflag = 'F'
            else:
                if GVPflag == 'T':
                    GVPseq = GVPseq + line.replace('\n',' ').replace('I', 'L')
                else:
                    rawseq = rawseq + line.replace('\n',' ').replace('I', 'L')
    return GVPinfodic, GVPseq, rawseq



def pick_GVP_result(PSMfile, pepCol, fdrCol, fdrCutoff):
    global GVPseq, GVPinfodic
    
    GVPresultdic = {}
    GVPfdrdic = {}
    with open(PSMfile, 'r') as file:
        for line in file:
            fdr = line.replace('\n','').split('\t')[fdrCol-1]
            if str.isdigit(fdr.replace('.','')):
                if float(fdr) < fdrCutoff:
                    pep = line.replace('\n','').split('\t')[pepCol-1].replace('I','L')
                    if pep in GVPseq:
                        for (GVP, info) in GVPinfodic.items():
                            if pep in GVP.replace('I','L'):
                                if pep in GVPresultdic.keys():
                                    GVPresultdic[GVP] = GVPresultdic[GVP] + ';' + pep
                                    GVPfdrdic[GVP] = GVPfdrdic[GVP] + ';' + fdr
                                else:
                                    GVPresultdic[GVP] = pep
                                    GVPfdrdic[GVP] = fdr
    return GVPresultdic, GVPfdrdic


def save2file(PSMfile):
    global GVPresultdic, GVPinfodic, rawseq
    with open(os.path.splitext(PSMfile)[0] + '.GVP', 'w') as file:
        head = 'GVP\tgene\trsID\tchr\tlocation\thgvs.c\thgvs.p\twildtype peptide\tdetected peptide\tfdr\tunique\n'
        file.write(head)
        for (GVP, pepIL) in GVPresultdic.items():
            line = GVP + '\t' + '\t'.join(GVPinfodic[GVP][:6]) + '\t' + GVPinfodic[GVP][-1] + '\t' + pepIL + '\t' + GVPfdrdic[GVP] + '\t'
            
            if pepIL in rawseq:
                line = line + 'N'
            else:
                line = line + 'Y'
            file.write(line + '\n')
        

GVPinfoFile = args.table
fastaFile = args.GVPdb
PSMfile = args.input
pepCol = eval(args.pepColumn)
fdrCol = eval(args.fdrColumn)
if args.fdr:
    fdrCutoff = eval(args.fdr)
else:
    fdrCutoff = 0.01

GVPinfodic, GVPseq, rawseq = read_GVPinfo(GVPinfoFile, fastaFile)
GVPresultdic, GVPfdrdic = pick_GVP_result(PSMfile, pepCol, fdrCol, fdrCutoff)
save2file(PSMfile)

