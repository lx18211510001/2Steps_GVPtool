# -*- coding: utf-8 -*-

import os
import os.path
import argparse

parser = argparse.ArgumentParser(description = '\n该脚本用于生成包含GVP的新的数据库', add_help = False, usage = '\npython3 db_format.py -db [reference_db.fasta] -t [GVP.table] -genelist [genelist] -o [reference_db_GVP.fasta]\npython3 db_format.py -database [reference_db.fasta] -table [GVP.table] -genelist [genelist] -output [reference_db_GVP.fasta]')
required = parser.add_argument_group('必选项')
optional = parser.add_argument_group('可选项')
required.add_argument('-db', '--database', metavar = '[reference_db.fasta]', help = '参考蛋白质序列数据库，fasta格式', required = True)
required.add_argument('-t', '--table', metavar = '[GVP.table]', help = 'GVP_generation.py生成的包含GVP信息的表格文件，以tab分隔', required = True)
required.add_argument('-o', '--output', metavar = '[reference_db_GVP.fasta]', help = '在原数据库中加入了GVP序列的新的数据库，fasta格式', required = True)
optional.add_argument('-g', '--genelist', metavar = '[genelist]', help = '若该参数存在，则在新建的数据库中仅包含这部分基因的GVP序列', required = False)
optional.add_argument('-h', '--help', action = 'help', help = '帮助信息')
args = parser.parse_args()


def read_refDB(refDBfile):
    refDBseq = ''
    with open(refDBfile, 'r') as file:
        for line in file:
            if line.startswith('>'):
                refDBseq = refDBseq + ' '
            else:
                refDBseq = refDBseq + line.replace('\n', '').replace('I', 'L')
    return refDBseq


def read_pepInfo(fn):
    genedic = {}
    gene_raw = ''
    pep_raw = ''
    with open(fn, 'r') as file:
        for line in file:
            gene = line.split(',')[0]
            pep = line.split(',')[1]
            if gene != '':
                genedic[gene] = {}
                genedic[gene][pep] = [line.replace('\n','').split(',')[2:]]
                gene_raw = gene
                pep_raw = pep
            elif gene == '' and pep != '':
                genedic[gene_raw][pep] = [line.replace('\n','').split(',')[2:]]
                pep_raw = pep
            else:
                genedic[gene_raw][pep_raw].append(line.replace('\n','').split(',')[2:])
    return genedic


def read_tablefile(tableFile):
    global refDBseq
    geneRsdic = {}
    gene_raw = ''
    with open(tableFile, 'r') as file:
        for line in file:
            gene = line.split('\t')[0]
            GVP = line.split('\t')[6]
            if GVP.replace('I','L') not in refDBseq:
                rsID = line.split('\t')[1]
                if gene != '':
                    geneRsdic[gene] = {}
                    gene_raw = gene
                geneRsdic[gene_raw][rsID] = line.split('\t')[1:]
    return geneRsdic
    

def genels_filter(geneRsdic, genelistFile):
    genels = []
    with open(genelistFile, 'r') as file:
        for line in file:
            genels.append(line.replace('\n',''))
            
    filter_geneRsdic = {}
    for (gene, rsdic) in geneRsdic.items():
        if gene in genels:
            filter_geneRsdic[gene] = rsdic
    return filter_geneRsdic


def db_generation(newDBfile):
    global refDBfile, geneRsdic
    with open(newDBfile, 'w') as file:
        
        for (gene, rsdic) in geneRsdic.items():
            pepls = []
            for (rsID, ls) in rsdic.items():
                pepls.append(ls[5])
            
            endPepls = []
            for i in range(len(pepls)):
                if pepls[i][-1] != 'K' and pepls[i][-1] != 'R':
                    endPepls.append(pepls[i])                    
            
            pepStopNumls = []
            if len(endPepls) == 0:
                newPepls = pepls
                pepStopNumls.append(len(pepls)-1)
            else:
                newPepls = []
                for i in range(len(pepls)):
                    if pepls[i] != endPepls[-1]:
                        newPepls.append(pepls[i])
                newPepls.append(endPepls[-1])
                
                for i in range(len(newPepls)):
                    if newPepls[i][-1] != 'K' and newPepls[i][-1] != 'R':
                        pepStopNumls.append(i)
            
            seq = ''
            m = 0
            start = 0
            for i in pepStopNumls:
                stop = i
                for j in range(start,stop+1):
                    seq = seq + newPepls[j]
                accesion = gene + '.m' + str(m)
                head = '>sp|' + accesion + '|' + accesion + '_HUMAN OS=Homo sapiens OX=9606 GN=' + accesion
                file.write(head + '\n' + seq + '\n')
                seq = ''
                m += 1
                start = stop+1
        
        with open(refDBfile, 'r') as reffn:
            file.write(reffn.read())


def saveGVPinfo2table(newDBfile):
    global geneRsdic
    head = '#gene\trsID\tchr\tlocation\thgvs.c\thgvs.p\tGVP\twildType peptide\n'
    with open(os.path.splitext(newDBfile)[0] + '_ref.table', 'w') as file:
        file.write(head)
        for (gene, rsdic) in geneRsdic.items():
            for (rsID, ls) in rsdic.items():
                line = gene + '\t' + '\t'.join(ls[:7]) + '\n'
                file.write(line)
        

refDBfile = args.database
tableFile = args.table
newDBfile = args.output

refDBseq = read_refDB(refDBfile)
geneRsdic = read_tablefile(tableFile)
if args.genelist:
    geneRsdic = genels_filter(geneRsdic, args.genelist)
db_generation(newDBfile)
saveGVPinfo2table(newDBfile)

