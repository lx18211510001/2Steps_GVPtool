#2Steps_GVPtool

workflow of toolkit

2Steps_GVPtool is a toolkit to construct database containing Genetically Variant Peptides(GVP), such as ComVarDB databse, and identify GVPs from MS proteomics spectra data, using the constomized database. ComVarDB combines entries from SwissProt and GVPs from dbSNP.  
  


**Before running**  
	* Install ANNOVAR or other SNP annotation software  
	* Install Python  
  
  
**Detailed pipeline inputs**  
	* VCF files input  -- a tab deliminated text file
	* PSM files input  -- a tab deliminated text file

 
**Prepare once**  
	* Register for download of ANNOVAR  
	* Download reference protein sequence database of SwissProt  
	* Download SNP data from dbSNP database  
```
# Get this repo
git clone https://github.com/lx18211510001/2Steps_GVPtool
cd 2Steps_GVPtool

# Get Annovar
cd /path/to/your/annovar
wget __link_you_get_from_annovar__
tar xvfz annovar.latest.tar.gz
## This creates a folder with annotate_variation.pl and more files, to be passed to the pipeline with --annovar_dir

# Download all the canonical and isoform entried of human of SwissProt entry on the Uniprot website.

# Get SNP data from the dbSNP database
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_all_20180418.vcf.gz
```


**Analyse your files with 2Steps_GVPtool**
```
# Example command using test files

## SNP annotation using ANNOVAR
perl annovar_dir/convert2annovar.pl -format vcf4 -allsample -includeinfo -withfreq test/test.vcf > test/test.avinput &
perl annovar_dir/annotate_variation.pl -out test/test -build hg38 test/test.avinput hg38db/ &
perl annovar_dir/coding_change.pl test/test.exonic_variant_function hg38db/hg38_refGene.txt hg38db/hg38_refGeneMrna.fa --includesnp -onlyAltering -alltranscript -out test/test.fasta &

## Find GVPs
python GVP_generation.py -e test.exonic_variant_function -f test.fasta -o test.table

## GVP-containing database construction 
python db_format.py -db hgdb/SwissProt_9606.fasta -t test.table -g genelist -o SwissProt_9606_GVP.fasta

## GVP identification from the search result file
python GVP_pickResult.py -t SwissProt_9606_GVP_ref.table -db SwissProt_9606_GVP.fasta -i PSM.txt -pcol 13 -fcol 28
```
