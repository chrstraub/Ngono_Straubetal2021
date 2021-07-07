# N. gonorrhoeae survey - analysis 
Please see Material & Methods for software version numbers and links to download the programs. 

Sequencing reads for this study have been deposited in the National Centre for Biotechnology Information (NCBI) database under BioProject number PRJNA738299.
For some analysis, additional raw sequencing reads from  an earlier NZ N. gonorrhoeae survey in 2014/2015 were used (Lee et al. 2018, https://doi.org/10.1093/jac/dkx405), available from BioProject number PRJNA394216.

## **In silico typing**
### MLST
``` bash
#!/bin/bash
mlst --scheme=neisseria --novel=novel.fa Ng2018_survey/assemblies/*
```

### NG-MAST
```bash
#!/bin/bash
#update DB
ngmaster --updatedb --db NGMAST/db

#run ngmaster on assemblies
ngmaster --db NGMAST/db --csv --printseq NGMAST_alleles_extracted.fasta  assemblies/*.fa > NGMAST.csv
```

### NG-STAR
```bash
#!/bin/bash
python3 pyngSTar.py -f -a -i assemblies/*.fa -p pyngSTarDB/ - ngstar.csv
```

## **Antimicrobial resistance**
### Plasmids - Abricate
``` bash
#!/bin/bash
module load abricate/0.8.7
dir=assemblies
out=ncbi_out

for file in ${dir}/*.fa; do
    echo $file
    id=$(basename $file .fa)
    echo $id
    abricate --db ncbi --quiet ${file} > ${out}/${id}.tab
done

#summary
for file in ${out}/*.tab; do
    var="$var $file"
done
echo $var

abricate --summary ${var} > ${out}/summary_ncbi.tab
```
### Plasmid types - isPCR
```bash
#!/bin/bash
cat assemblies/*.fa > db.fasta
/home/software/exonerate-2.2.0-x86_64/bin/ipcress --products --mismatch 1 primers_plasmid.txt db.fasta
/home/software/exonerate-2.2.0-x86_64/bin/ipcress --products --mismatch 1 primers_blaTEM.txt db.fasta

#a parser for ipcress results can be found here: https://github.com/chrstraub/iPCRess_parser
```

### Chromosomal mutations - Ariba
see `Ariba_gono_database/` for creating the custom database `db_custom_gono`

```bash
#create inputfile reads.tab with 2 columns (R1.fastq.gz  R2.fastq.gz) for all isolates
while read -r line
do
    R1=$(echo $line | awk '{print $1}')
    R2=$(echo $line | awk '{print $2}')
    id=$(( basename $R1 .fastq.gz)| awk -F '_' '{print $1}')
    ariba run db_custom_gono $R1 $R2 output/${id}_out
done < reads.tab
``

### 23S mutations
```bash
#!/bin/bash
module load snippy/4.3.6
module load samtools/1.9

#1 Run Snippy on all isolates against a reference where three alleles of 23S are masked
while IFS=$'\t' read -r id R1 R2; do
    snippy --cpus 20 --outdir $id --ref NG_NCCP11945_23S-masked.gbk --R1 $R1 --R2 $R2
done < nullarbor_input_final.tab

#2 Look at SNPs in 23S and count number of occurences
source activate py2_base #python 2 conda environment
while read -r line; do
    ids=`awk 'BEGIN {FS="\t"}; {print $1}'`
done < nullarbor_input_final.tab

echo ${ids}

python ng23S-mutations.py ${ids} > results_final.txt
cut -f1 -d$'\t' results_final.txt | sort | uniq -c | wc -l
```


## SNP calling & filtering
### **SNP calling workflow:**   
A snakemake workflow for mapping reads to a reference and SNP calling. 
see `https://github.com/chrstraub/SNP_calling_workflow`
1) index reference N. gonorrhoeae reference FA1090
2) map with bwa to reference
3) convert sam to bam
4) sort bam
6) mark duplicates (PICARD)
7) call SNPs with Freebayes (as a population, not individual)

**Raw SNP calls vcf file**
`pop_call_raw.vcf`

### Filtering steps
see `SNP_filter_popvcf.sh`
using conda environments for the various tools

#### Further info regarding masking repetitive regions (creation of bed files)
1. prophages - see `phage.bed`
2. repeat regions (MUMMER) - see `repeat_n20.bed`
3. recombinogenic regions (gubbins)

### **Prophages**
from doi:10.1186/1471-2180-7-66 Table 2
*ds DNA prophage sequences*
|Prophage|Sequence coordinates|Length of the DNA sequence (bp)|CDS annotations (Acc. No AE004969)|Number of CDS|
|---|---|---|---|---|
NgoΦ1|455173 – 498100|42 927|NGO0462-NGO0524|63
NgoΦ2|1044447 – 1078281|33 834|NGO1085-NGO1132|48
NgoΦ3|1583028 – 1599049|16 021|NGO1613-NGO1640|284
| |1606828 – 1609808|2 980|NGO1649-NGO1652|
NgoΦ4|972837 – 984008|11 171|NGO1000-NGO1020|21
NgoΦ5|721256 – 729870|8 614|NGO0720-NGO0732|13

*Filamentous ss DNA prophage sequences*
|Prophage|Sequence coordinates|Length of the DNA sequence (bp)|CDS annotations (Acc. No AE004969)|Number of CDS|
|---|---|---|---|---|
NgoΦ6|1080185 – 1088420|8 235|NGO1137 – NGO1146|13
NgoΦ7|1215967 – 1223383|7 416|NGO1262 – NGO1270|12
NgoΦ8|1103017 – 1109444|6 427|NGO1164 – NGO1170|8
NgoΦ9|1599378 – 1607537*|7 159|NGO1641 – NGO1648|9

### **Repeat regions - MUMMER**
```bash
source activate bcftools
conda install mummer
#installed version 3.23
ref=Ngono_FA1090.fna

#To find exact repeats of length 50 or greater in a single sequence seq.fasta, type:
repeat-match -t -n 20 -f ${ref} > temp.out   

echo "Tandem repeats"
sort -k1n -k2n temp.out | awk -f /home/anaconda3/envs/bcftools/opt/mummer-3.23/scripts/tandem-repeat.awk > tandem_repeats20.out
rm temp.out
```

```bash
#make bed file
#column 1 = chromosome
#column 2 = star
#column 3 = end

sed 1d tandem_repeats20.out | while IFS=$' ' read -r start extent len copy; 

do
    ref=NC_002946.2
    echo $start
    end=$(expr ${start} + ${extent})
    echo $end
    echo -e "${ref}\t${start}\t${end}" >> repeat_n20.bed
done
```

## Convert filtered SNP vcf to MSA for input for gubbins (output after running SNP_filter_popvcf.sh)
```bash
bgzip -c 2018_final_filt.vcf > 2018_final_filt.vcf.gz
tabix 2018_final_filt.vcf.gz 

#create fasta sequence from reference with SNPs substituted for a single sample
bcftools consensus -s 19AR0410z -f Ngono_FA1090.fna 2018_final_filt.vcf.gz > 19AR0410z.fa

source activate bcftools 
sample=$(bcftools query -l 2018_final_filt.vcf.gz)

for id in $sample;do 
id2=${id%z}
echo $id2
bcftools consensus -s $id -f Ngono_FA1090.fna 2018_final_filt.vcf.gz > SNP_fasta/${id2}_SNP.fasta
done

#Rename fasta headers
for file in *.fasta
do
        name=$(basename "$file" _SNP.fasta)
        awk -v name="$name" '/^>/{print ">" name; next}{print}' $file > rename_$name.fasta
done;

cat rename_*.fasta > 2018_all_SNP.fasta
cat 2018_all_SNP.fasta Ngono_FA1090.fna >> 2018_all_SNP.aln
```

## Gubbins

```bash
#!/bin/bash
in_dir=gubbins/in

source activate gubbins &&
run_gubbins.py ${in_dir}/2018_all_SNP.aln --tree_builder raxml -i 30 -u -v -p gubbins_all

```
## Phylogenetic tree
Gubbins output:  
- gubbins_all.filtered_polymorphic_sites.fasta --> 13,003 nucleotides
- gubbins_1815_40it.filtered_polymorphic_sites.fasta --> 14,602 nucleotides

use .filtered_polymorphic_sites.fasta (contains Ns), remove Ns and build a new tree with IQtree

```bash
#filter out sites other than ACGT
snp-sites -c gubbins_all.filtered_polymorphic_sites.fasta > gubbins_clean_variant.fasta #5,059 ACGT sites, 3,871 informative sites
snp-sites -c gubbins_1815_40it.filtered_polymorphic_sites.fasta > gubbins_1815_40it_clean_variant.fasta #14,602 ACGT sites,  informative sites

iqtree -s gubbins_clean_variant.fasta -B 1000 -T 4
```

Analysis for hierarchical clustering (hierBAPS), as well as all analysis based on pairwise SNP-distances was done in R - please see separate `R_analysis folder`.