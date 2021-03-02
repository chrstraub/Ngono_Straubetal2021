# Ariba - create custom DBs for Neisseria gonorrhoeae
### **For variants**
Input files:

- alignment files in: `prepare_ref/aln2meta_input/*.aln`
- metadata filesin: `prepare-ref/aln2meta_input/*.tsv`

Alignment file contains all different alleles (nucleotide sequences) - choose simple fasta headers. Must be alignment - all same length (including gaps).
For each gene, create a unique alignment and metadata file.

The metadata (.tsv) file contains the SNP information for each gene, e.g.: 
```bash
gyrA.WHO_F_00668c       S91F    gyrA.91F        Resistance to fluoroquinolones
gyrA.WHO_F_00668c       D95N    gyrA.95N        Resistance to fluoroquinolones
gyrA.WHO_F_00668c       D95G    gyrA.95G        Resistance to fluoroquinolones
gyrA.WHO_F_00668c       D95Y    gyrA.95Y        Resistance to fluoroquinolones
gyrA.WHO_F_00668c       D95A    gyrA.95A        Resistance to fluoroquinolones
```

Run **aln2meta**   and **prepareref**
You can specify coding/non-coding in your command. If you have a batch of coding and non-coding, you have to run the command separately. Prefix for out directory. 
ARIBA will sanity check the specified SNPs against the sequences and output warnings if necessary.

```bash
#non coding
for x in 16S 23S
do
    ariba aln2meta --variant_only prepare_ref/aln2meta_input/${x}.aln \
    prepare_ref/aln2meta_input/${x}_in.tsv noncoding prepare_ref/aln2meta_output/${x}
done

#coding
for x in folP gyrA mtrR parC parE penA pilQ ponA porB1b rplD rpoB rpsE rpsJ
do
    echo ${x} &&
    ariba aln2meta --variant_only prepare_ref/aln2meta_input/${x}.aln \
    prepare_ref/aln2meta_input/${x}_in.tsv coding prepare_ref/aln2meta_output/${x} &&
    echo ${x} finished
done

#and run **prepareref**:  
#Before using prepareref, the files need to be concatenated.
cat prepare_ref/aln2meta_output/out*fa > prepare_ref/aln2meta_output/all_variant.fa
cat prepare_ref/aln2meta_output/out*tsv > prepare_ref/aln2meta_output/all_variant.tsv
cat prepare_ref/aln2meta_output/out*cluster > prepare_ref/aln2meta_output/all_variant.cluster

ariba prepareref -f prepare_ref/aln2meta_output/all_variant.fa -m prepare_ref/aln2meta_output/all_variant.tsv --cdhit_clusters prepare_ref/aln2meta_output/all_variant.cluster db_variant
```

### **For presence/absence**
Input files: 
- `prepare_ref/pres_abs/pres_abs_all.fa`
- `prepare_ref/pres_abs/pres_abs_meta.tsv`

FASTA files contains nucleotide sequences (e.g. promoter seqs or genes) for rplV, porA, mtrR promoter, macAB promoter and norM promoter:

```bash
rplV.NEIS0137_allele1   1       0       .       .       .
rplV.NEIS0137_allele10  1       0       .       .       .
porA.WHO_U_01087c       1       0       .       .       .
porA.NCCP11945  1       0       .       .       .
mtrR_promoter.WHO_F     0       0       .       .       .
macAB_promoter.strain35_02      0       0       .       .       .
macAB_promoter.FA10190  0       0       .       .       .
norM_promoter.FA10190   0       0       .       .       .

```

Run **prepareref** to prepare the presence/absence database:  
```bash
ariba prepareref -f prepare_ref/pres_abs/pres_abs_all.fa -m prepare_ref/pres_abs/pres_abs_meta db_pres_abs
```

Use the two folders containing the custom gono DBs:  
`db_variant/`   
`db_pres_abs/`   
for ARIBA analysis!