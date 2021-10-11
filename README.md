# StV

HOR Structrual Variant (StV) prediction using CHM13 hor-monomer annotation

Usage: ./stv monomers.bed

Config: insert path to dir with program into stv.sh and scripts/mons2hors.py
___

Version 12

Changes:

• New AS-SF annotation

___

Files:

• t2t-chm13.v1.1-AS-SF-HORs-annotation.170921.bed - monomeric track

• t2t_cenAnnotation.v3.061021.CHM13v1.1.bed - annotation

• cenAnnotation_live.bed - coordinates of live HOR domain derived from file t2t_cenAnnotation.v3.070721.bed 
```shell 
grep 'L)' t2t_cenAnnotation.v3.061021.CHM13v1.1.bed | awk '{if ($3-$2>30000 && $3-$2!=36169) {print $0} }' > cenAnnotation_live.bed
```

• stv12.bed - resulting file which can be put in the browser

• stv_row.bed - same as stv12.bed but doesn't contain stv numbering, colors, the first description line

• stats.tsv - contains the number of each stv in each chr
___

Metfod description (from Supplementary Material of https://www.biorxiv.org/content/10.1101/2021.07.12.452052v1):

HOR Structrual Variant (StV) prediction using CHM13 hor-monomer annotation 

The StV track was derived from the AS HOR annotation track (HOR-track; AS_HOR_Annot.bed) by using python scripts. The HOR-track shows the coordinates of every monomer in the assembly and its name which indicates to which HOR the monomer belongs and what is the number of this monomer in a HOR (defined as a standard master HOR with a fixed cyclic shift, see the HOR-track description and Uralsky 2019). So, each array is identified by a HOR name (e.g. S3CXH1L for the live array of the X chromosome) and each monomer by its number (e.g. S3CXH1L.1 for monomer #1 of the HOR). Hybrid monomers are indicated with a “/” (e.g. S2C8H1L.4/7, where the first part of the monomer is S2C8H1L.4 and the second is S2C8H1L.4/7).

The main idea of generating the StV map was iterating through AS HOR annotation bed file and cutting after the last monomer in a HOR (the one with maximal number in a given array) or before the first monomer (the one with number #1) if the last monomer is not present in a given StV. For arrays where AS was on reverse strand, a cut was performed after the first monomer or before the last monomer if the first monomer is not present in a given StV. Note that about a half of the live centromeres have AS on direct strand and another half on reverse strand. Some arrays may include an inversion (e.g. chr1:124130688-125857749).

StV naming. The first part of the StV name is the HOR name (e.g. S2C8H1L). Then, after “.” follows listing of monomer numbers included in StV separated by “_”. If monomer numbers go in natural order (or descending natural order for the reverse strand HORs), we join the first and the last monomers in this string with “-” to indicate an interval and do not show the numbers in between (e.g. S3C11H1L.1_2_3 would appear S3C11H1L.1-3). If AS is on reverse strand, the monomer numbers go in reverse (e.g. S3CXH1L.1-12 would appear as S3CXH1L.12-1). Hybrids are always flanked by “_” (e.g. S2C8H1L.1-3_4/7_8-11). Rarely occurring monomers identified as monomers of another HOR or as SF class monomers (usually due to misclassification) are shown by their own name (e.g. S2C15H1L.1-4_S5C1qH6d.2_5-11 or S3C11H1L.1-2_W3_5).

The coordinates of the live HORs were taken from the centromere annotation track (t2t_cenAnnotation.v2.021921.bed). The resulting full-length HORs and StVs were counted and listed in the stats table. Reverse strand HOR and StV names were reversed in the stats table to make it easier to read (e.g. S3CXH1L.12-11_8-1 would become S3CXH1L.1-8_11-12). In centromeres 1, 5 and 19, there are long chains of "_6/4_5” repeats which appear within a HOR (HOR formula S1C1/5/19H1L.1-5(_6/4_5){n}-6). In StV names these repeats were shown as “(_6/4_5){n}” (e.g. S1C1/5/19H1L.1-5(_6/4_5){7}-6 where {7} indicates the number of repeats in a given HOR). In the stats table all numbers in “{}” were replaced by “{n}” because we consider them as one StV, however the statistics for {n} is collected by the script and can be retrieved. The StVs and full-length HORs in every centromere were numbered. A number is at the beginning of StV’s name after “#” sign. After a number there is “:” sign (e.g. #117:S3CXH1L.12-11_8-1). Numbering goes from p- to q-side of a chromosome regardless of the AS direction. The StV track was colored randomly but rare StVs (<10% in a centromere) are in bright colors to increase visibility and common StVs and full-length HORs are in dim colors.
