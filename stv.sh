#!/bin/bash
# Pipeline creates StV maps and StV stats table from monomeric bed file

if [[ $# -ne 1 ]]
then
    echo "usage: ./stv <monomeric track>"
    exit 1
fi

# put here path to dir with program
path='/home/fedor/Programs/my/stv'

# monomers to StVs
python3 $path/scripts/mons2hors.py $1 > stv.bed

# manualy replace wrong StV around cen1 inversion
left='chr1\t124128489\t124130518\tS1C1/5/19H1L.1-5(_6/4_5){3}_S1C16H1L.4\t0\t+\t124128489\t124130518\t0,0,0\nchr1\t124130518\t124131251\tS1C1/5/19H1L.5-1\t0\t-\t124130518\t124131251\t0,0,0'
right='chr1\t125857239\t125858698\tS1C1/5/19H1L.6-(5_6/4_){3}5-4_6/4\t0\t-\t125857239\t125858698\t0,0,0\nchr1\t125858698\t125859450\tS1C1/5/19H1L.6/4_5(_6/4_5){1}-6\t0\t+\t125858698\t125859450\t0,0,0'
awk -v corrected_stv=$left '{if ($1=="chr1" && $2==124128489 && $3==124131251) {print corrected_stv} else {print $0}}' stv.bed > stv_corrected.bed
awk -v corrected_stv=$right '{if ($1=="chr1" && $2==125857239 && $3==125859450) {print corrected_stv} else {print $0}}' stv_corrected.bed > stv_row.bed

# stats file
python3 $path/scripts/bed2stat.py stv_row.bed > stats.tsv

# coloring
python3 $path/scripts/coloring.py stv_row.bed > stv_colored.bed

# numbering
python3 $path/scripts/numbering.py stv_colored.bed > stv.bed

# add descriotion line
sed -i "1 i\track name=\"StV\" description=\"Structural Variants\" itemRgb=\"On\" visibility=\"1\"" stv.bed

rm stv_colored.bed
rm stv_corrected.bed
