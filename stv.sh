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
left='chr1\t124128658\t124130687\tS1C1/5/19H1L.1-5(_6/4_5){3}_S1C16H1L.4\t0\t+\t124128658\t124130687\t0,0,0\nchr1\t124130687\t124131420\tS1C1/5/19H1L.5-1\t0\t-\t124130687\t124131420\t0,0,0'
right='chr1\t125857408\t125858868\tS1C1/5/19H1L.6-(5_6/4_){3}5-4_6/4\t0\t-\t125857408\t125858868\t0,0,0\nchr1\t125858868\t125859619\tS1C1/5/19H1L.6/4_5(_6/4_5){1}-6\t0\t+\t125858868\t125859619\t0,0,0'
awk -v corrected_stv=$left '{if ($1=="chr1" && $2==124128658 && $3==124131420) {print corrected_stv} else {print $0}}' stv.bed > stv_corrected.bed
awk -v corrected_stv=$right '{if ($1=="chr1" && $2==125857408 && $3==125859619) {print corrected_stv} else {print $0}}' stv_corrected.bed > stv_row.bed

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
