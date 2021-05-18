# Usage: python3 bed2stat.py <input>.bed > <stat>.tsv
from collections import Counter
from re import sub, search
from sys import argv


def count_uniq_stv(stv_bed):
    all_stv = [sub(r'{\d+}', '{n}', line[3]) for line in stv_bed]
    return dict(Counter(all_stv).most_common())


# parse input bed
input_bed = []
with open(argv[1]) as bed:
    for line in bed:
        line = line.split()
        if line[0] != 'track':
            input_bed.append(line)

input_bed.append(['chr24'] + input_bed[0][1:])  # idk wtf is with last chr. Let's add shit for chrX not to be the last chr. 
current_chr = input_bed[0][0]
current_bed = []
start = input_bed[0][1]
end = 0
for line in input_bed:
    chr = line[0]
    if chr == current_chr:
        current_bed.append(line)
        end = line[2]
    else:
        stv_stat = count_uniq_stv(current_bed)
        print('{}:{}-{}'.format(current_chr, start, end))
        print('total', sum(stv_stat.values()), sep='\t')
        for i in stv_stat:
            print(i, stv_stat[i], sep='\t')
        current_chr = chr
        current_bed = []
        start = line[1]
