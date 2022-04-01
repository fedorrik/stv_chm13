# Convert monomeric bed file into the HOR bed file
# Usage: python3 mons2hors.py <mons>.bed > <hors>.bed
# Author: Fedor Ryabov 
from collections import Counter
from random import choice, randint
from re import search, sub
from sys import argv


# path to dir with program
path_to_dir = '/home/fedor/Programs/my/stv/'
input_bed_path = argv[1]


def get_max_mon(live_mons_bed):
    all_numbers = [mon[3].split('.')[1] for mon in live_mons_bed]
    # delete hybrids which errors in max()
    only_int = []
    for n in all_numbers:
        if n.isdigit():
            only_int.append(int(n))
    max_mon = str(max(only_int))
    return max_mon

def stv_namer(live_stv_name, mons_numbers, strand):
    if strand == '+':
        stv_name = mons_numbers[0]
        if mons_numbers[0].isdigit():
            i_prev = mons_numbers[0]
        elif '/' in mons_numbers[0] and 'S' not in mons_numbers[0]: # first mon is hybrid
            i_prev = 25
        else: # 8&12 in chr18
            i_prev = 25
        status = 'Closed'
        for i in mons_numbers[1:]:
            if i.isdigit():
                if int(i) == int(i_prev) + 1:
                    i_prev = i
                    status = 'ToBeClosed'
                else:
                    if status == 'ToBeClosed':
                        stv_name += '-{}_{}'.format(i_prev, i)
                    else:
                        stv_name += '_{}'.format(i)
                    i_prev = i
                    status = 'Closed'
            # hybrid like '4/7'
            elif '/' in i and 'S' not in i: 
                if status == 'ToBeClosed':
                    stv_name += '-{}_{}'.format(i_prev, i)
                else:
                    stv_name += '_{}'.format(i)
                i_prev = 25
                status = 'Closed'
            # 8&12 in chr18
            else: 
                if status == 'ToBeClosed':
                    stv_name += '-{}_{}'.format(i_prev, i)
                else:
                    stv_name += '_{}'.format(i)
                status = 'Closed'
                i_prev = '20'
        if status == 'ToBeClosed':
            stv_name += '-{}'.format(i)
    # reversed. strand == '-'
    else: 
        i_prev = 25
        stv_name = mons_numbers[0]
        if mons_numbers[0].isdigit():
            i_prev = mons_numbers[0]
        elif '/' in mons_numbers[0] and 'S' not in mons_numbers[0]: # first mon is hybrid
            i_prev = 25
        status = 'Closed'
        for i in mons_numbers[1:]:
            if i.isdigit():
                if int(i) == int(i_prev) - 1:
                    i_prev = i
                    status = 'ToBeClosed'
                else:
                    if status == 'ToBeClosed':
                        stv_name += '-{}_{}'.format(i_prev, i)
                    else:
                        stv_name += '_{}'.format(i)
                    i_prev = i
                    status = 'Closed'
            # hybrid like '4/7'
            elif '/' in i and 'S' not in i: 
                if status == 'ToBeClosed':
                    stv_name += '-{}_{}'.format(i_prev, i)
                else:
                    stv_name += '_{}'.format(i)
                i_prev = 25
                status = 'Closed'
            else:
                if status == 'ToBeClosed':
                    stv_name += '-{}_{}'.format(i_prev, i)
                else:
                    stv_name += '_{}'.format(i)
                status = 'Closed'
                i_prev = '20'
        if status == 'ToBeClosed':
            stv_name += '-{}'.format(i)
    stv_name = '{}.{}'.format(live_stv_name, stv_name)
    return(stv_name)

def cen18_correcter(stv_bed):
    corrected_bed = []
    for line in stv_bed:
        stv_name = line[3]
        if stv_name != 'S2C18H1L.8&12':
            corrected_bed.append(line)
        else:
            corrected_stv = corrected_bed.pop()
            corrected_stv[2] = line[2]
            corrected_stv[7] = line[2]
            corrected_stv[3] += '_8&12'
            corrected_bed.append(corrected_stv)
    return corrected_bed

def cen1_5_19_compressor(stv_bed):
    # func is so fat because there are two case: normal (_6/4_5) and inversion in cen1 (5_6/4)
    corrected_bed = []
    # first loop _6/4_5_6/4_5... --> (X){n}, n > 1
    for line in stv_bed:
        if line[0] == 'chr1' and line[5] == '-':
            if search(r'[-_](5_6/4_){2,}', line[3]):
                old_name = line[3]
                new_name = inv_repeat_squizzer(old_name)
                line[3] = new_name
                corrected_bed.append(line)
            else:
                corrected_bed.append(line)
        else:
            if search(r'(_6/4_5){2,}[-_]', line[3]):
                old_name = line[3]
                new_name = repeat_squizzer(old_name)
                line[3] = new_name
                corrected_bed.append(line)
            else:
                corrected_bed.append(line)
    # one more loop to _6/4_5 --> (X){1}
    recorrected_bed = []
    for line in corrected_bed:
        if line[0] == 'chr1' and line[5] == '-':
            search_result = search(r'([-_])5_6/4_', line[3])
            if search_result:
                connector = search_result.groups()[0]
                old_name = line[3]
                new_name = sub(r'[-_]5_6/4_', '{}(X){}1{}'.format(connector, '{', '}'), old_name)
                line[3] = new_name
                recorrected_bed.append(line)
            else:
                recorrected_bed.append(line)
        else:
            search_result = search(r'_6/4_5([-_])', line[3])
            if search_result:
                connector = search_result.groups()[0]
                old_name = line[3]
                new_name = sub(r'_6/4_5[-_]', '(X){}1{}{}'.format('{', '}', connector), old_name)
                line[3] = new_name
                recorrected_bed.append(line)
            else:
                recorrected_bed.append(line)
    # (X) --> (_6/4_5)
    for line in recorrected_bed:
        if line[0] == 'chr1' and line[5] == '-':
            line[3] = line[3].replace('(X)', '(5_6/4_)')
        else:
            line[3] = line[3].replace('(X)', '(_6/4_5)')
    return recorrected_bed

def repeat_squizzer(name):
    result = search(r'(_6/4_5){2,}[-_]', name)
    all_repeats = result.group()
    n_repeats = int((len(all_repeats)-1)/6)
    connector = all_repeats[-1]
    name = name.replace(all_repeats, '(X){}{}{}{}'.format('{', n_repeats, '}', connector), 1)
    if search(r'(_6/4_5){2,}[-_]', name):
        name = repeat_squizzer(name)
    return name

def inv_repeat_squizzer(name):
    result = search(r'[-_](5_6/4_){2,}', name)
    all_repeats = result.group()
    n_repeats = int((len(all_repeats)-1)/6)
    connector = all_repeats[0]
    name = name.replace(all_repeats, '{}(X){}{}{}'.format(connector, '{', n_repeats, '}'), 1)
    if search(r'[-_](5_6/4_){2,}', name):
        name = inv_repeat_squizzer(name)
    return name

def write_stat_file(uniq_stv_counts, file):
    uniq_stv_counts['total'] = sum(uniq_stv_counts.values())
    with open(file, 'w') as f:
        for key in uniq_stv_counts:
            f.write('{}\t{}\n'.format(key, uniq_stv_counts[key]))

def add_stv_number(stv_bed):
    cnt = 0
    for line in stv_bed:
        cnt += 1
        line[3] = '#{}:{}'.format(cnt, line[3])
    return stv_bed

def print_bed(bed):
    for line in bed:
        print('\t'.join(line))


# parse input bed
input_bed = []
with open(input_bed_path) as bed:
    for line in bed:
        if line[:5] != 'track': # skip header
            input_bed.append(line.split())

# dict chrN: [[mon_names], [strands]]
chr_dict = {}
for line in input_bed:
    chr = line[0]
    mon = line[3].split('.')[0]
    strand = line[5]
    if chr in chr_dict:
        chr_dict[chr][0].append(mon)
        chr_dict[chr][1].append(strand)
    else:
        chr_dict[chr] = [[mon], [strand]]
# dict chrN: [live_mon, strand]
for chr in chr_dict:
    chr_dict[chr] = [Counter(chr_dict[chr][0]).most_common(1)[0][0], Counter(chr_dict[chr][1]).most_common(1)[0][0]]
# dict with live mons coordinats from t2t_cenAnnotation.v2.021921.bed
live_mons_coords = {}
with open('{}cenAnnotation_live.bed'.format(path_to_dir)) as f:
    for line in f:
        line = line.split()
        chr = line[0]
        start = line[1]
        end = line[2]
        if chr not in live_mons_coords:
            live_mons_coords[chr] = [[start, end]]
        else:
            live_mons_coords[chr].append([start, end])


# MAIN LOOP
#del chr_dict['chr1']
stv_bed = []
for chr in chr_dict:
    #print(chr)
    live_mon_name = chr_dict[chr][0]
    # BUG cen21 live nome is not the most common now
    if chr == 'chr21':
        live_mon_name = 'S2C13/21H1L'
    # get live mons
    # find start and end of live cen and take everything between them 
    live_mons = []
    for line in input_bed:
        line_chr = line[0]
        start = int(line[1])
        end = int(line[2])
        for live_region in live_mons_coords[chr]:
            if line_chr == chr and start >= int(live_region[0]) and end <= int(live_region[1]):
                live_mons.append(line)
    #print(live_mons_coords[chr])
    # get max mon
    only_live_mons = []
    for line in live_mons:
        if line[3].split('.')[0] == live_mon_name:
            only_live_mons.append(line)
    max_mon = get_max_mon(only_live_mons)
    only_live_mons = []

    # fill stv list
    stvs = []
    mons_numbers = []
    n_prev = max_mon
    start = live_mons[0][1]
    end = live_mons[0][2]
    is_prev_max = True
    for line in live_mons:
        #print(line)
        strand = line[5]
        if line[3].split('.')[0] != live_mon_name: # non live mon
            mons_numbers.append(line[3])
            continue
        n = line[3].split('.')[1]
        mons_numbers.append(n)
        if strand == '+':
            # max mon (or max mon last in hybrid) THAN cut after it
            if n == max_mon or n[-2:] == '/{}'.format(max_mon) or n[-3:] == '/{}'.format(max_mon):
                end = line[2]
                stv_name = stv_namer(live_mon_name, mons_numbers, strand)
                stvs.append([chr, start, end, stv_name, '0', strand, start, end, '0,0,0'])
                start = line[2]
                stv_name = []
                mons_numbers = []
                n_prev = n
                is_prev_max = True
                continue
            # first mon (or first mon is 1st in hybrid) AND previous wasn't max (or hybrid with max) OR big gap before THAN cut before it
            elif ((n == '1' or n[:2] == '1/') and is_prev_max == False) or int(line[1])-int(end) > 800:
                mons_numbers.pop()
                stv_name = stv_namer(live_mon_name, mons_numbers, strand)
                stvs.append([chr, start, end, stv_name, '0', strand, start, end, '0,0,0'])
                start = line[1]
                stv_name = []
                mons_numbers = [n]
            end = line[2]
            n_prev = n
            is_prev_max = False
        else: # strand == '-'
            # first mon (or 1st is the last in hybrid) THAN cut after it
            if n == '1' or n[-2:] == '/1':
                end = line[2]
                stv_name = stv_namer(live_mon_name, mons_numbers, strand)
                stvs.append([chr, start, end, stv_name, '0', strand, start, end, '0,0,0'])
                start = line[2]
                stv_name = []
                mons_numbers = []
                n_prev = n
                is_prev_max = True
                continue
            # fix BUG with cen1 inversion 
            elif chr == 'chr1' and n == max_mon and is_prev_max == False:
                mons_numbers.pop()
                stv_name = stv_namer(live_mon_name, mons_numbers, strand)
                stvs.append([chr, start, end, stv_name, '0', strand, start, end, '0,0,0'])
                start = line[1]
                stv_name = []
                if n == max_mon:
                    mons_numbers = [max_mon]
                else:
                    mons_numbers = [n]
            # max mon (or max is 1st in hybrid) AND previous wasn't the first
            elif ((n == max_mon or n[:2] == '{}/'.format(max_mon) or n[-3:] == '{}/'.format(max_mon)) and is_prev_max == False and chr != 'chr1') or int(line[1])-int(end) > 200:
                mons_numbers.pop()
                stv_name = stv_namer(live_mon_name, mons_numbers, strand)
                stvs.append([chr, start, end, stv_name, '0', strand, start, end, '0,0,0'])
                start = line[1]
                stv_name = []
                if n == max_mon:
                    mons_numbers = [max_mon]
                else:
                    mons_numbers = [n]
            end = line[2]
            n_prev = n
            is_prev_max = False
    if len(mons_numbers) > 0:
        stv_name = stv_namer(live_mon_name, mons_numbers, strand)
        stvs.append([chr, start, end, stv_name, '0', strand, start, end, '0,0,0'])

    stvs_corrected = cen18_correcter(cen1_5_19_compressor(stvs)) 
    stv_bed += stvs_corrected

print_bed(stv_bed)

