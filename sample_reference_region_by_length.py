






import os
import sys
import random


NR_SAMPLE_INTERVALS = 1000
REGION_LENGTH = 10000 # 20kb

# centromeres = { '1':[121500000, 128900000], \
#                '2':[90500000, 96800000], \
#                '3':[87900000, 93900000], \
#                '4':[48200000,52700000], \
#                '5':[46100000,50700000], \
#                '6':[58700000, 63300000], \
#                '7':[58000000, 61700000], \
#                '8':[43100000, 48100000], \
#                '9':[47300000, 50700000], \
#                '10':[38000000, 42300000], \
#                '11':[51600000, 55700000], \
#                '12':[33300000, 38200000], \
#                '13':[16300000, 19500000], \
#                '14':[16100000, 19100000], \
#                '15':[15800000, 20700000], \
#                '16':[34600000, 38600000], \
#                '17':[22200000, 25800000], \
#                '18':[15400000, 19000000], \
#                '19':[24400000, 28600000], \
#                '20':[25600000, 29400000], \
#                '21':[10900000, 14300000], \
#                '22':[12200000, 17900000]}

centromeres = { '1':[], \
                '2':[], \
                '3':[], \
                '4':[], \
                '5':[], \
                '6':[], \
                '7':[], \
                '8':[], \
                '9':[], \
                '10':[], \
                '11':[], \
                '12':[], \
                '13':[], \
                '14':[], \
                '15':[], \
                '16':[], \
                '17':[], \
                '18':[], \
                '19':[], \
                '20':[], \
                '21':[], \
                '22':[]}

bad_intervals = { '1':[], \
                '2':[], \
                '3':[], \
                '4':[], \
                '5':[], \
                '6':[], \
                '7':[], \
                '8':[], \
                '9':[], \
                '10':[], \
                '11':[], \
                '12':[], \
                '13':[], \
                '14':[], \
                '15':[], \
                '16':[], \
                '17':[], \
                '18':[], \
                '19':[], \
                '20':[], \
                '21':[], \
                '22':[]}


###
def read_formated_bad_intervals_file(_bad_intervals_file_name, _bad_intervals, _centromeres):
    _bad_intervals_file = open(_bad_intervals_file_name, 'r')
    _bad_intervals_header = _bad_intervals_file.readline()
    for _line in _bad_intervals_file:
        _toks = _line.strip().split('\t')
        if _toks[-1].split('_')[-1] == 'centromere':
            _centromeres[toks[0]] = [int(_toks[1]), int(_toks[2])]
            continue
        _bad_intervals[_toks[0]].append([int(_toks[1]), int(_toks[2])])
    
    return [_bad_intervals_header, _bad_intervals, _centromeres] 

###
def reduce_interval_list(_list): ### takes a list of intervals and merges adjacent intervals that either overlap or stitch perfectly at the borders
    _list = sorted(_list, key=lambda x: (x[0],x[1]), reverse=False)
    i = 0
    while len(_list) > 1:
        if i == len(_list) - 1:
            break
        if _list[i][1] >= _list[i + 1][0]:
            _new_i = [min(_list[i][0], _list[i+1][0]), max(_list[i][1], _list[i+1][1]) ]
            _list.pop(i)
            _list.pop(i)
            _list.insert(i, _new_i)
            continue
        i += 1
    return _list[:]        
###
def sample_region_by_length(_centromeres, _bad_intervals, _chromosome_ends, _length):
    _sampled_chr = random.choice(_centromeres.keys())
    _sampled_arm = random.choice([True, False]) # True = p-arm; False = q-arm
    if _sampled_arm: # p-arm case
        if _centromeres[_sampled_chr][0] - _length - _chromosome_ends[_sampled_chr][0] <= 0: # sort of pointless check
            sys.stderr.write('### Sampling attempt failed(chr_len): sampled_chr={} sampled_chromosome_arm={} arm_start={} arm_end={} desired_interval_length={} \n'.format(_sampled_chr, _sampled_arm, _chromosome_ends[_sampled_chr][0], _centromeres[_sampled_chr][0], _length))
            return []
        _interval_start = -1
        _flag = -1
        for i in range(0,100):
            _interval_start = random.randrange(_chromosome_ends[_sampled_chr][0], _centromeres[_sampled_chr][0] - _length)
            if sum([0 if _interval_start + _length < _it[0] or _interval_start > _it[1] else 1 for _it in _bad_intervals[_sampled_chr]]) == 0:
                _flag = 1
                break
        if _flag == -1:
            sys.stderr.write('### Sampling attempt failed: sampled_chr={} sampled_chromosome_arm={} arm_start= arm_end={} desired_interval_length={} \n'.format(_sampled_chr, _sampled_arm, _chromosome_ends[_sampled_chr][0], _centromeres[_sampled_chr][0], _length))
            return []
    else:
        if _chromosome_ends[_sampled_chr][1] - _length - _centromeres[_sampled_chr][1] <= 0:
            sys.stderr.write('### Sampling attempt failed(chr_len): sampled_chr={} sampled_chromosome_arm={} arm_start={} arm_end={} desired_interval_length={} \n'.format(_sampled_chr, _sampled_arm, _centromeres[_sampled_chr][1], _chromosome_ends[_sampled_chr][1], _length))
            return []
        _interval_start = -1
        _flag = -1
        for i in range(0,100):
            _interval_start = random.randrange(_centromeres[_sampled_chr][1], _chromosome_ends[_sampled_chr][1] - _length)
            if sum([ 0 if _interval_start + _length < _it[0] or _interval_start > _it[1] else 1 for _it in _bad_intervals[_sampled_chr] ]) == 0:
                _flag = 1
                break
        if _flag == -1:
            sys.stderr.write('### Sampling attempt failed: sampled_chr={} sampled_chromosome_arm={} arm_start={} arm_end={} desired_interval_length={} \n'.format(_sampled_chr, _sampled_arm, _centromeres[_sampled_chr][1], _chromosome_ends[_sampled_chr][1], _length))                             
            return []
    return [_sampled_chr, _interval_start, _interval_start + _length -1]
###



##################################################
### void main (int argc, char **argv)



i = 1
out_file_name = ''
while i < len(sys.argv):
    if sys.argv[i] == '--exclude':
        bad_intervals_file_name = sys.argv[i + 1]
        i += 2
        continue
    if sys.argv[i] == '--out':
        out_file_name = sys.argv[i + 1]
        i += 2
        continue
    if sys.argv[i] == '--count':
        NR_SAMPLE_INTERVALS = int(sys.argv[i + 1])
        i += 2
        continue
    if sys.argv[i] == '--len':
        REGION_LENGTH = int(sys.argv[i + 1])
        i += 2
        continue


bad_intervals_header, bad_intervals, centromeres = read_formated_bad_intervals_file(bad_intervals_file_name, dict(bad_intervals), dict(centromeres))


for _chr in bad_intervals.keys():
    bad_intervals[_chr] = reduce_interval_list(bad_intervals[_chr])
chromosome_ends = dict([ (i, [min([j[1] for j in bad_intervals[i]]), max([j[0] for j in bad_intervals[i]])] ) for i in bad_intervals.keys() ])



# print centromeres
# print "*" * 40
# print chromosome_ends
# print "*" * 40
# print bad_intervals
# exit(1)

## SAMPLE!
intervals = []
for i in range(0,NR_SAMPLE_INTERVALS):
    if i % 100 == 0:
        sys.stderr.write( '...sampled {} intervals\n'.format(i) )
    sampled_interval = []
    while len(sampled_interval) == 0:
        sampled_interval = sample_region_by_length(centromeres, bad_intervals, chromosome_ends, REGION_LENGTH)
    intervals.append(sampled_interval[:])
    
if out_file_name != '':
    out_file = open(out_file_name, 'w')
    out_file.write( '\n'.join([ '{0[0]}:{0[1]}-{0[2]}'.format(i) for i in intervals ]) )
    out_file.close()
else:
    print '\n'.join([ '{0[0]}:{0[1]}-{0[2]}'.format(i) for i in intervals ])




