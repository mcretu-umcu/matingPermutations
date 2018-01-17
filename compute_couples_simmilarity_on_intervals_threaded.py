


import os
import sys
import random


NR_SAMPLE_INTERVALS = 1000
REGION_LENGTH = 6000000 # 6MB
NR_SNPS = 45000
LENGTH_BUFFER_MULTIPLIER = 0.5
SAMPLING_METHOD = 1
MEGA_THREADED = False

# deprecated !
centromeres_dep = { '1':[121500000, 128900000], \
                '2':[90500000, 96800000], \
                '3':[87900000, 93900000], \
                '4':[48200000,52700000], \
                '5':[46100000,50700000], \
                '6':[58700000, 63300000], \
                '7':[58000000, 61700000], \
                '8':[43100000, 48100000], \
                '9':[47300000, 50700000], \
                '10':[38000000, 42300000], \
                '11':[51600000, 55700000], \
                '12':[33300000, 38200000], \
                '13':[16300000, 19500000], \
                '14':[16100000, 19100000], \
                '15':[15800000, 20700000], \
                '16':[34600000, 38600000], \
                '17':[22200000, 25800000], \
                '18':[15400000, 19000000], \
                '19':[24400000, 28600000], \
                '20':[25600000, 29400000], \
                '21':[10900000, 14300000], \
                '22':[12200000, 17900000]}
# initialize the dictionaries (not sure this is needed)
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
def safe_drop_header(_vcf):
    while True:
        line = _vcf.readline()
        if line[0:2] == '#C':
            _vcf.seek(-1 * len(line), 1)
            return
###
# a chromosome hash is : { 'chromosome' : [ file_index_of_first_snp, first_snp_position_on_chr, file_index_of_last_snp, last_snp_position_on_chr, N_markers_on_chromosome_p_arm, N_markers_on_chromosome_q_arm ] }
def hash_vcf_chromosomes(_vcf_file):
    _pos = _vcf_file.tell()
    _line = _vcf_file.readline()
    _toks = _line.strip().split('\t')[0:9]
    if _toks[0] != '1':
        print '%%% First chromosome is not chromosome 1: {}'.format(_toks[0])
        exit(1)
    _hash_chr = {_toks[0]:[int(_pos), int(_toks[1])]}
    
    _prev_toks = _toks[:]
    _prev_pos  = _vcf_file.tell()
    if int(_toks[1]) <= centromeres[_toks[0]][0]:
        _marker_count_p = 1
        _marker_count_q = 0
    if int(_toks[1]) >= centromeres[_toks[0]][1]:
        _marker_count_p = 0
        _marker_count_q = 1
    
    for _line in iter(_vcf_file.readline,''):
        _toks = _line.strip().split('\t')
        _pos = _vcf_file.tell()
        if _toks[0] == _prev_toks[0]:
            if int(_toks[1]) <= centromeres[_toks[0]][0]:
                _marker_count_p += 1
            if int(_toks[1]) >= centromeres[_toks[0]][1]:
                _marker_count_q += 1
            _prev_toks = _toks[:]
            _prev_pos = _vcf_file.tell()
            continue
        _hash_chr[_prev_toks[0]].extend([_prev_pos, int(_prev_toks[1]), _marker_count_p, _marker_count_q])
        _hash_chr[_toks[0]] = [_pos, int(_toks[1])]
        print '\t### Finished chromosome {} with entry {}'.format(_prev_toks[0], _hash_chr[_prev_toks[0]])
        _prev_toks = _toks[:]
        _prev_pos = _vcf_file.tell()
        if int(_toks[1]) <= centromeres[_toks[0]][0]:
            _marker_count_p = 1
            _marker_count_q = 0
        if int(_toks[1]) >= centromeres[_toks[0]][1]:
            _marker_count_p = 0
            _marker_count_q = 1
    _hash_chr[_prev_toks[0]].extend([_pos, int(_toks[1]), _marker_count_p, _marker_count_q])
    
    return _hash_chr
### DEPRECATED!!!
def sample_region_by_length_dep(_hash_chr):
    _sampled_chr = random.choice(_hash_chr.keys())
    _sample_chromosome_arm = random.choice([True, False]) # True = p; False = q
    if _sample_chromosome_arm:
        if centromeres[_sampled_chr][0] - REGION_LENGTH - _hash_chr[_sampled_chr][1] <= 0:
            print '### Sampling attempt failed: sampled_chr={} sampled_chromosome_arm={} first_snp_pos={} telomere_start={} desired_interval_length={} '.format(_sampled_chr, str(_sample_chromosome_arm), _hash_chr[_sampled_chr][1], centromeres[_sampled_chr][0], REGION_LENGTH)
            return []
        _interval_start = random.randrange( _hash_chr[_sampled_chr][1], centromeres[_sampled_chr][0] - REGION_LENGTH )
        #return [ _sampled_chr, _interval_start, _interval_start + REGION_LENGTH ]
    else:
        if _hash_chr[_sampled_chr][3] - REGION_LENGTH - centromeres[_sampled_chr][1] <= 0:
            print '### Sampling attempt failed: sampled_chr={} sampled_chromosome_arm={} first_snp_pos={} telomere_start={} desired_interval_length={} '.format(_sampled_chr, str(_sample_chromosome_arm), _hash_chr[_sampled_chr][3], centromeres[_sampled_chr][1], REGION_LENGTH)
            return []
        _interval_start = random.randrange( centromeres[_sampled_chr][1], _hash_chr[_sampled_chr][3] - REGION_LENGTH )
    return [ _sampled_chr, _interval_start, _interval_start + REGION_LENGTH ]
###
def sample_region_by_markers(_hash_chr):
    _sampled_chr = random.choice(_hash_chr.keys())
    _sample_chromosome_arm = random.choice([True, False]) # True = p; False = q
    if _sample_chromosome_arm:
        #if centromeres[_sampled_chr][0] - REGION_LENGTH * LENGTH_BUFFER_MULTIPLIER - _hash_chr[_sampled_chr][1] <= 0:
        if _hash_chr[_sampled_chr][-2] < NR_SNPS :    
            print '### Sampling attempt failed: sampled_chr={} sampled_chromosome_arm={} first_snp_pos={} telomere_start={} desired_number_of_markers={} '.format(_sampled_chr, str(_sample_chromosome_arm), _hash_chr[_sampled_chr][1], centromeres[_sampled_chr][0], NR_SNPS)
            return []
        _interval_start_snp = random.randrange( 1, _hash_chr[_sampled_chr][-2] - NR_SNPS )
        #return [ _sampled_chr, _interval_start, _interval_start + REGION_LENGTH ]
    else:
        if _hash_chr[_sampled_chr][3] - REGION_LENGTH * LENGTH_BUFFER_MULTIPLIER - centromeres[_sampled_chr][1] <= 0:
            print '### Sampling attempt failed: sampled_chr={} sampled_chromosome_arm={} first_snp_pos={} telomere_start={} desired_number_of_markers={} '.format(_sampled_chr, str(_sample_chromosome_arm), _hash_chr[_sampled_chr][3], centromeres[_sampled_chr][1], NR_SNPS)
            return []
        _interval_start_snp = random.randrange( 1, _hash_chr[_sampled_chr][-1] - NR_SNPS )
    
    return [ _sampled_chr, _interval_start_snp, _interval_start_snp + NR_SNPS ] # last element not really needed
###
def read_formated_bad_intervals_file(_bad_intervals_file_name, _bad_intervals, _centromeres):
    _bad_intervals_file = open(_bad_intervals_file_name, 'r')
    _bad_intervals_header = _bad_intervals_file.readline()
    for _line in _bad_intervals_file:
        _toks = _line.strip().split('\t')
        if _toks[-1].split('_')[-1] == 'centromere':
            _centromeres[_toks[0]] = [int(_toks[1]), int(_toks[2])]
            continue
        _bad_intervals[_toks[0]].append([int(_toks[1]), int(_toks[2])])
    
    return _bad_intervals_header, _bad_intervals, _centromeres

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
def map_ped_trios_to_vcf(_trios, _vcf_header):
    
    _indiv_mapping = {}
    _found_pairs = []
    for _it in _trios:
        if _it[2] in _vcf_header and _it[3] in _vcf_header:
            _indiv_mapping[_it[2]] = _vcf_header.index(_it[2])
            _indiv_mapping[_it[3]] = _vcf_header.index(_it[3])
            _found_pairs.append(_it[:])
    
    return _found_pairs, _indiv_mapping
###
###
###
####################################################
### void main(int argc, char **argv)



i = 1
vcf_file_name = ''
comb_file_name = ''
police_job = ''
while i < len(sys.argv):
    if sys.argv[i] == '--vcf':
        vcf_file_name = sys.argv[i + 1]
        i += 2
        continue
    if sys.argv[i] == '--ped':
        ped_file_name = sys.argv[i + 1]
        i += 2
        continue
    if sys.argv[i] == '--out-file':
        out_file_name = sys.argv[i + 1]
        i += 2
        continue
    if sys.argv[i] == '--exclude':
        bad_intervals_file_name = sys.argv[i + 1]
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
    if sys.argv[i] == '--nr-snps':
        NR_SNPS = int(sys.argv[i + 1])
        i += 2
        continue
    if sys.argv[i] == '--sampling-method':
        SAMPLING_METHOD = int(sys.argv[i + 1])
        i += 2
        continue
    if sys.argv[i] == '--horse-dirs':
        tempdir = sys.argv[i + 1]
        logs = sys.argv[i + 2]
        i += 3
        continue
    if sys.argv[i] == '--horse-scripts':
        array_job = sys.argv[i + 1]
        horse_script = sys.argv[i + 2]
        i += 3
        continue
    if sys.argv[i] == '--combinationList':
        comb_file_name = sys.argv[i + 1]
        i += 2
        continue
    if sys.argv[i] == '--job-params':
        jobName = sys.argv[i + 1]
        jobRt = sys.argv[i + 2]
        i += 3
        continue
    if sys.argv[i] == '--police':
        police_job = sys.argv[i + 1]
        i += 2
        continue
    if sys.argv[i] == '--mega-threaded':
        MEGA_THREADED = True
        i += 1
        continue
    print '### ERROR: found unspecified parameter: {}'.format(sys.argv[i])
    exit(1)
    

print '*** Parsed input parameters ...'

bad_intervals_header, bad_intervals, centromeres = read_formated_bad_intervals_file(bad_intervals_file_name, dict(bad_intervals), dict(centromeres))
# pass through the VCF and hash the start and end positions for each chromosome's data so that the computing script can find location faster; if the actual computing turns out to be fast we can take this part out, as it takes a while...
vcf_file = open(vcf_file_name, 'r')
safe_drop_header(vcf_file)
vcf_header = vcf_file.readline().strip().split('\t')

hash_chr = hash_vcf_chromosomes(vcf_file)


for _chr in bad_intervals.keys():
    bad_intervals[_chr] = reduce_interval_list(bad_intervals[_chr])
chromosome_ends = dict([ (i, [min([j[1] for j in bad_intervals[i]]), max([j[0] for j in bad_intervals[i]])] ) for i in bad_intervals.keys() ])

# Build a set of trios to be compared
print '*** Identifying trios ...'
ped_file = open(ped_file_name, 'r')
trios = [line.strip().split('\t') for line in ped_file if line.strip().split('\t')[1].lower()[-1] == 'c' ] # or line.strip().split('\t')[1].lower()[-1] == 'd']
trios = map(lambda x : map(lambda y : y.lower(),x), trios)
trios, indiv_mapping = map_ped_trios_to_vcf(trios, vcf_header)

# SAMPLE! ...and build the inputs for the horse script
print '*** Sampling intervals ...'
inputs = []
for i in range(1,NR_SAMPLE_INTERVALS+1):
    sampled_interval = []
    if SAMPLING_METHOD == 1:
        while sampled_interval == []:
            sampled_interval = sample_region_by_length(centromeres, bad_intervals, chromosome_ends, REGION_LENGTH)
    elif SAMPLING_METHOD == 2:
        while sampled_interval == []:
            sampled_interval = sample_region_by_markers(hash_chr)
    
    inputs.append([ sampled_interval[0], sampled_interval[1], sampled_interval[2], hash_chr[sampled_interval[0]][0], hash_chr[sampled_interval[0]][2], SAMPLING_METHOD ])
#    
if comb_file_name != '':
    comb_file = open(comb_file_name, 'w')
    comb_file.write('\n'.join([ '\t'.join(map(str,i[:])) for i in inputs ]))
    comb_file.close()
else:
    print '\n'.join([ '\t'.join(map(str,i[:])) for i in inputs ])

out_header_fields = ['chrom', 'firstPos', 'lastPos', 'nSites']
out_header_fields.extend([ it[0] for it in trios ])
out_header = '\t'.join(out_header_fields)

# submit the actual computation
if MEGA_THREADED:
    print '*** Threading everything out ...'
    command_array = 'qsub -N ' + jobName + ' -l h_rt=' + jobRt + ',h_vmem=6G -cwd -pe threaded 1 -t 1-' + str(len(inputs)) + ' ' + array_job + ' ' + ' '.join([ comb_file_name, vcf_file_name, ped_file_name, tempdir, horse_script ])                                                                        
    os.chdir(logs + '/')
    os.system(command_array)
    print '\t---> submitted all computation jobs ...'
    
    if police_job != '': # because some jobs misteriously evaporate on the hpc
        command_police = 'qsub -N police.' + jobName + ' -hold_jid ' + jobName + ' -l h_rt=04:00:00' + ',h_vmem=6G -cwd -pe threaded 1 '  + police_job + ' ' + ' '.join([ comb_file_name, vcf_file_name, ped_file_name, tempdir, horse_script, logs, jobName, jobRt ])
        os.system(command_police)
        print '\t---> the po\'lice iz on the streets ...'
    
    out_file = open(out_file_name, 'w') 
    out_file.write(out_header + '\n')
    out_file.close()
    exit(0)
    


































