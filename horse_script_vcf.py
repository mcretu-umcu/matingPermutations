

import os
import sys
from string import maketrans

NO_CALL = '.'
PED_SEP = '\t'

###
# the sequential ifs in this method are not (all) mutually exclusive -- so order matters!!!
def compare_genotypes(_a, _b):
    # _ret = [total_nc, total_c, identical_nc, identical_c]
    _ret = [0.0,0.0,0.0,0.0]
    if _a[0] == NO_CALL and _a[1] == NO_CALL or _b[0] == NO_CALL and _b[1] == NO_CALL :
        return _ret
    #
    if (_a[0] == _b[0] and _a[1] == _b[1]) or (_a[1] == _b[0] and _a[0] == _b[1]): # alleles are identical (with at most one no-call allele each of them)
        if _a[0] == _a[1]:
            _ret = [1.0,1.0,1.0,1.0]
        else:
            if _a[0] == NO_CALL or _a[1] == NO_CALL:
                # print "### DBG: gt counts should not be equal"
                # print '\t a={}; b={}'.format(_a, _b)
                _ret = [1.0,0.0,0.5,0.0]
            else:
                _ret = [1.0,1.0,0.5,0.5] # this is the HET - HET identical case -> currently consider it 0.5 identical
        return _ret
    #
    if _a[0] != NO_CALL and _a[1] != NO_CALL and _b[0] != NO_CALL and _b[1] != NO_CALL : # at least one allele is different in the 2 GTs by now
        if _a[0] == _b[0] or _a[0] == _b[1] or _a[1] == _b[0] or _a[1] == _b[1]:
            return [1.0,1.0,0.5,0.5]
        return [1.0,1.0,0.0,0.0]
    #
    if _a[0] != NO_CALL and (_a[0] == _b[0] or _a[0] == _b[1]) or _a[1] != NO_CALL and (_a[1] == _b[0] or _a[1] == _b[1]): # there is at least one null allele in at least one of the GTs by now
        print "### DBG: gt counts should not be equal 2"
        print '\t a={}; b={}'.format(_a, _b)
        return [1.0,0.0,0.5,0.0]
    #
    if (_a[0] != NO_CALL or _a[1] != NO_CALL) and (_b[0] != NO_CALL or _b[1] != NO_CALL) : # this if is the else of the first if
        # print "### DBG: Something is wrong, reached the useless if-statement"
        # print '\t a={}; b={}'.format(_a, _b)
        return [1.0,0.0,0.0,0.0]
    else:
        print '### DBG: Very wrong; reached dead code'
    
###
def compare_alleles(_a, _b):
    # _ret = [total_c_al, _identical_c_al]
    _ret = [0,0]
    _ret[0] = min( len([ _i for _i in _a if _i != NO_CALL ]), len([ _i for _i in _b if _i != NO_CALL ]) )
    # at least one of the genotypes has only NO_CALL alleles -- not comparable
    if _ret[0] == 0 :
        return [0,0]
    # there is only 1 allele that is not a NO_CALL for BOTH genotypes
    if _ret[0] == 1 :
        if _a[0] != NO_CALL and (_a[0] == _b[0] or _a[0] == _b[1]) or _a[1] != NO_CALL and ( _a[1] == _b[0] or _a[1] == _b[1]) :
            return [1,1]
        else:
            return [1,0]
    # both alleles, for both genotypes are called
    if _ret[0] == 2 :
        if _a[0] == _b[0]:
            _ret[1] += 1
            if _a[1] == _b[1]:
                _ret[1] += 1
        else:
            if _a[0] == _b[1]:
                _ret[1] += 1
                if _a[1] == _b[0]:
                    _ret[1] += 1
            else:
                if _a[1] == _b[0]:
                    _ret[1] += 1
                else:
                    if _a[1] == _b[1]:
                        _ret[1] += 1
    # print _ret
    if _ret[1] > 2:
        print '### BUG BUG BUG: found more than 2 identical alleles between two genotypes: {} - {}'.format(_a, _b)
    return _ret
                
###
def compare_alleles_no_phase(_a, _b):
    # _ret = [total_c_al, _identical_c_al]
    _ret = [0,0]
    _ret[0] = min( len([ _i for _i in _a if _i != NO_CALL ]), len([ _i for _i in _b if _i != NO_CALL ]) )
    # at least one of the genotypes has only NO_CALL alleles -- not comparable
    if _ret[0] == 0 :
        return [0,0]
    # there is only 1 allele that is not a NO_CALL for BOTH genotypes
    if _ret[0] == 1 :
        if _a[0] != NO_CALL and (_a[0] == _b[0] or _a[0] == _b[1]) or _a[1] != NO_CALL and ( _a[1] == _b[0] or _a[1] == _b[1]) :
            return [1,1]
        else:
            return [1,0]
    # both alleles, for both genotypes are called
    if _ret[0] == 2 :
        if _a[0] == _b[0]:
            _ret[1] += 1
            if _a[1] == _b[1]:
                _ret[1] += 1
        else:
            if _a[0] == _b[1]:
                _ret[1] += 1
                if _a[1] == _b[0]:
                    _ret[1] += 1
            else:
                if _a[1] == _b[0]:
                    _ret[1] += 1
                else:
                    if _a[1] == _b[1]:
                        _ret[1] += 1
    # print _ret 
    return _ret
    
###
def safe_drop_header(_vcf):
    _header = ''
    _cols = []
    while True:
        _line = _vcf.readline()
        if _line == '':
            break
        if _line[0:2] == '##':
            _header += _line
            continue
        if _line[0:2] == '#C':
            _cols = _line.strip().split('\t')
            break
        print '### WARNING: something is weird in VCF header: ' + _line
    return _header, _cols, _vcf
        
###
def get_gt_as_list (_toks):
    _gt = _toks.split(':')[0].translate(maketrans('/','|'), '')
    return _gt.split('|')
###




##########################################################
####### void main (int argc, char **argv) ###
vcff = sys.argv[1]
outf = sys.argv[2]
pedf = sys.argv[3]
i = int(sys.argv[4])
j = int(sys.argv[5])
indv_i = sys.argv[6]
indv_j = sys.argv[7]
#

vcf = open(vcff, 'r')
header, cols, records = safe_drop_header(vcf)
#
_errs = 0
if cols[i].lower() != indv_i.lower():
    print '### ERROR ERROR BAD BAD: i\'th individual in the vcf: {} does not match the expected individual there: {} --- case insensitive comparison'.format( cols[i], indv_i )
    exit(1)

if cols[j].lower() != indv_j.lower():
    print '### ERROR ERROR BAD BAD: j\'th individual in the vcf: {} does not match the expected individual there: {} --- case insensitive comparison'.format( cols[j], indv_j )
    exit(1)
#
ped = open(pedf, 'r')
_found = 0
info_i = []
info_j = []
for line in ped:
    if _found == 2:
        break
    toks_info = line.strip().split(PED_SEP)
    if toks_info[1].lower() == indv_i.lower():
        info_i = toks_info[0:5]
        _found += 1
        continue
    if toks_info[1].lower() == indv_j.lower():
        info_j = toks_info[0:5]
        _found += 1
        continue
if info_i == []:
    info_i = ['-', indv_i, '-', '-', '-']
if info_j == []:
    info_j = ['-', indv_j, '-', '-', '-']

#
stats = [0,0,0,0,0,0]
while True:
    line = records.readline()
    if line == '':
        break
    toks = line.strip().split('\t')
    _als_i = get_gt_as_list(toks[i])
    _als_j = get_gt_as_list(toks[j])
    if len(_als_i) != 2 or len(_als_j) != 2:
        print '### skipping some weird genotypes: {}:{}, {}:{}'.format(indv_i, toks[i], indv_j, toks[j])
        continue
    gt_stats = compare_genotypes(_als_i, _als_j)
    al_stats = compare_alleles  (_als_i, _als_j)
    gt_stats.extend(al_stats)
    
    stats = [stats[_ii] + gt_stats[_ii] for _ii in range(len(stats))]
#
if stats[0] == 0 or stats[1] == 0 or stats[4] == 0:
    print '### ERROR: some counts were 0 for combination {} <--> {} : '.format(indv_i, indv_j) + ' '.join(map(str,stats[:]))
    exit(1)
#
_ret = info_i[:]
_ret.extend(info_j)
_ret.extend(stats)
_ret.extend([ float(stats[2])/stats[0], float(stats[3])/stats[1], float(stats[5])/stats[4] ])
#
out = open(outf, 'w')
out.write('\t'.join(map(str,_ret)) + '\n')
out.close()
# print 'computed combination: ' + '\t'.join([ indivs[_p][_i][0][1], indivs[_p][_j][0][1] ])
  








