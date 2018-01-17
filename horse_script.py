

import os
import sys



###
# the sequential ifs in this method are not (all) mutually exclusive -- so order matters!!!
def compare_genotypes(_a, _b):
    # _ret = [total_nc, total_c, identical_nc, identical_c]
    _ret = [0.0,0.0,0.0,0.0]
    if _a[0] == '0' and _a[1] == '0' or _b[0] == '0' and _b[1] == '0' :
        return _ret
    
    if (_a[0] == _b[0] and _a[1] == _b[1]) or (_a[1] == _b[0] and _a[0] == _b[1]):
        if _a[0] == _a[1]:
            _ret = [1.0,1.0,1.0,1.0]
        else:
            if _a[0] == '0' or _a[1] == '0':
                print "### DBG: gt counts should not be equal"
                _ret = [1.0,0.0,0.5,0.0]
            else:
                _ret = [1.0,1.0,0.5,0.5]
        return _ret
    
    
    if _a[0] != '0' and _a[1] != '0' and _b[0] != '0' and _b[1] != '0' :
        if _a[0] == _b[0] or _a[0] == _b[1] or _a[1] == _b[0] or _a[1] == _b[1]:
            return [1.0,1.0,0.5,0.5]
        return [1.0,1.0,0.0,0.0]
    
    if _a[0] != '0' and (_a[0] == _b[0] or _a[0] == _b[1]) or _a[1] != '0' and (_a[1] == _b[0] or  _a[1] == _b[1]):
        print "### DBG: gt counts should not be equal 2"
        return [1.0,0.0,0.5,0.0]
    
    if (_a[0] != '0' or _a[1] != '0') and (_b[0] != '0' or _b[1] != '0') : # this if is useless
        print "### DBG: Something is wrong, reached the useless if-statement"
        return [1.0,0.0,0.0,0.0]
###
def compare_alleles(_a, _b):
    # _ret = [total_c_al, _identical_c_al]
    _ret = [0,0]
    _ret[0] = min( len([ _i for _i in _a if _i != '0' ]), len([ _i for _i in _b if _i != '0' ]) )
    # at least one of the genotypes has only NO_CALL alleles -- not comparable
    if _ret[0] == 0 :
        return [0,0]
    # there is only 1 allele that is not a NO_CALL for BOTH genotypes
    if _ret[0] == 1 :
        if _a[0] != '0' and (_a[0] == _b[0] or _a[0] == _b[1]) or _a[1] != '0' and ( _a[1] == _b[0] or _a[1] == _b[1]) :
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
def compare_alleles_no_phase(_a, _b):
    # _ret = [total_c_al, _identical_c_al]
    _ret = [0,0]
    _ret[0] = min( len([ _i for _i in _a if _i != '0' ]), len([ _i for _i in _b if _i != '0' ]) )
    # at least one of the genotypes has only NO_CALL alleles -- not comparable
    if _ret[0] == 0 :
        return [0,0]
    # there is only 1 allele that is not a NO_CALL for BOTH genotypes
    if _ret[0] == 1 :
        if _a[0] != '0' and (_a[0] == _b[0] or _a[0] == _b[1]) or _a[1] != '0' and ( _a[1] == _b[0] or _a[1] == _b[1]) :
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





##########################################################
####### void main (int argc, char **argv) ###
pedf = sys.argv[1]
outf = sys.argv[2]
i = int(sys.argv[3])
j = int(sys.argv[4])
indv_i = sys.argv[5]
indv_j = sys.argv[6]



ped = open(pedf, 'r')
_idx = -1
found = 0
toks_i = []
toks_j = []
for line in ped:
    _idx += 1
    if _idx == i:
        found += 1
        toks_i = line.strip().split(' ')
        if toks_i[1].lower() != indv_i.lower():
            print '### ERROR ERROR BAD BAD: individual at line ' + str(i) + ' different than what we expected there(found - expected): ' + '\t'.join([toks_i[1], indv_i])
            exit(1)
    if _idx == j:
        found += 1
        toks_j = line.strip().split(' ')
        if toks_j[1].lower() != indv_j.lower():
            print '### ERROR ERROR BAD BAD: individual at line ' + str(j) + ' different than what we expected there(found - expected): ' + '\t'.join([toks_j[1], indv_j])
            exit(1)
    if found == 2:
        break
if found != 2:
    print '### ERROR ERROR BAD BAD: did not manage to find both individuals in ped_file, this error should never actually come up ' + str(len(toks_i)) + ' ' + str(len(toks_j))

_stats = [0,0,0,0,0,0]
for _it in range(6,len(toks_i[:]),2):
    _gt_stats = compare_genotypes([ toks_i[_it], toks_i[_it + 1] ], [ toks_j[_it], toks_j[_it + 1] ] )
    _al_stats = compare_alleles([ toks_i[_it], toks_i[_it + 1] ], [ toks_j[_it], toks_j[_it + 1] ] )
    _it_stats = _gt_stats[:]
    _it_stats.extend(_al_stats[:])
    _stats = [ _stats[_ii] + _it_stats[_ii] for _ii in range(0, len(_stats)) ]
        
# _ret = [_p]
_ret = []
_ret.extend(toks_i[0:5])
_ret.extend(toks_j[0:5])
_ret.extend(_stats)
if _stats[0] == 0 or _stats[1] == 0 or _stats[4] == 0:
    print '### ERROR: some counts were 0: ' + ' '.join(map(str,_stats[:]))
    exit(1)
_ret.extend([ float(_stats[2])/_stats[0], float(_stats[3])/_stats[1], float(_stats[5])/_stats[4] ])

out = open(outf, 'w')
out.write('\t'.join(map(str,_ret)) + '\n')
out.close()
# print 'computed combination: ' + '\t'.join([ indivs[_p][_i][0][1], indivs[_p][_j][0][1] ])
  



