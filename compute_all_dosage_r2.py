


import os
import sys
import math
import multiprocessing as mp
from string import maketrans


DOSAGE_OFFSET = 3


###
def read_params():
    global args
    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == '--straight-only':
            args['OOS'] = True
            i += 1
            continue
        if sys.argv[i] == '--dosages':
            args['dosages'] = sys.argv[i + 1]
            i += 2
            continue
        if sys.argv[i] == '--r2':
            args['r2'] = sys.argv[i + 1]
            i += 2
            continue
        if sys.argv[i] == '--r2-score':
            args['r2score'] = float(sys.argv[i + 1])
            i += 2
            continue
        if sys.argv[i] == '--fam' :
            args['fam'] = sys.argv[i + 1]
            i += 2
            continue
        if sys.argv[i] == '--only-founders':
            args['OF'] = True
            i += 1
            continue
        if sys.argv[i] == '--only-hla-markers':
            args['OHM'] = True
            i += 1
            continue
        if sys.argv[i] == '--simple-threads' or sys.argv[i] == '-t':
            args['NR_PROCESSES'] = int(sys.argv[i + 1])
            i += 2
            continue
        if sys.argv[i] == '--out-f' or sys.argv[i] == '-o':
            args['out_file'] = sys.argv[i + 1]
            i += 2
            continue
        print "### ERROR: unrecognized argument " + sys.argv[i]
        exit(1)
    if not args.has_key('out_file') :
        print "### ERROR: output files prefix not supplied"
        exit(1)

###
def load_fam():
    _fam_file = open(args['fam'])
    _fam = [line.strip().split(' ') for line in _fam_file]
    
    return _fam
###
def load_dosages():
    _dosages_file = open(args['dosages'], 'r')
    _r2_file = open(args['r2'], 'r')
    
    _dosages = []
    while True:
        _line_d = _dosages_file.readline().strip()
        _line_r = _r2_file.readline().strip()
        if _line_d == '' or _line_r == '':
            if _line_d != '' or _line_r != '':
                print '### ERROR: dosage and r2 files did not have the same number of lines ... stopping'
                quit(1)
            break
        _toks_d = _line_d.split('\t')
        _toks_r = _line_r.split('\t')
        if _toks_d[0] != _toks_r[0]:
            print '### ERROR: marker names do not match somewhere in the two files (dosage and r2 file respectively): dosage marker=' + _toks_d[0] + '\t r2 marker=' + _toks_r[0]
            quit(1)
        if _toks_r[1].lower() == "nan":
            continue
        if float(_toks_r[1]) < args['r2score']:
            continue
        if args['OHM']:
            if 'HLA' not in _toks_d[0] and 'AA' not in _toks_d[0] and 'SNP' not in _toks_d[0] and 'INS' not in _toks_d[0]:
                continue
        _dosages.append([ float(_toks_d[_i]) if _i >= DOSAGE_OFFSET else _toks_d[_i] for _i in range(len(_toks_d)) ])
        
    return _dosages
###
def mean(_l):
    return float(sum(_l))/float(len(_l))
###
def corr(_x, _y):
    _mean_y = mean(_y)
    _mean_x = mean(_x)
    
    _cov_reduced = 0
    _var_x_reduced = 0
    _var_y_reduced = 0
    for _i in range(len(_x)):
        _cov_reduced += (_x[_i] - _mean_x) * (_y[_i] - _mean_y)
        _var_x_reduced += (_x[_i] - _mean_x) * (_x[_i] - _mean_x)
        _var_y_reduced += (_y[_i] - _mean_y) * (_y[_i] - _mean_y)
        
    return _cov_reduced / math.sqrt(_var_x_reduced * _var_y_reduced)
    
    
###
def match_indivs(_i, _j):
    _corr = corr([ dosages[_ii][_i] for _ii in range(len(dosages)) ], [ dosages[_ii][_j] for _ii in range(len(dosages)) ])
    
    _ret =  fam[_i - DOSAGE_OFFSET][0:5]
    _ret.extend(fam[_j - DOSAGE_OFFSET][0:5])
    _ret.append(len(dosages))
    _ret.append(_corr)
    return _ret
###


###################################################
### void main(int argc, char **argv)            ###


args = {}
args['OOS'] = False
args['OF']  = False
args['OHM']  = False
args['NR_PROCESSES'] = 1
read_params()
#
fam = load_fam()
print '\t *** read all individuals ..'
dosages = load_dosages()
print '\t *** loaded all dosage data ..'
#
inputs = []
for i in range(len(fam) - 1):
    for j in range(i + 1, len(fam)):
        if args['OOS'] and fam[i][-2] == fam[j][-2]:
            continue
        if args['OF'] and ( fam[i][2] != '0' or fam[i][3] != '0' or fam[j][2] != '0' or fam[j][3] != '0' ):
            continue
        inputs.append([i + DOSAGE_OFFSET, j + DOSAGE_OFFSET])
#
print '\t *** built all input combinations, starting the parallel computation ..'
pool = mp.Pool(processes=args['NR_PROCESSES'])
results = [pool.apply(match_indivs, args=x) for x in inputs]
#
header = '\t'.join(["fam1", "ind1", "father1", "mother1", "sex1", "fam2", "ind2", "father2", "mother2", "sex2", "Nmarkers","corr"])
out_file = open(args['out_file'], 'w')
out_file.write(header + '\n')
out_file.write('\n'.join([ '\t'.join(map(str,i)) for i in results ]))
out_file.close()















