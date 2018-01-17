

import os
import sys
import subprocess
import multiprocessing as mp

OF = False
MT = False
OOS = False
NR_PROCESSES = 1
OMIT_NC = True
POP_SEP = '\t'

###
def get_population(_pop_name, _pop_f):
    _pop = open(_pop_f, 'r')
    return [_pop_name, [_line.strip().split(POP_SEP) for _line in _pop]]
###
def indiv_in_population(_id, _pop, _of):
    for _i in _pop[1]:
        if _id.lower() == _i[1].lower():
            if _of and (_i[2] != '0' or _i[3] != '0'):
                return False
            return True
    return False
    
###
def get_sex (_ind, _ped):
    for _i in _ped[1]:
        if _ind == _i[1]:
            return _i[4]
    return '-1'
###
def parse_vcf(_vcf):
    _header = ''
    _cols = []
    _records = [] 
    for _line in _vcf:
        if len(_line) < 2:
            print '$$$ Weird input VCF line ...'
            break
        if _line[0:2] == '##':
            _header += _line
            continue
        if _line[0:2] == '#C':
            _cols = _line.strip().split('\t')
            continue
        break
        _records.append( _line.strip().split('\t') )
    
    return _header, _cols, _records
###

# the sequential ifs in this method are not (all) mutually exclusive -- remember order
def compare_genotypes(_a, _b):
    # _ret = [total_nc, total_c, identical_nc, identical_c]
    _ret = [0,0,0,0]
    if _a[0] == '0' and _a[1] == '0' or _b[0] == '0' and _b[1] == '0' :
        return _ret
    
    if (_a[0] == _b[0] and _a[1] == _b[1]) or (_a[1] == _b[0] and _a[0] == _b[1]):
        if _a[0] == '0' or _a[1] == '0':
            _ret = [1,0,1,0]
        else:
            _ret = [1,1,1,1]
        return _ret
    
    if _a[0] != '0' and _a[1] != '0' and _b[0] != '0' and _b[1] != '0' :
        return [1,1,0,0]
    if _a[0] != '0' and (_a[0] == _b[0] or _a[0] == _b[1]) or _a[1] != '0' and (_a[1] == _b[0] or  _a[1] == _b[1]):
        return [1,0,1,0]
    
    if (_a[0] != '0' or _a[1] != '0') and (_b[0] != '0' or _b[1] != '0') : # this if is useless
        return [1,0,0,0]
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
    print _ret 
    return _ret
                
###
# uses the global variables
def match_indivs(_p, _i, _j):
    # _stats = [total_gt_nc, total_gt_c, identical_gt_nc, identical_gt_c, total_al_c, identical_al_c]
    _stats = [0,0,0,0,0,0]
    for _it in range(6,len(indivs[_p][_i][0][:]),2):
        _gt_stats = compare_genotypes([ indivs[_p][_i][0][_it], indivs[_p][_i][0][_it + 1] ], [ indivs[_p][_j][0][_it], indivs[_p][_j][0][_it + 1] ] )
        _al_stats = compare_alleles([ indivs[_p][_i][0][_it], indivs[_p][_i][0][_it + 1] ], [ indivs[_p][_j][0][_it], indivs[_p][_j][0][_it + 1] ] )
        _it_stats = _gt_stats[:]
        _it_stats.extend(_al_stats[:])
        _stats = [ _stats[_ii] + _it_stats[_ii] for _ii in range(0, len(_stats)) ]
        
    _ret = [_p]
    _ret.extend(indivs[_p][_i][0][0:5])
    _ret.extend(indivs[_p][_j][0][0:5])
    _ret.extend(_stats)
    _ret.extend([ float(_stats[2])/_stats[0], float(_stats[3])/_stats[1], float(_stats[5])/_stats[4] ])
    
    print 'computed combination: ' + '\t'.join([ indivs[_p][_i][0][1], indivs[_p][_j][0][1] ])
    return _ret   
###
def read_params():
    global args
    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == '--straight_only':
            args['OOS'] = True
            i += 1
            continue
        if sys.argv[i] == '--mega_threaded':
            args['MT'] = True
            i += 1
            continue
        if sys.argv[i] == '--tmp_dirs':
            args['tempdir'] = sys.argv[i + 1]
            args['logs'] = sys.argv[i + 2]
            i += 3
            continue
        if sys.argv[i] == '--comb_file':
            args['combination_list'] = sys.argv[i + 1]
            i += 2
            continue
        if sys.argv[i] == '--carry_horse':
            args['array_job'] = sys.argv[i + 1]
            args['horse_script'] = sys.argv[i + 2]
            i += 3
            continue
        if sys.argv[i] == '--job_params':    
            args['jobName'] = sys.argv[i + 1]
            args['jobRt'] = sys.argv[i + 2]
            i += 3
            continue
        if sys.argv[i] == '--police_job':
            args['police_job'] = sys.argv[i + 1]
            i += 2
            continue
        if sys.argv[i] == '--vcf':
            args['vcf_f'] = sys.argv[i + 1]
            i += 2
            continue
        if sys.argv[i] == '--population' :
            args['population'].append([sys.argv[i + 1], sys.argv[i + 2]])
            # populations.append(get_population(sys.argv[i + 1], sys.argv[i + 2]))
            i += 3
            continue
        if sys.argv[i] == '--only_founders':
            args['OF'] = True
            i += 1
            continue
        if sys.argv[i] == '--simple_threads' or sys.argv[i] == '-t':
            args['NR_PROCESSES'] = int(sys.argv[i + 1])
            i += 2
            continue
        if sys.argv[i] == '--out_f' or sys.argv[i] == '-o':
            args['out_prefix'] = sys.argv[i+1]
            i += 2
            continue
        print "### ERROR: unrecognized argument " + sys.argv[i]
        exit(1)
    if not args.has_key('out_prefix') :
        print "### ERROR: output files prefix not supplied"
        exit(1)

###





###################################################
### void main(int argc, char **argv)            ###

# ped_f = sys.argv[1]
populations = []
args = {'population' : []}
args ['OF'] = False

print ' ----> reading params ...'
read_params()
print str(args)


# populations = [ [POP_NAME, [IND1, IND2,...] ], ... ]
for i in args['population']:
    populations.append(get_population(i[0], i[1]))

    

print ' ----> processed params ...'

if args.has_key('vcf_f') == False :
    print '### ERROR: No input VCF file (containing sample genotypes) passed to the script ...'
    quit(1)
vcf = open(args['vcf_f'], 'r')

# indivs = [ [[IND1_pop1_name, IND1_pop1_colIndex], [IND2_pop1,IND2_pop1_colIndex]], ... ]
indivs = []
for _i in range(len(populations)): # because no one's perfect
    indivs.append([])


vcf_header, vcf_colnames, vcf_data = parse_vcf(vcf)
col_index = 9


# print '$$$ Processing {} population(s) independently: {}'.format( len(populations), populations )

for col in vcf_colnames[9:]: # i.e.: each col = one individual
    col_idx = vcf_colnames.index(col)
    
    for i in range(len(populations)):
        if indiv_in_population(col, populations[i], args['OF']):
            indivs[i].append([col, col_idx])
            print str(populations[i][0]) + '<--- ' + col
            break
print " ----> distributed individuals by population; we have " + str(len(indivs)) + " populations with " + ','.join((map(str,map(len,indivs[:])))) + " individuals respectively..."
#
inputs = []
if args.has_key('OOS') and args['OOS'] == True:
    print " ----> only computing similarity between diferent sex individuals ..."
for p in range(len(indivs)):
    if len(indivs[p]) == 0:
        print '### Warning: we found 0 individuals in the input file belonging to population {}'.format(indivs[p][0])
        continue
    for i in range(len(indivs[p])):
        for j in range(i,len(indivs[p])):
            if args.has_key('OOS') and args['OOS'] == True:
                if get_sex(indivs[p][i][0], populations[p]) != get_sex(indivs[p][j][0], populations[p]): 
                    inputs.append([p,i,j])
            else:
                inputs.append([p,i,j])
print ' ----> computed all individual combinations needed({}) ... '.format( str(len(inputs)) )
#

### TODO!!!::: many things are hard coded across 3 scripts here --- clean it up at some point!!!!!
if args.has_key('tempdir') and args.has_key('logs'):
    print ' ----> going MEGA THREADED mode ...'
    # jobName = 'mateGoNLwg'
    vcf_name = args['vcf_f'].split('/')[-1]
    comb_list = open(args['combination_list'], 'w')
    for it in inputs:
        comb_list.write(' '.join([ str(indivs[it[0]][it[1]][1]), str(indivs[it[0]][it[2]][1]), str(indivs[it[0]][it[1]][0]), str(indivs[it[0]][it[2]][0]) ]) + '\n')
    comb_list.close()
    # -o ' + logs + '/' + ped_f_name +  '.o -e '  + logs + '/' + ped_f_name + '.e
    cmd = 'qsub -N ' + args['jobName'] + ' -l h_rt=' + args['jobRt'] + ',h_vmem=6G -cwd -pe threaded 1 -t 1-' + str(len(inputs)) + ' ' + args['array_job'] + ' ' + ' '.join([ args['combination_list'], args['vcf_f'], args['tempdir'], args['horse_script'], args['population'][0][1] ]) # this last argument is a hack; I got lost in my own stuff with it
    #os.chdir('/hpc/cog_bioinf/data/mcretu/mating/logs/temp/')
    print ' ----> Horse command: ' + cmd + '\n'
    os.chdir(args['logs'] + '/')
    os.system(cmd)
    print ' ----> submitted all computation jobs ...'
    if args.has_key('police_job'):
        cmd_police = 'qsub -N police.' + args['jobName'] + ' -hold_jid ' + args['jobName'] + ' -l h_rt=04:00:00' + ',h_vmem=6G -cwd -pe threaded 1 '  + args['police_job'] + ' ' + ' '.join([ args['combination_list'], args['vcf_f'], args['tempdir'], args['horse_script'], args['logs'], args['jobName'], args['jobRt'], args['population'][0][1] ])
        os.system(cmd_police)
    
    
    out_file_name = args['out_prefix'] # + '.' + populations[0][0] + '.txt'
    out_file = open(args['out_prefix'], 'w') # + '.' + populations[0][0] + '.txt', 'w')
    header = '\t'.join(["fam1", "ind1", "father1", "mother1", "sex1", "fam2", "ind2", "father2", "mother2", "sex2", "totalNCgt", "totalCgt", "identicalNCgt", "identicalCgt", "totalCal", "identicalCal", "qNCgt", "qCgt", "qCal" ])
    out_file.write(header + '\n')
    out_file.close()
    
    # cmd2 = 'echo "cat ' + tempdir + '/tmp.combination.file* >> ' + out_file_name + ' " | qsub -N gather -hold_jid ' + jobName + ' -l h_rt=1:00:00,h_vmem=5G -cwd -pe threaded 1 -o ' + logs + '/gather.o -e ' + logs + '/gather.e'  
    # os.system(cmd2)
    
    quit(0)
    


pool = mp.Pool(processes=NR_PROCESSES)
results = [pool.apply(match_indivs, args=x) for x in inputs]

out_files = []
header = '\t'.join(["fam1", "ind1", "father1", "mother1", "sex1", "fam2", "ind2", "father2", "mother2", "sex2", "totalNCgt", "totalCgt", "identicalNCgt", "identicalCgt", "totalCal", "identicalCal", "qNCgt", "qCgt", "qCal" ])
for _t in range(len(populations)):
    _temp = open(out_prefix + '.' + populations[_t][0] + '.txt', 'w')
    _temp.write(header + '\n')
    out_files.append(_temp)

for _r in results:
    out_files[_r[0]].write('\t'.join(map(str,_r[1:])) + '\n')
for _f in out_files:
    #_f.write("TEST")
    _f.close()

    



    
    
#print out
#print "*" * 50
#print err




