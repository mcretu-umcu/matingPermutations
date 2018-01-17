


import os
import sys
from string import maketrans

###
def safe_drop_header(_vcf):
    while True:
        line = _vcf.readline()
        if line[0:2] == '#C':
            _vcf.seek(-1 * len(line), 1)
            return
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
def extract_sites_subset_length(_vcf_file, _chrom, _start, _stop, _hash_chr_start, _hash_chr_stop):
    _vcf_file.seek(_hash_chr_start, 0)
    _ret_sites = []
    
    _prev_toks = _vcf_file.readline().strip().split('\t')
    if _prev_toks[0] != _chrom:
        print '### ERROR: Found diferent than expected chromosome at hashed position in the vcf: expected={}, found={}'.format(_chrom, prev_toks[0])
        exit(1)
    if int(_prev_toks[1]) >= start and int(_prev_toks[1]) <= _stop:
        _ret_sites.append(_prev_toks[:])
    _flag = False
    for _line in _vcf_file :
        _toks = _line.strip().split('\t')
        if _toks[0] != _chrom:
            break
        if int(_toks[1]) >= _start and int(_toks[1]) <= _stop:
            _ret_sites.append(_toks[:])
            continue
        if int(_toks[1]) > _stop:
            _flag = True
            break
        # print "### WARNING (sample by length) : Something dubious is happening when reading the desired interval: chrom={}, start={}, stop={}, last evaluated position is: \n\t{}".format(_chrom, _start, _stop, '\t'.join(_toks[0:8]))
    if not _flag:
        print "### WARNING (sample by length) : Site selection stopped for another reason than reaching the end of the interval: chrom={}, start={}, stop={}, last evaluated position is: \n\t{}".format(_chrom, _start, _stop, '\t'.join(_toks[0:8]))
    return _ret_sites[:]
###
def extract_sites_subset_count(_vcf_file, _chrom, _start, _stop, _hash_chr_start, _hash_chr_stop):
    _vcf_file.seek(_hash_chr_start, 0)
    _ret_sites = []
    
    _i = 0
    _prev_toks = _vcf_file.readline().strip().split('\t')
    if _prev_toks[0] != _chrom:
        print '### ERROR: Found diferent than expected chromosome at hashed position in the vcf: expected={}, found={}'.format(_chrom, _prev_toks[0])
        exit(1)
    _i += 1
    
    if _i >= _start and _i < _stop:
        _ret_sites.append(_prev_toks[:])
    _flag = False
    for _line in _vcf_file:
        _toks = _line.strip().split('\t')
        if _toks[0] != _chrom:
            break
        _i += 1
        if _i >= _start and _i < _stop:
            _ret_sites.append(_toks[:])
            continue
        if _i >= _stop:
            _flag = True
            break
        # print "### WARNING (sample by count) : Something dubious is happening when reading the desired interval: chrom={}, start={}, stop={}, last evaluated position is: \n\t{}".format(_chrom, _start, _stop, '\t'.join(_toks[0:8]))
    if not _flag:
        print "### WARNING (sample by count) : Site selection stopped for another reason than reaching the end of the interval: chrom={}, start={}, stop={}, last evaluated position is: \n\t{}".format(_chrom, _start, _stop, '\t'.join(_toks[0:8]))
    return _ret_sites[:]
###
# the sequential ifs in this method are not (all) mutually exclusive -- so order matters!!!
def compare_genotypes(_a, _b, _null_allele):
    # _ret = [total_nc, total_c, identical_nc, identical_c]
    # print '{} {} -> {}'.format(_a,_b,_null_allele)
    _ret = [0.0,0.0,0.0,0.0]
    if _a[0] == _null_allele and _a[1] == _null_allele or _b[0] == _null_allele and _b[1] == _null_allele :
        return _ret
    
    if (_a[0] == _b[0] and _a[1] == _b[1]) or (_a[1] == _b[0] and _a[0] == _b[1]):
        if _a[0] == _a[1]:
            _ret = [1.0,1.0,1.0,1.0]
        else:
            if _a[0] == _null_allele or _a[1] == _null_allele:
                print "### DBG: gt counts should not be equal"
                _ret = [1.0,0.0,0.5,0.0]
            else:
                _ret = [1.0,1.0,0.5,0.5]
        return _ret
    
    
    if _a[0] != _null_allele and _a[1] != _null_allele and _b[0] != _null_allele and _b[1] != _null_allele :
        if _a[0] == _b[0] or _a[0] == _b[1] or _a[1] == _b[0] or _a[1] == _b[1]:
            return [1.0,1.0,0.5,0.5]
        return [1.0,1.0,0.0,0.0]
    
    if _a[0] != _null_allele and (_a[0] == _b[0] or _a[0] == _b[1]) or _a[1] != _null_allele and (_a[1] == _b[0] or  _a[1] == _b[1]):
        print "### DBG: gt counts should not be equal 2"
        return [1.0,0.0,0.5,0.0]
    
    if (_a[0] != _null_allele or _a[1] != _null_allele) and (_b[0] != _null_allele or _b[1] != _null_allele) : # this if is useless
        print "### DBG: Something is wrong, reached the useless if-statement"
        return [1.0,0.0,0.0,0.0]
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
def get_index_of_tag(_format, _tag):
    _toks = _format.split(':')
    _idx_list = [ _i for _i in range(len(_toks)) if _toks[_i] == _tag ]
    if len(_idx_list) == 0:
        return -1
    if len(_idx_list) > 1:
        print '### WARNING: found some tag multiple times in VCF format field (returning first occurance): {} : {}'.format(_tag, _idx_list)
    return _idx_list[0]
###
def get_symbolic_GT_as_list(_data, _idx):
    _toks = _data.split(':')
    if len(_toks) <= _idx:
        print '### ERROR: alleged index of GT tag does not match the data: idx={}; data={}'.format(_idx, _data)
        exit(1)
    _gt = _toks[_idx]
    if len(_gt) == len(_gt.translate(None, '|/')):
        print '### ERROR: genotype not well formed: {}'.format(_gt)
        exit(1)
    return _gt.translate(maketrans('/|', '  '),'').split(' ')
###

####################################
### void main(int argc, char **argv)


vcf_file_name = sys.argv[1]
out_file_name = sys.argv[2]
ped_file_name = sys.argv[3]
chrom = sys.argv[4]
start = int(sys.argv[5])
stop  = int(sys.argv[6])
hash_chr_start = int(sys.argv[7])
hash_chr_stop  = int(sys.argv[8])
sample_type = int(sys.argv[9])


vcf_file = open(vcf_file_name, 'r')
safe_drop_header(vcf_file)
vcf_header = vcf_file.readline().strip().split('\t')
vcf_header = map (lambda x : x.lower(), vcf_header)

# build the set of couples to evaluate, from the intersection of the ped file and the found VCF samples
ped_file = open(ped_file_name, 'r')
trios = [line.strip().split('\t') for line in ped_file if line.strip().split('\t')[1].lower()[-1] == 'c' ] # or line.strip().split('\t')[1].lower()[-1] == 'd']
trios = map(lambda x : map(lambda y : y.lower(),x), trios)
trios, indiv_mapping = map_ped_trios_to_vcf(trios, vcf_header)

# extract all the region of interest from the VCF
if sample_type == 1:
    sites = extract_sites_subset_length(vcf_file, chrom, start, stop, hash_chr_start, hash_chr_stop)
elif sample_type == 2:
    sites = extract_sites_subset_count(vcf_file, chrom, start, stop, hash_chr_start, hash_chr_stop)


# count identical genotype proportions
interval_stats = [0.0] * len(trios)
gt_tag_idx = get_index_of_tag(vcf_header[8],'GT')
for s in range(len(sites)):
    site_sums = []
    for ind in range(len(trios)):
        _gt_1 = get_symbolic_GT_as_list(sites[s][indiv_mapping[trios[ind][2]]], gt_tag_idx)
        _gt_2 = get_symbolic_GT_as_list(sites[s][indiv_mapping[trios[ind][3]]], gt_tag_idx)
        _site_couple_stats = compare_genotypes(_gt_1, _gt_2, '.')
        site_sums.append(_site_couple_stats[-1])
    interval_stats = [ interval_stats[ii] + site_sums[ii] for ii in range(len(interval_stats)) ]

Nsites = len(sites)

interval_stats = [ interval_stats[ii]/float(Nsites) for ii in range(0,len(interval_stats)) ]

# write output 
out_list = [sites[0][0], sites[0][1], sites[-1][1], len(sites)]
out_list.extend(interval_stats)

out_file = open(out_file_name, 'w')
out_file.write( '\t'.join(map(str, out_list)) + '\n' )
out_file.close()





























