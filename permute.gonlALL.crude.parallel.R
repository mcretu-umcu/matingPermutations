library(doSNOW)
library(foreach)

###
colour_couples_from_ped <- function (ped, data, colour,colour_base)
{
  temp1 = paste(data$ind1,data$ind2, sep=" ", collapse = NULL)
  temp2 = paste(data$ind2,data$ind1, sep=" ", collapse = NULL)
  couples = paste(ped[,3], ped[,4], sep=" ", collapse=NULL)
  
  data$colour = rep(colour_base,length(data[,1]))
  data$colour[temp1 %in% couples] = colour
  # print (temp1)
  # print (temp1 %in% couples)
  # print (couples)
  data$colour[temp2 %in% couples] = colour
  
  return(data)
}
###
# Because R stands for "R"etarded!
typeset_cols_PCs <- function(dataset, pcs_used)
{
  dataset$ind1 = as.character(dataset$ind1)
  dataset$ind2 = as.character(dataset$ind2)
  dataset$fam1 = as.character(dataset$fam1)
  dataset$fam2 = as.character(dataset$fam2)
  if ( match("colour", colnames(dataset), nomatch = -1) != -1 )
  {
    dataset$colour = as.character(dataset$colour)
  }
  dataset$qCgt = as.numeric(as.character(dataset$qCgt))
  dataset$qCal = as.numeric(as.character(dataset$qCal))
  
  if (length(pcs_used) == 0)
  {
    return (dataset)
  }
  for (i in pcs_used)
  {
    dataset[,paste("i1",i, sep="", collapse=NULL)] = as.numeric(as.character(dataset[,paste("i1",i, sep="", collapse=NULL)]))
    dataset[,paste("i2",i, sep="", collapse=NULL)] = as.numeric(as.character(dataset[,paste("i2",i, sep="", collapse=NULL)]))
    dataset[,paste("couple",i,"avg", sep="", collapse=NULL)] = as.numeric(as.character(dataset[,paste("couple",i,"avg", sep="", collapse=NULL)]))
    dataset[,paste("couple",i,"se", sep="", collapse=NULL)] = as.numeric(as.character(dataset[,paste("couple",i,"se", sep="", collapse=NULL)]))
  }
  
  return (dataset)
}
###
sample_couples2 <- function (males, females, true_pairs, NallowedPairs = 1)
{
  males = as.character(males)
  females = as.character(females)
  true_pairs$m = as.character(true_pairs$m)
  true_pairs$f = as.character(true_pairs$f)
  
  natural = 1:length(males)
  # randomized = sample(1:length(males), length(males))
  # test = natural - randomized
  # print (length(test[test == 0 ]))
  # new_females = females[randomized]
  # return ( data.frame(m=males, f=new_females) )
  # }
  randomized = c()
  while (T)
  {
    randomized = sample(1:length(males), length(males))
    test = natural - randomized
    if ( length(test[test == 0]) <= NallowedPairs )
      break
  }
  new_females = females[randomized]
  # print (length(test[test == 0]))
  return ( data.frame(m=males, f=new_females) )
}
###
add_average_of_permuted_columns <- function (permutationFrame, allData, columns, recoupleG, rest)
{
  
  
  subsetG = allData[ match( as.character(allData$ind1), as.character(recoupleG$m) ) ==  match( as.character(allData$ind2), as.character(recoupleG$f) ) | match( as.character(allData$ind2), as.character(recoupleG$m) ) ==  match( as.character(allData$ind1), as.character(recoupleG$f) ) , ]
  missedCouples = length(which(is.na(subsetG)))
  # print (length(subsetG[,1]))
  subsetG = na.omit(subsetG)
  # print (length(subsetG[,1]))
  # print ("###############")
  
  means = c()
  for (i in 1:length(columns))
  {
    # print (columns[i])
    subsetG[,columns[i]] = as.numeric(as.character(subsetG[,columns[i]]))
    means = c( means, mean(subsetG[,columns[i]]) )
    # print (subsetG[,columns[i]])
    # print (mean(subsetG[,columns[i]]))
  }
  one.point = c( means, rest )
  to.add = t(data.frame(one.point))
  colnames(to.add) = colnames(permutationFrame)
  permutationFrame = rbind(permutationFrame, to.add)
  
  return (permutationFrame)
}
###
perform_one_permutation = function ( columns_to_permute.var, gonl.kids.var,ii )
{
  gonl.permuted.averages.var = as.data.frame( matrix( , 0, 4 + length(columns_to_permute), dimnames = list( c(), c(columns_to_permute, "colour", "data", "region", "freq") )) )
  
  for (iii in 1:ii)
  {
    if (iii %% 100 == 0)
    {print (iii)}
    
    recoupledG.var = sample_couples2( as.character(gonl.kids.var$V3), as.character(gonl.kids.var$V4), data.frame(m=as.character(gonl.kids.var$V3), f=as.character(gonl.kids.var$V4)), 2 )
    recoupledG.var$m = as.character(recoupledG.var$m)
    recoupledG.var$f = as.character(recoupledG.var$f)
    # print (recoupledG.var)
    # print (length(recoupledG.var$m))
  
  
  
    ### ALL 
    gonl.permuted.averages.var = add_average_of_permuted_columns(gonl.permuted.averages.var, gonl.wg.all.mhc.000, columns_to_permute.var, recoupledG.var, c("rand", "wholeGenomeAll", "MHC", "maf0"))
    gonl.permuted.averages.var = add_average_of_permuted_columns(gonl.permuted.averages.var, gonl.wg.all.mhc.005, columns_to_permute.var, recoupledG.var, c("rand", "wholeGenomeAll", "MHC", "maf0.5"))
    gonl.permuted.averages.var = add_average_of_permuted_columns(gonl.permuted.averages.var, gonl.wg.all.mhcX.005, columns_to_permute.var, recoupledG.var, c("rand", "wholeGenomeAll", "xMHC", "maf0.5"))
    gonl.permuted.averages.var = add_average_of_permuted_columns(gonl.permuted.averages.var, gonl.wg.all.wg.005, columns_to_permute.var, recoupledG.var, c("rand", "wholeGenomeAll", "GenomeWide", "maf0.5"))
    ## nor filter SNPs
    gonl.permuted.averages.var = add_average_of_permuted_columns(gonl.permuted.averages.var, gonl.wg.snps.mhc.000, columns_to_permute.var, recoupledG.var, c("rand", "wholeGenomeSNPs", "MHC", "maf0"))
    gonl.permuted.averages.var = add_average_of_permuted_columns(gonl.permuted.averages.var, gonl.wg.snps.mhcX.000, columns_to_permute.var, recoupledG.var, c("rand", "wholeGenomeSNPs", "xMHC", "maf0"))
    gonl.permuted.averages.var = add_average_of_permuted_columns(gonl.permuted.averages.var, gonl.wg.snps.wg.000, columns_to_permute.var, recoupledG.var, c("rand", "wholeGenomeSNPs", "GenomeWide", "maf0"))
    ## rare
    gonl.permuted.averages.var = add_average_of_permuted_columns(gonl.permuted.averages.var, gonl.wg.snps.mhc.005, columns_to_permute.var, recoupledG.var, c("rand", "wholeGenomeSNPs", "MHC", "maf0.5"))
    gonl.permuted.averages.var = add_average_of_permuted_columns(gonl.permuted.averages.var, gonl.wg.snps.mhcX.005, columns_to_permute.var, recoupledG.var, c("rand", "wholeGenomeSNPs", "xMHC", "maf0.5"))    
    gonl.permuted.averages.var = add_average_of_permuted_columns(gonl.permuted.averages.var, gonl.wg.snps.wg.005, columns_to_permute.var, recoupledG.var, c("rand", "wholeGenomeSNPs", "GenomeWide", "maf0.5"))
    ## low and common
    gonl.permuted.averages.var = add_average_of_permuted_columns(gonl.permuted.averages.var, gonl.wg.snps.mhc.050, columns_to_permute.var, recoupledG.var, c("rand", "wholeGenomeSNPs", "MHC", "maf5"))
    gonl.permuted.averages.var = add_average_of_permuted_columns(gonl.permuted.averages.var, gonl.wg.snps.wg.050, columns_to_permute.var, recoupledG.var, c("rand", "wholeGenomeSNPs", "GenomeWide", "maf5"))
    gonl.permuted.averages.var = add_average_of_permuted_columns(gonl.permuted.averages.var, gonl.wg.snps.mhc.200, columns_to_permute.var, recoupledG.var, c("rand", "wholeGenomeSNPs", "MHC", "maf20"))
    gonl.permuted.averages.var = add_average_of_permuted_columns(gonl.permuted.averages.var, gonl.wg.snps.wg.200, columns_to_permute.var, recoupledG.var, c("rand", "wholeGenomeSNPs", "GenomeWide", "maf20"))

  }
  return (gonl.permuted.averages.var)
}
#######
args = commandArgs(trailingOnly=TRUE)

maf = "0.01"
hardPath = '/hpc/cog_bioinf/kloosterman/users/mcretu/mating/distances/gonl/'
gonl.wg.all.mhc.000 = read.table(paste(hardPath, 'wg/gonl.all.wg.maf0.000.mhc.allele.txt', sep = "", collapse=NULL), header = T, sep = '\t')
gonl.wg.all.mhc.005 = read.table(paste(hardPath, 'wg/gonl.all.wg.maf0.005.mhc.allele.txt', sep = "", collapse=NULL), header = T, sep = '\t')
gonl.wg.all.mhcX.005 = read.table(paste(hardPath, 'wg/gonl.all.wg.maf0.005.mhcX.allele.txt', sep = "", collapse=NULL), header = T, sep = '\t')
gonl.wg.all.wg.005 = read.table(paste(hardPath, 'wg/gonl.all.wg.maf0.005.wg.allele.txt', sep = "", collapse=NULL), header = T, sep = '\t')

gonl.wg.snps.mhc.000 = read.table(paste(hardPath, 'wg/gonl.snps.wg.maf0.000.mhc.allele.txt', sep = "", collapse=NULL), header = T, sep = '\t')
gonl.wg.snps.mhcX.000 = read.table(paste(hardPath, 'wg/gonl.snps.wg.maf0.000.mhcX.allele.txt', sep = "", collapse=NULL), header = T, sep = '\t')
gonl.wg.snps.wg.000 = read.table(paste(hardPath, 'wg/gonl.snps.wg.maf0.000.wg.allele.txt', sep = "", collapse=NULL), header = T, sep = '\t')

gonl.wg.snps.mhc.005 = read.table(paste(hardPath, 'wg/gonl.snps.wg.maf0.005.mhc.allele.txt', sep = "", collapse=NULL), header = T, sep = '\t')
gonl.wg.snps.mhcX.005 = read.table(paste(hardPath, 'wg/gonl.snps.wg.maf0.005.mhcX.allele.txt', sep = "", collapse=NULL), header = T, sep = '\t')
gonl.wg.snps.wg.005 = read.table(paste(hardPath, 'wg/gonl.snps.wg.maf0.005.wg.allele.txt', sep = "", collapse=NULL), header = T, sep = '\t')

gonl.wg.snps.mhc.050 = read.table(paste(hardPath, 'wg/gonl.snps.wg.maf0.050.mhc.allele.txt', sep = "", collapse=NULL), header = T, sep = '\t')
gonl.wg.snps.wg.050 = read.table(paste(hardPath, 'wg/gonl.snps.wg.maf0.050.wg.allele.txt', sep = "", collapse=NULL), header = T, sep = '\t')
gonl.wg.snps.mhc.200 = read.table(paste(hardPath, 'wg/gonl.snps.wg.maf0.200.mhc.allele.txt', sep = "", collapse=NULL), header = T, sep = '\t')
gonl.wg.snps.wg.200 = read.table(paste(hardPath, 'wg/gonl.snps.wg.maf0.200.wg.allele.txt', sep = "", collapse=NULL), header = T, sep = '\t')



gonl.kids = read.table('/hpc/cog_bioinf/kloosterman/users/mcretu/mating/res/gonl.public.kids.ped', header = F, sep = '\t')
# eliminate twins; i.e.: make sure each pair of parents appears only once
gonl.kids = gonl.kids[substr(as.character(gonl.kids$V2), nchar(as.character(gonl.kids$V2)), nchar(as.character(gonl.kids$V2))) != "d",]
related_couples = c("gonl-127", "gonl-9", "gonl-173", "gonl-186", "gonl-18", "gonl-81", "gonl-83", "gonl-120", "gonl-199")
# related_couples = c('gonl-245','gonl-244','gonl-247','gonl-246','gonl-128','gonl-240','gonl-243','gonl-242','gonl-125','gonl-126','gonl-127','gonl-120','gonl-121','gonl-122','gonl-123','gonl-60','gonl-63','gonl-64','gonl-66','gonl-67','gonl-68','gonl-69','gonl-139','gonl-137','gonl-134','gonl-132','gonl-131','gonl-130','gonl-91','gonl-232','gonl-233','gonl-234','gonl-235','gonl-237','gonl-239','gonl-59','gonl-213','gonl-241','gonl-56','gonl-51','gonl-53','gonl-39','gonl-38','gonl-36','gonl-35','gonl-34','gonl-33','gonl-102','gonl-103','gonl-100','gonl-106','gonl-104','gonl-200','gonl-222','gonl-188','gonl-189','gonl-226','gonl-225','gonl-224','gonl-182','gonl-183','gonl-229','gonl-181','gonl-186','gonl-187','gonl-184','gonl-248','gonl-49','gonl-41','gonl-46','gonl-47','gonl-249','gonl-29','gonl-25','gonl-26','gonl-27','gonl-20','gonl-21','gonl-22','gonl-23','gonl-115','gonl-114','gonl-116','gonl-111','gonl-110','gonl-112','gonl-216','gonl-217','gonl-214','gonl-215','gonl-9','gonl-198','gonl-211','gonl-195','gonl-194','gonl-197','gonl-196','gonl-191','gonl-190','gonl-193','gonl-192','gonl-11','gonl-16','gonl-209','gonl-90','gonl-95','gonl-94','gonl-97','gonl-99','gonl-98','gonl-207','gonl-206','gonl-160','gonl-163','gonl-164','gonl-165','gonl-166','gonl-167','gonl-168','gonl-169','gonl-208','gonl-86','gonl-212','gonl-84','gonl-83','gonl-80','gonl-81','gonl-173','gonl-171','gonl-176','gonl-175','gonl-55','gonl-178','gonl-7','gonl-105','gonl-218','gonl-145','gonl-142','gonl-141','gonl-148','gonl-149','gonl-221','gonl-158','gonl-250','gonl-151','gonl-153','gonl-152','gonl-154','gonl-157','gonl-156','gonl-72','gonl-71','gonl-70','gonl-77','gonl-76','gonl-75','gonl-79')

gonl.kids$V1 = as.character(gonl.kids$V1)
gonl.kids$V3 = as.character(gonl.kids$V3)
gonl.kids$V4 = as.character(gonl.kids$V4)
gonl.kids = gonl.kids[ match(gonl.kids$V1, related_couples, nomatch=-1) == -1, ]

### colour mates -- non-mates
gonl.wg.all.mhc.000 = colour_couples_from_ped(gonl.kids, gonl.wg.all.mhc.000, "mates", "rand")
gonl.wg.all.mhc.005 = colour_couples_from_ped(gonl.kids, gonl.wg.all.mhc.005, "mates", "rand")
gonl.wg.all.mhcX.005 = colour_couples_from_ped(gonl.kids, gonl.wg.all.mhcX.005, "mates", "rand")
gonl.wg.all.wg.005 = colour_couples_from_ped(gonl.kids, gonl.wg.all.wg.005, "mates", "rand")

gonl.wg.snps.mhc.000 = colour_couples_from_ped(gonl.kids, gonl.wg.snps.mhc.000, "mates", "rand")
gonl.wg.snps.mhcX.000 = colour_couples_from_ped(gonl.kids, gonl.wg.snps.mhcX.000, "mates", "rand")
gonl.wg.snps.wg.000 = colour_couples_from_ped(gonl.kids, gonl.wg.snps.wg.000, "mates", "rand")
gonl.wg.snps.mhc.005 = colour_couples_from_ped(gonl.kids, gonl.wg.snps.mhc.005, "mates", "rand")
gonl.wg.snps.mhcX.005 = colour_couples_from_ped(gonl.kids, gonl.wg.snps.mhcX.005, "mates", "rand")
gonl.wg.snps.wg.005 = colour_couples_from_ped(gonl.kids, gonl.wg.snps.wg.005, "mates", "rand")

gonl.wg.snps.mhc.050 = colour_couples_from_ped(gonl.kids, gonl.wg.snps.mhc.050, "mates", "rand")
gonl.wg.snps.wg.050 = colour_couples_from_ped(gonl.kids, gonl.wg.snps.wg.050, "mates", "rand")
gonl.wg.snps.mhc.200 = colour_couples_from_ped(gonl.kids, gonl.wg.snps.mhc.200, "mates", "rand")
gonl.wg.snps.wg.200 = colour_couples_from_ped(gonl.kids, gonl.wg.snps.wg.200, "mates", "rand")

# View(gonl.mhc1[gonl.mhc1$colour == "mates",])

################## Don't use anymore, for the moment
# gonl.wg.mhc$qCgtDif = gonl.wg.wg$qCgt[match(paste(gonl.wg.mhc$ind1,gonl.wg.mhc$ind2, sep=" ", collapse = NULL), paste(gonl.wg.wg$ind1,gonl.wg.wg$ind2, sep=" ", collapse = NULL))] - gonl.wg.mhc$qCgt
# gonl.wg.mhcX$qCgtDif = gonl.wg.wg$qCgt[match(paste(gonl.wg.mhcX$ind1,gonl.wg.mhcX$ind2, sep=" ", collapse = NULL), paste(gonl.wg.wg$ind1,gonl.wg.wg$ind2, sep=" ", collapse = NULL))] - gonl.wg.mhcX$qCgt
# gonl.hm.mhc$qCgtDif = gonl.hm.wg$qCgt[match(paste(gonl.hm.mhc$ind1,gonl.hm.mhc$ind2, sep=" ", collapse = NULL), paste(gonl.hm.wg$ind1,gonl.hm.wg$ind2, sep=" ", collapse = NULL))] - gonl.hm.mhc$qCgt
# gonl.hm.mhcX$qCgtDif = gonl.hm.wg$qCgt[match(paste(gonl.hm.mhcX$ind1,gonl.hm.mhcX$ind2, sep=" ", collapse = NULL), paste(gonl.hm.wg$ind1,gonl.hm.wg$ind2, sep=" ", collapse = NULL))] - gonl.hm.mhcX$qCgt
# gonl.wg.wg$qCgtDif = gonl.wg.wg$qCgt
# gonl.hm.wg$qCgtDif = gonl.hm.wg$qCgt
###################
# pl = ggplot()
# pl = pl + geom_histogram(aes(gonl.wg.mhc$qCgt[gonl.wg.mhc$colour == "mates"]),colour="blue")
# pl = pl + geom_histogram(aes(gonl.wg.mhc$qCgt[gonl.wg.mhc$colour == "rand"]),colour="green")
# pl


gonl.wg.all.mhc.000$data = rep("MHC", length(gonl.wg.all.mhc.000[,1]))
gonl.wg.all.mhc.005$data = rep("MHC", length(gonl.wg.all.mhc.005[,1]))
gonl.wg.all.mhcX.005$data = rep("xMHC", length(gonl.wg.all.mhcX.005[,1]))
gonl.wg.all.wg.005$data = rep("GenomeWide", length(gonl.wg.all.wg.005[,1]))

gonl.wg.snps.mhc.000$data = rep("MHC", length(gonl.wg.snps.mhc.000[,1]))
gonl.wg.snps.mhcX.000$data = rep("xMHC", length(gonl.wg.snps.mhcX.000[,1]))
gonl.wg.snps.wg.000$data = rep("GenomeWide", length(gonl.wg.snps.wg.000[,1]))
gonl.wg.snps.mhc.005$data = rep("MHC", length(gonl.wg.snps.mhc.005[,1]))
gonl.wg.snps.mhcX.005$data = rep("xMHC", length(gonl.wg.snps.mhcX.005[,1]))
gonl.wg.snps.wg.005$data = rep("GenomeWide", length(gonl.wg.snps.wg.005[,1]))

gonl.wg.snps.mhc.050$data = rep("MHC", length(gonl.wg.snps.mhc.050[,1]))
gonl.wg.snps.wg.050$data = rep("GenomeWide", length(gonl.wg.snps.wg.050[,1]))
gonl.wg.snps.mhc.200$data = rep("MHC", length(gonl.wg.snps.mhc.200[,1]))
gonl.wg.snps.wg.200$data = rep("GenomeWide", length(gonl.wg.snps.wg.200[,1]))
##
gonl.wg.all.mhc.000$freq = rep("maf0", length(gonl.wg.all.mhc.000[,1]))
gonl.wg.all.mhc.005$freq = rep("maf0.5", length(gonl.wg.all.mhc.005[,1]))
gonl.wg.all.mhcX.005$freq = rep("maf0.5", length(gonl.wg.all.mhcX.005[,1]))
gonl.wg.all.wg.005$freq = rep("maf0.5", length(gonl.wg.all.wg.005[,1]))

gonl.wg.snps.mhc.000$freq = rep("maf0", length(gonl.wg.snps.mhc.000[,1]))
gonl.wg.snps.mhcX.000$freq = rep("maf0", length(gonl.wg.snps.mhcX.000[,1]))
gonl.wg.snps.wg.000$freq = rep("maf0", length(gonl.wg.snps.wg.000[,1]))
gonl.wg.snps.mhc.005$freq = rep("maf0.5", length(gonl.wg.snps.mhc.005[,1]))
gonl.wg.snps.mhcX.005$freq = rep("maf0.5", length(gonl.wg.snps.mhcX.005[,1]))
gonl.wg.snps.wg.005$freq = rep("maf0.5", length(gonl.wg.snps.wg.005[,1]))

gonl.wg.snps.mhc.050$freq = rep("maf5", length(gonl.wg.snps.mhc.050[,1]))
gonl.wg.snps.wg.050$freq = rep("maf5", length(gonl.wg.snps.wg.050[,1]))
gonl.wg.snps.mhc.200$freq = rep("maf20", length(gonl.wg.snps.mhc.200[,1]))
gonl.wg.snps.wg.200$freq = rep("maf20", length(gonl.wg.snps.wg.200[,1]))


### typset columns to be sure nothing gets treated as factors
gonl.wg.all.mhc.000 = typeset_cols_PCs(gonl.wg.all.mhc.000, c())
gonl.wg.all.mhc.005 = typeset_cols_PCs(gonl.wg.all.mhc.005, c())
gonl.wg.all.mhcX.005 = typeset_cols_PCs(gonl.wg.all.mhcX.005, c())
gonl.wg.all.wg.005 = typeset_cols_PCs(gonl.wg.all.wg.005, c())

gonl.wg.snps.mhc.000 = typeset_cols_PCs(gonl.wg.snps.mhc.000, c())
gonl.wg.snps.mhcX.000 = typeset_cols_PCs(gonl.wg.snps.mhcX.000, c())
gonl.wg.snps.wg.000 = typeset_cols_PCs(gonl.wg.snps.wg.000, c())
gonl.wg.snps.mhc.005 = typeset_cols_PCs(gonl.wg.snps.mhc.005, c())
gonl.wg.snps.mhcX.005 = typeset_cols_PCs(gonl.wg.snps.mhcX.005, c())
gonl.wg.snps.wg.005 = typeset_cols_PCs(gonl.wg.snps.wg.005, c())

gonl.wg.snps.mhc.050 = typeset_cols_PCs(gonl.wg.snps.mhc.050, c())
gonl.wg.snps.wg.050 = typeset_cols_PCs(gonl.wg.snps.wg.050, c())
gonl.wg.snps.mhc.200 = typeset_cols_PCs(gonl.wg.snps.mhc.200, c())
gonl.wg.snps.wg.200 = typeset_cols_PCs(gonl.wg.snps.wg.200, c())

######## Permute for p-value
NR_PERM_G = 1000
columns_to_permute = c("qCgt", "qCal")
couples_prefix = "coupleAvg"
pval_prefix = "pVal"

columns_to_plot = columns_to_permute

# cluster = makeCluster(25,type="SOCK")
# registerDoSNOW(cluster)
gonl.permuted.averages = as.data.frame( matrix( , 0, 4 + length(columns_to_permute), dimnames = list( c(), c(columns_to_permute, "colour", "data", "region", "freq") )) )

# results = foreach (i = 1:NR_PERM_G, .combine=rbind) %do% perform_one_permutation(gonl.wg.mhc, gonl.wg.mhcX, gonl.wg.wg, gonl.hm.mhc, gonl.hm.mhcX, gonl.hm.wg, columns_to_permute, gonl.kids,1000)
results = perform_one_permutation( columns_to_permute, gonl.kids,1000)

colnames (results) = colnames(gonl.permuted.averages)
write.table(results, file=paste('/hpc/cog_bioinf/kloosterman/users/mcretu/mating/permutations/wg/perm.gonl.full.1M.',as.character(args[1]),'.txt',sep='',collapse=NULL), sep='\t', quote=FALSE, row.names=F, col.names=F)
# stopCluster(cluster)



