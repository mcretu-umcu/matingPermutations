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
  hm.permuted.averages.var = as.data.frame( matrix( , 0, 4 + length(columns_to_permute), dimnames = list( c(), c(columns_to_permute, "colour", "data", "region", "freq") )) )
  
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
    hm.permuted.averages.var = add_average_of_permuted_columns(hm.permuted.averages.var, hm.ceu.snps.mhc.005, columns_to_permute.var, recoupledG.var, c("rand", as.character(args[3]), "MHC", "maf0.5"))
    ## nor filter SNPs
    # hm.permuted.averages.var = add_average_of_permuted_columns(hm.permuted.averages.var, hm.yri.snps.mhc.005, columns_to_permute.var, recoupledG.var, c("rand", "HapMapYRI", "MHC", "maf0.5"))

  }
  return (hm.permuted.averages.var)
}
#######
args = commandArgs(trailingOnly=TRUE)

maf = "0.01"
hardPath = '/hpc/cog_bioinf/kloosterman/users/mcretu/mating/distances/hapmap/'
hm.ceu.snps.mhc.005 = read.table(paste(hardPath, as.character(args[1]), sep = "", collapse=NULL), header = T, sep = '\t')
# hm.yri.snps.mhc.005 = read.table(paste(hardPath, 'distances.yri.mhc.005.txt', sep = "", collapse=NULL), header = T, sep = '\t')



ceu.kids = read.table(as.character(args[2]), header = F, sep = '\t')
# yri.kids = read.table('/hpc/cog_bioinf/kloosterman/users/mcretu/mating/res/yri.kids.subset.ped', header = F, sep = '\t')
# eliminate twins; i.e.: make sure each pair of parents appears only once

ceu.kids$V1 = as.character(ceu.kids$V1)
ceu.kids$V3 = as.character(ceu.kids$V3)
ceu.kids$V4 = as.character(ceu.kids$V4)
# yri.kids$V1 = as.character(yri.kids$V1)
# yri.kids$V3 = as.character(yri.kids$V3)
# yri.kids$V4 = as.character(yri.kids$V4)


### colour mates -- non-mates
hm.ceu.snps.mhc.005 = colour_couples_from_ped(ceu.kids, hm.ceu.snps.mhc.005, "mates", "rand")
# hm.yri.snps.mhc.005 = colour_couples_from_ped(yri.kids, hm.yri.snps.mhc.005, "mates", "rand")




hm.ceu.snps.mhc.005$data = rep("MHC", length(hm.ceu.snps.mhc.005[,1]))
# hm.yri.snps.mhc.005$data = rep("MHC", length(hm.yri.snps.mhc.005[,1]))

##
hm.ceu.snps.mhc.005$freq = rep("maf0.5", length(hm.ceu.snps.mhc.005[,1]))
# hm.yri.snps.mhc.005$freq = rep("maf0.5", length(hm.yri.snps.mhc.005[,1]))

### typset columns to be sure nothing gets treated as factors
hm.ceu.snps.mhc.005 = typeset_cols_PCs(hm.ceu.snps.mhc.005, c())
# hm.yri.snps.mhc.005 = typeset_cols_PCs(hm.yri.snps.mhc.005, c())

######## Permute for p-value
NR_PERM_G = 1000
columns_to_permute = c("qCgt", "qCal")
couples_prefix = "coupleAvg"
pval_prefix = "pVal"

columns_to_plot = columns_to_permute

# cluster = makeCluster(25,type="SOCK")
# registerDoSNOW(cluster)
hm.permuted.averages = as.data.frame( matrix( , 0, 4 + length(columns_to_permute), dimnames = list( c(), c(columns_to_permute, "colour", "data", "region", "freq") )) )

# results = foreach (i = 1:NR_PERM_G, .combine=rbind) %do% perform_one_permutation(gonl.wg.mhc, gonl.wg.mhcX, gonl.wg.wg, gonl.hm.mhc, gonl.hm.mhcX, gonl.hm.wg, columns_to_permute, gonl.kids,1000)
results = perform_one_permutation( columns_to_permute, ceu.kids,1000)

colnames (results) = colnames(hm.permuted.averages)
write.table(results, file=paste('/hpc/cog_bioinf/kloosterman/users/mcretu/mating/permutations/hapmap/',as.character(args[3]), '/perm.hm.full.1M.',as.character(args[4]),'.txt',sep='',collapse=NULL), sep='\t', quote=FALSE, row.names=F, col.names=F)
# stopCluster(cluster)



