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
  dataset$corr = as.numeric(as.character(dataset$corr))
  # dataset$qCal = as.numeric(as.character(dataset$qCal))
  
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
  
  
  
    ### ALL 
    gonl.permuted.averages.var = add_average_of_permuted_columns(gonl.permuted.averages.var, gonl.wg.all.mhc.imputed, columns_to_permute.var, recoupledG.var, c("rand", "imputedAll", "MHC", "maf0"))
    gonl.permuted.averages.var = add_average_of_permuted_columns(gonl.permuted.averages.var, gonl.wg.ohm.mhc.imputed, columns_to_permute.var, recoupledG.var, c("rand", "imputedOHM", "MHC", "maf0"))

  }
  return (gonl.permuted.averages.var)
}
#######
args = commandArgs(trailingOnly=TRUE)

maf = "0.01"
hardPath = '/hpc/cog_bioinf/kloosterman/users/mcretu/mating/distances/gonl/'
gonl.wg.all.mhc.imputed = read.table(paste(hardPath, 'gonl.all.snp2hla.txt', sep = "", collapse=NULL), header = T, sep = '\t')
gonl.wg.ohm.mhc.imputed = read.table(paste(hardPath, 'gonl.all.snp2hla.OHM.txt', sep = "", collapse=NULL), header = T, sep = '\t')


gonl.kids = read.table('/hpc/cog_bioinf/kloosterman/users/mcretu/mating/res/gonl.public.kids.ped', header = F, sep = '\t')
related_couples = c("gonl-127", "gonl-9", "gonl-173", "gonl-186", "gonl-18", "gonl-81", "gonl-83", "gonl-120", "gonl-199")
# eliminate twins; i.e.: make sure each pair of parents appears only once



gonl.kids$V1 = as.character(gonl.kids$V1)
gonl.kids$V3 = as.character(gonl.kids$V3)
gonl.kids$V4 = as.character(gonl.kids$V4)
gonl.kids = gonl.kids[ match(gonl.kids$V1, related_couples, nomatch=-1) == -1, ]
gonl.kids = gonl.kids[substr(as.character(gonl.kids$V2), nchar(as.character(gonl.kids$V2)), nchar(as.character(gonl.kids$V2))) != "d",]

### colour mates -- non-mates
gonl.wg.all.mhc.imputed = colour_couples_from_ped(gonl.kids, gonl.wg.all.mhc.imputed, "mates", "rand")
gonl.wg.ohm.mhc.imputed = colour_couples_from_ped(gonl.kids, gonl.wg.ohm.mhc.imputed, "mates", "rand")

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


gonl.wg.all.mhc.imputed$data = rep("MHC", length(gonl.wg.all.mhc.imputed[,1]))
gonl.wg.all.mhc.imputed$freq = rep("maf0", length(gonl.wg.all.mhc.imputed[,1]))

gonl.wg.ohm.mhc.imputed$data = rep("MHC", length(gonl.wg.ohm.mhc.imputed[,1]))
gonl.wg.ohm.mhc.imputed$freq = rep("maf0", length(gonl.wg.ohm.mhc.imputed[,1]))



### typset columns to be sure nothing gets treated as factors
gonl.wg.all.mhc.imputed = typeset_cols_PCs(gonl.wg.all.mhc.imputed, c())
gonl.wg.ohm.mhc.imputed = typeset_cols_PCs(gonl.wg.ohm.mhc.imputed, c())

######## Permute for p-value
NR_PERM_G = 1000
columns_to_permute = c("corr")
couples_prefix = "coupleAvg"
pval_prefix = "pVal"

columns_to_plot = columns_to_permute

# cluster = makeCluster(25,type="SOCK")
# registerDoSNOW(cluster)
gonl.permuted.averages = as.data.frame( matrix( , 0, 4 + length(columns_to_permute), dimnames = list( c(), c(columns_to_permute, "colour", "data", "region", "freq") )) )

# results = foreach (i = 1:NR_PERM_G, .combine=rbind) %do% perform_one_permutation(gonl.wg.mhc, gonl.wg.mhcX, gonl.wg.wg, gonl.hm.mhc, gonl.hm.mhcX, gonl.hm.wg, columns_to_permute, gonl.kids,1000)
results = perform_one_permutation( columns_to_permute, gonl.kids, 1000)

colnames (results) = colnames(gonl.permuted.averages)
write.table(results, file=paste('/hpc/cog_bioinf/kloosterman/users/mcretu/mating/permutations/imputedAll/perm.gonl.imputed.full.1M.',as.character(args[1]),'.txt',sep='',collapse=NULL), sep='\t', quote=FALSE, row.names=F, col.names=F)
# stopCluster(cluster)




