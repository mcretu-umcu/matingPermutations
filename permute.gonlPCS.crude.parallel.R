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
add_pcs <- function (regress.table, gonl.pcs.public, pcs_to_use, ind1Prefix, ind2Prefix)
{
  for (pc in pcs_to_use)
  {
    regress.table[,paste(ind1Prefix,pc, sep="", collapse=NULL)] = gonl.pcs.public[match(regress.table$ind1, gonl.pcs.public[,"V2"]), pc]
    regress.table[,paste(ind2Prefix,pc, sep="", collapse=NULL)] = gonl.pcs.public[match(regress.table$ind2, gonl.pcs.public[,"V2"]), pc]
  }
  
  return (regress.table)
}
###
add_couple_pcs <- function (regress.table, pcs_to_use, ind1Prefix, ind2Prefix, couplePrefix)
{
  for (pc in pcs_to_use)
  {
    regress.table[,paste(couplePrefix,pc,"avg", sep="", collapse=NULL)] = ( regress.table[,paste(ind1Prefix,pc, sep="", collapse=NULL)] + regress.table[,paste(ind2Prefix,pc, sep="", collapse=NULL)] ) / 2.0    
    regress.table[,paste(couplePrefix,pc,"se", sep="", collapse=NULL)] = ( regress.table[,paste(ind1Prefix,pc, sep="", collapse=NULL)] - regress.table[,paste(ind2Prefix,pc, sep="", collapse=NULL)] ) ** 2.0    
    
  }
  return (regress.table)
}
###
compute_regression_residuals <- function(regress, pcs_to_use, couplePrefix, lm1, pc_combination, target_column)
{
  coeffs = as.data.frame(t(lm1[[1]]))
  idx_data = match( paste(couplePrefix,pcs_to_use,pc_combination, sep="", collapse=NULL), colnames(regress) )
  idx_coeffs = match( paste(couplePrefix,pcs_to_use,pc_combination, sep="", collapse=NULL), colnames(coeffs) )
  idx_target = match(target_column, colnames(regress))
  
  regress[,paste("qCgtRez",pc_combination, sep="", collapse=NULL)] = rep(-1,length(regress[,1]))
  idx_rez = match(paste("qCgtRez",pc_combination, sep="", collapse=NULL),colnames(regress))
  #print (colnames(regress))
  
  test = apply(regress,1,function(x){ estimate=coeffs[[1]]; for(j in 1:length(idx_data)){estimate = estimate + as.numeric(x[idx_data[j]]) * as.numeric(coeffs[[idx_coeffs[j]]]);}; x[idx_rez] = as.numeric(as.character(x[idx_target])) - estimate; return(x) } )
  # test = apply(regress.table,1,function(x){ estimate=0.0; for(j in 1:length(idx_data)){print (as.numeric(x[idx_data[j]])) } } )
  test1 = as.data.frame(t(test))
  test1 = typeset_cols_PCs(test1, pcs_to_use)
  
  return (test1)
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
    gonl.permuted.averages.var = add_average_of_permuted_columns(gonl.permuted.averages.var, gonl.wg.all.mhc.regress, columns_to_permute.var, recoupledG.var, c("rand", "wholeGenomeAll", "MHC", "maf0.5"))

  }
  return (gonl.permuted.averages.var)
}
#######
args = commandArgs(trailingOnly=TRUE)

maf = "0.01"
hardPath = '/hpc/cog_bioinf/kloosterman/users/mcretu/mating/'
gonl.wg.all.mhc.005 = read.table(paste(hardPath, 'distances/gonl/wg/gonl.all.wg.maf0.005.mhc.allele.txt', sep = "", collapse=NULL), header = T, sep = '\t')
gonl.kids = read.table(paste(hardPath,'res/gonl.public.kids.ped',sep="", collapse=NULL), header = F, sep = '\t')
related_couples = c("gonl-127", "gonl-9", "gonl-173", "gonl-186", "gonl-18", "gonl-81", "gonl-83", "gonl-120", "gonl-199")

gonl.pcs = read.table(paste(hardPath,'res/gonl.genome.release4.all.PCA.v2.pcs.txt', sep="", collapse=NULL), sep=" ", header=T)

regions = read.table(paste(hardPath,'res/GoNL.provinces.FINAL.txt',sep="", collapse=NULL), header=T, sep=" ")
id.map  = read.table(paste(hardPath,'res/gonl.release5.public_ids_mapping.txt',sep="", collapse=NULL), header=F, sep="\t")
regions.public = merge(regions, id.map, by.x=c("IID2"), by.y=c("V1"))

gonl.pcs.public = merge( gonl.pcs, id.map, by.x = c("IID"), by.y=c("V1") )
# eliminate twins; i.e.: make sure each pair of parents appears only once
gonl.kids$V1 = as.character(gonl.kids$V1)
gonl.kids$V3 = as.character(gonl.kids$V3)
gonl.kids$V4 = as.character(gonl.kids$V4)
gonl.kids = gonl.kids[ match(gonl.kids$V1, related_couples, nomatch=-1) == -1, ]

gonl.kids = gonl.kids[substr(as.character(gonl.kids$V2), nchar(as.character(gonl.kids$V2)), nchar(as.character(gonl.kids$V2))) != "d",]
gonl.kids.regress = gonl.kids[gonl.kids$V3 %in% gonl.pcs.public$V2 & gonl.kids$V4 %in% gonl.pcs.public$V2,]


# test = gonl.pcs.public[substr( as.character(gonl.pcs.public$V2), nchar(as.character(gonl.pcs.public$V2)), nchar(as.character(gonl.pcs.public$V2)) ) != "c",]

## Keep only families where we have computed PCs for both parents !!! (keeping it trace-able!)
# pc_families = unique(as.character(substr( as.character(gonl.pcs.public$V2), 1, nchar(as.character(gonl.pcs.public$V2)) - 1 )))
# pc_families = pc_families[paste(pc_families,"a", sep="", collapse=NULL) %in% as.character(gonl.pcs.public$V2) & paste(pc_families,"b", sep="", collapse=NULL) %in% as.character(gonl.pcs.public$V2)]

######################## regress PCs out, permute and plot : 
# gonl.wg.all.mhc.005 = read.table(paste(hardPath, '/gonl.all.wg.maf0.005.mhc.allele.txt', sep = "", collapse=NULL), header = T, sep = '\t')

gonl.wg.all.mhc.regress = gonl.wg.all.mhc.005[ gonl.wg.all.mhc.005$ind1 %in% gonl.kids.regress$V3 & gonl.wg.all.mhc.005$ind2 %in% gonl.kids.regress$V4 | gonl.wg.all.mhc.005$ind2 %in% gonl.kids.regress$V3 & gonl.wg.all.mhc.005$ind1 %in% gonl.kids.regress$V4,]           
# gonl.leftovers = gonl.wg.all.mhc.005[!(substr(as.character(gonl.wg.all.mhc.005$ind1), 1, nchar(as.character(gonl.wg.all.mhc.005$ind1)) - 1 ) %in% pc_families & substr(as.character(gonl.wg.all.mhc.005$ind2), 1, nchar(as.character(gonl.wg.all.mhc.005$ind2)) - 1 ) %in% pc_families),]           
# gonl.leftovers = colour_couples_from_ped(gonl.kids, gonl.leftovers, "mates", "rand")
# nrow(gonl.leftovers[gonl.leftovers$colour == "mates",])

gonl.wg.all.mhc.regress = colour_couples_from_ped(gonl.kids.regress, gonl.wg.all.mhc.regress, "mates", "rand")

pcs_to_use = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
ind1Prefix = "i1"
ind2Prefix = "i2"
couplePrefix = "couple"
target_column = "qCgt"

############## get PCs
gonl.wg.all.mhc.regress = add_pcs(gonl.wg.all.mhc.regress, gonl.pcs.public, pcs_to_use, ind1Prefix, ind2Prefix)

############## compute PC values for couples
gonl.wg.all.mhc.regress = add_couple_pcs(gonl.wg.all.mhc.regress, pcs_to_use, ind1Prefix, ind2Prefix, couplePrefix)


############## regress
# regress.table$class = rep(0, length(regress.table[,1]))
# regress.table$class[regress.table$colour == "mates"] = 1

# +couplePC2avg+couplePC4avg+couplePC10avg
# +couplePC10se+couplePC4se+couplePC2se
formula1 = paste(target_column, "~couplePC1avg+couplePC2avg+couplePC3avg+couplePC4avg+couplePC5avg+couplePC6avg+couplePC7avg+couplePC8avg+couplePC9avg+couplePC10avg", sep="", colapse=NULL)
formula2 = paste(target_column, "~couplePC1se+couplePC2se+couplePC3se+couplePC4se+couplePC5se+couplePC6se+couplePC7se+couplePC8se+couplePC9se+couplePC10se", sep="", colapse=NULL)

###### sample a permutation for training:
recoupleG = sample_couples2( as.character(gonl.kids.regress$V3), as.character(gonl.kids.regress$V4), data.frame(m=as.character(gonl.kids.regress$V3), f=as.character(gonl.kids.regress$V4)), 1 )
recoupleG$m = as.character(recoupleG$m)
recoupleG$f = as.character(recoupleG$f)

# select one random permutation for training
gonl.wg.all.mhc.regress.train = na.omit(gonl.wg.all.mhc.regress[match(gonl.wg.all.mhc.regress$ind1, recoupleG$m) == match(gonl.wg.all.mhc.regress$ind2, recoupleG$f) | match(gonl.wg.all.mhc.regress$ind2, recoupleG$m) == match(gonl.wg.all.mhc.regress$ind1, recoupleG$f), ])



# use first method for combining PCs
regress.formula = formula1
lm.wg.all.mhc = glm(regress.formula , data=gonl.wg.all.mhc.regress.train, family=gaussian)
############## add column containing regression residuals
gonl.wg.all.mhc.regress = compute_regression_residuals(gonl.wg.all.mhc.regress, pcs_to_use, couplePrefix, lm.wg.all.mhc, "avg", "qCgt")


### use second method for combining PCs
regress.formula = formula2
lm.wg.all.mhc = glm(regress.formula , data=gonl.wg.all.mhc.regress.train, family=gaussian)
############## add column containing regression residuals
gonl.wg.all.mhc.regress = compute_regression_residuals(gonl.wg.all.mhc.regress, pcs_to_use, couplePrefix, lm.wg.all.mhc, "se", "qCgt")



######## Permute for p-value
NR_PERM_G = 1000
columns_to_permute = c("qCgtRezavg","qCgtRezse")
couples_prefix = "coupleAvg"
pval_prefix = "pVal"
columns_to_plot = columns_to_permute

# cluster = makeCluster(25,type="SOCK")
# registerDoSNOW(cluster)
gonl.permuted.averages = as.data.frame( matrix( , 0, 4 + length(columns_to_permute), dimnames = list( c(), c(columns_to_permute, "colour", "data", "region", "freq") )) )

# results = foreach (i = 1:NR_PERM_G, .combine=rbind) %do% perform_one_permutation(gonl.wg.mhc, gonl.wg.mhcX, gonl.wg.wg, gonl.hm.mhc, gonl.hm.mhcX, gonl.hm.wg, columns_to_permute, gonl.kids,1000)
results = perform_one_permutation( columns_to_permute, gonl.kids.regress, 1000)

colnames (results) = colnames(gonl.permuted.averages)
write.table(results, file=paste('/hpc/cog_bioinf/kloosterman/users/mcretu/mating/permutations/pcs/perm.gonl.pcs.full.1M.',as.character(args[1]),'.txt',sep='',collapse=NULL), sep='\t', quote=FALSE, row.names=F, col.names=F)
# stopCluster(cluster)




