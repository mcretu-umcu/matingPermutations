


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
add_couple_averages_and_pvals_lower <- function(permutationFrame, allData, columns, colPrefix, pvalPrefix, dataVal, regionVal, freqVal)
{
  
  for ( i in 1:length(columns) )
  {
    matesCol = paste(colPrefix, columns[i], sep="", collapse=NULL)
    allData[,columns[i]] = as.numeric(as.character(allData[,columns[i]]))
    matesVal = mean(allData[ allData$colour == "mates", columns[i]])
    permutationFrame[ permutationFrame$data == dataVal & permutationFrame$region == regionVal, matesCol ] = matesVal
    ##
    dataVector = as.numeric(as.character(permutationFrame[permutationFrame$data == dataVal & permutationFrame$region == regionVal, columns[i]]))
    # matesVal = unique(permutationFrame[ permutationFrame$data == dataVal & permutationFrame$region == regionVal, matesCol])
    # if (length(matesVal) != 1)
    # {
    #   print('Something went wrong with computing the 2 sided pvalue...stop')
    # }
    dataVector = as.numeric(as.character(c(dataVector, matesVal)))
    # dataVectorMean = mean(dataVector)
    # dataVectorSd = sd(dataVector)
    
    # dataVector = abs((dataVector - dataVectorMean) / dataVectorSd)
    # matesValNorm = abs((matesVal - dataVectorMean) / dataVectorSd)
    # dataVectorMean = mean(dataVector)
    # print (dataVectorMean)
    
    # dataVector[dataVector > dataVectorMean] = rep(2* dataVectorMean, length(dataVector[dataVector > dataVectorMean])) - dataVector[dataVector > dataVectorMean]
    ##
    # print (mean(dataVector))
    # print (dataVector[1:100])
    # ggplot() + geom_histogram (aes(dataVector))
    # ggsave(paste("/Users/mcretu/mating/permutations/test.pdf", sep=""), width=12, height = 7 )
    # pvalue = length(dataVector[dataVector <= matesVal]) / length(dataVector)
    pvalue = length(dataVector[dataVector <= matesVal]) / length(dataVector)
    ##
    
    permutationFrame[ permutationFrame$data == dataVal & permutationFrame$region == regionVal, paste(pvalPrefix, columns[i], sep="", collapse=NULL) ] = pvalue
  }
  return (permutationFrame)
}



################################################################
### int main (int argc, char **argv)

args = commandArgs(trailingOnly=TRUE)




hardPath = '/hpc/cog_bioinf/kloosterman/users/mcretu/mating/distances/gonl/'
gonl.wg.all.mhc.005 = read.table(paste(hardPath, 'wg/gonl.all.wg.maf0.005.mhc.allele.txt', sep = "", collapse=NULL), header = T, sep = '\t')


gonl.kids = read.table('/hpc/cog_bioinf/kloosterman/users/mcretu/mating/res/gonl.public.kids.ped', header = F, sep = '\t')
related_couples = c("gonl-127", "gonl-9", "gonl-173", "gonl-186", "gonl-18", "gonl-81", "gonl-83", "gonl-120", "gonl-199")
# eliminate twins; i.e.: make sure each pair of parents appears only once
gonl.kids$V1 = as.character(gonl.kids$V1)
gonl.kids$V3 = as.character(gonl.kids$V3)
gonl.kids$V4 = as.character(gonl.kids$V4)
gonl.kids = gonl.kids[ match(gonl.kids$V1, related_couples, nomatch=-1) == -1, ]
gonl.kids = gonl.kids[substr(as.character(gonl.kids$V2), nchar(as.character(gonl.kids$V2)), nchar(as.character(gonl.kids$V2))) != "d",]

gonl.wg.all.mhc.005 = colour_couples_from_ped(gonl.kids, gonl.wg.all.mhc.005, "mates", "rand")


############## !!!!!!!!!!!!! testing each region subset
regions = read.table('/hpc/cog_bioinf/kloosterman/users/mcretu/mating/res/GoNL.provinces.FINAL.txt', header=T, sep=" ")
id.map  = read.table('/hpc/cog_bioinf/kloosterman/users/mcretu/mating/res/gonl.release5.public_ids_mapping.txt', header=F, sep="\t")
regions.public = merge(regions, id.map, by.x=c("IID2"), by.y=c("V1"))
provinces = as.vector(unique(as.character(regions.public$Province)))


NR_PERM_G = 1000
columns_to_permute_G = c("qCgt", "qCal")
rest_G =  c("colour", "region", "data", "group")
couples_prefix_G = "coupleAvg"
pval_prefix_G = "pVal"
p_vals_reg = c()

# columns_to_plot_G  = columns_to_permute_G



### Special case the Admixd Dutch""
p="Admixed"
print (p)
gonl.wg.all.mhc.005.reg = gonl.wg.all.mhc.005[regions.public$Province[match(gonl.wg.all.mhc.005$ind1, regions.public$V2)] !=  regions.public$Province[match(gonl.wg.all.mhc.005$ind2, regions.public$V2)] ,]

gonl.kids.reg = gonl.kids[gonl.kids$V3%in%gonl.wg.all.mhc.005.reg$ind1 | gonl.kids$V3%in%gonl.wg.all.mhc.005.reg$ind2,]
gonl.permuted.averages.rg = as.data.frame( matrix( , 0, length(columns_to_permute_G) + 4, dimnames = list( c(), c(columns_to_permute_G, "colour", "data", "region", "group") )) )
for (i in 1:NR_PERM_G)
{
  if (i %% 100 == 0)
  {print (i)}
  recoupledG = sample_couples2( as.character(gonl.kids.reg$V3), as.character(gonl.kids.reg$V4), data.frame(m=as.character(gonl.kids.reg$V3), f=as.character(gonl.kids.reg$V4)), 1 )
  recoupledG$m = as.character(recoupledG$m)
  recoupledG$f = as.character(recoupledG$f)
  ### permute the MHC
  gonl.permuted.averages.rg = add_average_of_permuted_columns(gonl.permuted.averages.rg, gonl.wg.all.mhc.005.reg, columns_to_permute_G, recoupledG, c("rand", "wholeGenomeAll", "MHC", p))
}

###
# for (i in 1:length(columns_to_permute_G))
# {
#   gonl.permuted.averages.rg[,columns_to_permute_G[i] ] = as.numeric(as.character(gonl.permuted.averages.rg[,columns_to_permute_G[i]]))
#   gonl.permuted.averages.rg[,paste(couples_prefix_G, columns_to_permute_G[i], sep = '', collapse=NULL)] = rep(0.0,length(gonl.permuted.averages.rg[,1]))
#   gonl.permuted.averages.rg[,paste(pval_prefix_G, columns_to_permute_G[i], sep = '', collapse=NULL)] = rep(0.0,length(gonl.permuted.averages.rg[,1]))
# }

### add the averages from the real combination to the table
# gonl.permuted.averages.rg = add_couple_averages_and_pvals_lower(gonl.permuted.averages.rg, gonl.wg.all.mhc.005.reg, columns_to_permute_G, couples_prefix_G,  pval_prefix_G, "wholeGenomeAll", "MHC")


gonl.permuted.averages.provinces = gonl.permuted.averages.rg

###
# gonl.permuted.averages.provinces = as.data.frame( matrix( , 0, length(columns_to_permute_G) + 4, dimnames = list( c(), c(columns_to_permute_G, "colour", "data", "region", "group") )) )
for (p in provinces)
{
  print (p)
  gonl.wg.all.mhc.005.reg = gonl.wg.all.mhc.005[regions.public$Province[match(gonl.wg.all.mhc.005$ind1, regions.public$V2)] == p & regions.public$Province[match(gonl.wg.all.mhc.005$ind2, regions.public$V2)] == p ,]
  
  gonl.kids.reg = gonl.kids[gonl.kids$V3%in%gonl.wg.all.mhc.005.reg$ind1 | gonl.kids$V3%in%gonl.wg.all.mhc.005.reg$ind2,]
  gonl.permuted.averages.rg = as.data.frame( matrix( , 0, length(columns_to_permute_G) + 4, dimnames = list( c(), c(columns_to_permute_G, "colour", "data", "region", "group") )) )
  for (i in 1:NR_PERM_G)
  {
    if (i %% 100 == 0)
    {print (i)}
    recoupledG = sample_couples2( as.character(gonl.kids.reg$V3), as.character(gonl.kids.reg$V4), data.frame(m=as.character(gonl.kids.reg$V3), f=as.character(gonl.kids.reg$V4)), 1 )
    recoupledG$m = as.character(recoupledG$m)
    recoupledG$f = as.character(recoupledG$f)
    ### permute the MHC
    gonl.permuted.averages.rg = add_average_of_permuted_columns(gonl.permuted.averages.rg, gonl.wg.all.mhc.005.reg, columns_to_permute_G, recoupledG, c("rand", "wholeGenomeAll", "MHC", p))
  }

  ###
  # for (i in 1:length(columns_to_permute_G))
  # {
  #   gonl.permuted.averages.rg[,columns_to_permute_G[i] ] = as.numeric(as.character(gonl.permuted.averages.rg[,columns_to_permute_G[i]]))
  #   gonl.permuted.averages.rg[,paste(couples_prefix_G, columns_to_permute_G[i], sep = '', collapse=NULL)] = rep(0.0,length(gonl.permuted.averages.rg[,1]))
  #   gonl.permuted.averages.rg[,paste(pval_prefix_G, columns_to_permute_G[i], sep = '', collapse=NULL)] = rep(0.0,length(gonl.permuted.averages.rg[,1]))
  # }
  
  ### add the averages from the real combination to the table
  # gonl.permuted.averages.rg = add_couple_averages_and_pvals_lower(gonl.permuted.averages.rg, gonl.wg.all.mhc.005.reg, columns_to_permute_G, couples_prefix_G,  pval_prefix_G, "wholeGenomeAll", "MHC")
  
  # print (colnames(gonl.permuted.averages.rg))
  gonl.permuted.averages.provinces = rbind(gonl.permuted.averages.provinces, gonl.permuted.averages.rg)
  
}

############## !!!!!!!!!!!!! redistributing regions in 3 groups on the NS axis

north  = c("Friesland", "Groningen", "Drenthe", "Overijssel", "Noord_Holland")
center = c("Zuid_Holland", "Utrecht", "Gelderland")
south  = c("Zeeland", "Noord_Brabant", "Limburg")
groups = list(north, center, south)
groups.names = c("north", "center", "south")
print(groups)
print(groups[[2]])


print(provinces)
print(groups)

columns_to_plot_G = columns_to_permute_G
gonl.permuted.averages.groups = as.data.frame( matrix( , 0, length(columns_to_permute_G) + 4, dimnames = list( c(), c(columns_to_permute_G, "colour", "data", "region", "group") )) )
for (p in 1:length(groups.names)) 
{
  gonl.wg.all.mhc.005.gr = gonl.wg.all.mhc.005[regions.public$Province[match(gonl.wg.all.mhc.005$ind1, regions.public$V2)] %in% groups[[p]]  & regions.public$Province[match(gonl.wg.all.mhc.005$ind2, regions.public$V2)] %in% groups[[p]] ,]
  

  gonl.kids.gr = gonl.kids[gonl.kids$V3%in%gonl.wg.all.mhc.005.gr$ind1 | gonl.kids$V3%in%gonl.wg.all.mhc.005.gr$ind2,]
  gonl.permuted.averages.gr = as.data.frame( matrix( , 0, length(columns_to_permute_G) + 4, dimnames = list( c(), c(columns_to_permute_G, "colour", "data", "region", "group") )) )
  
  print (gonl.wg.all.mhc.005.gr[1:5,])
  print (gonl.kids.gr)
  print (nrow(gonl.wg.all.mhc.005.gr))
  print (nrow( gonl.kids.gr ))
  
  for (i in 1:NR_PERM_G)
  {
    if (i %% 100 == 0)
    {print (i)}
    recoupleG = sample_couples2( as.character(gonl.kids.gr$V3), as.character(gonl.kids.gr$V4), data.frame(m=as.character(gonl.kids.gr$V3), f=as.character(gonl.kids.gr$V4)), 1 )
    recoupleG$m = as.character(recoupleG$m)
    recoupleG$f = as.character(recoupleG$f)
    ### permute the MHC?
    gonl.permuted.averages.gr = add_average_of_permuted_columns(gonl.permuted.averages.gr, gonl.wg.all.mhc.005.gr, columns_to_permute_G, recoupleG, c("rand", "wholeGenomeAll", "MHC", groups.names[p]))
    
  }
  
  
  # for (i in 1:length(columns_to_permute_G))
  # {
  #   gonl.permuted.averages.gr[,columns_to_permute_G[i] ] = as.numeric(as.character(gonl.permuted.averages.gr[,columns_to_permute_G[i]]))
  #   gonl.permuted.averages.gr[,paste(couples_prefix_G, columns_to_permute_G[i], sep = '', collapse=NULL)] = rep(0.0,length(gonl.permuted.averages.gr[,1]))
  #   gonl.permuted.averages.gr[,paste(pval_prefix_G, columns_to_permute_G[i], sep = '', collapse=NULL)] = rep(0.0,length(gonl.permuted.averages.gr[,1]))
  # }
  
  ### add the averages from the real combination to the table
  # gonl.permuted.averages.gr = add_couple_averages_and_pvals_lower(gonl.permuted.averages.gr, gonl.wg.all.mhc.005.gr, columns_to_permute_G, couples_prefix_G,  pval_prefix_G, "wholeGenomeAll", "MHC")
  
  
  gonl.permuted.averages.provinces = rbind(gonl.permuted.averages.provinces, gonl.permuted.averages.gr)
  # print (colnames(gonl.permuted.averages.gr))
}
####

print (colnames(gonl.permuted.averages.provinces))
write.table(gonl.permuted.averages.provinces, file=paste('/hpc/cog_bioinf/kloosterman/users/mcretu/mating/permutations/',as.character(args[1]),'/perm.gonl.reg.1M.',as.character(args[2]),'.txt',sep='',collapse=NULL), sep='\t', quote=FALSE, row.names=F, col.names=F)




























