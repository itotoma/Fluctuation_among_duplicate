raw_count_data <- read.table(file.path(getwd(),"103118_SS_Data.tsv"), header=T, sep="\t")

raw_count_sum_by_deal = list()
genotype = unique(raw_count_data$Genotype_Group)
for(i in c(1:length(genotype))){
  print(genotype[i])
  raw_count_data_by_condition = raw_count_data %>% filter(Genotype_Group == genotype[i]) %>% 
    group_by(Condition) %>% group_split()
  raw_count_sum_by_condition = raw_count_data_by_condition %>%  lapply(
    function(x){
      deal = paste(genotype[i], "_",unique(x$Condition), sep="")
      print(deal)
      x %>% select(-KANMX, -NATMX, -Genotype, -Genotype_Group, -Replicate, -Condition, -tenXBarcode) %>%
        colSums() %>% as.data.frame() %>% set_names(c(deal)) %>% return()}
  ) %>% do.call('cbind', .)
  raw_count_sum_by_deal = c(raw_count_sum_by_deal, list(raw_count_sum_by_condition))
}

raw_count_sum = do.call('cbind', raw_count_sum_by_deal)

#raw_count_mean_over_deal
Eod = rowMeans(as.matrix(raw_count_sum))
names(Eod) = rownames(raw_count_sum)

#raw_count_mean_over_deal
Vod = rowVars(as.matrix(raw_count_sum))
names(Vod) = rownames(raw_count_sum)


#downloaded from Kevin B et al., Genome Res. 2005 Oct; 15(10): 1456–1461. doi: 10.1101/gr.3672305
ScerOhnologs = read.csv(file.path(getwd(),"genores_gr.3672305_Byrne_Supp_Table2.csv"))
ScerFlucOverOhnolog = ScerOhnologs %>% apply(MARGIN=1, function(row){
  Ohnolog1 = row[2]
  Ohnolog2 = row[3]
  ExpectedCV = sqrt((Vod[Ohnolog1] + Vod[Ohnolog2]))/(Eod[Ohnolog1] + Eod[Ohnolog2])
  OhnologJoitData = c(as.numeric(raw_count_sum[Ohnolog1,]), as.numeric(raw_count_sum[Ohnolog2,]))
  RealCV = sqrt(var(OhnologJoitData))/mean(OhnologJoitData)
  data.frame(id = row[1], O1 = Ohnolog1, O2 = Ohnolog2, ExpectedCV = ExpectedCV, RealCV = RealCV) %>% 
  return()
}) %>% do.call('rbind',.)
ScerFlucOverOhnolog %>%  na.omit() %>% ggplot(aes(x=ExpectedCV, y=RealCV)) + geom_point()

randompair = data.frame(No. = c(1:554), O1 = sample(ScerOhnologs$Gene.1, 554), O2 = sample(ScerOhnologs$Gene.2, 554))
ScerFlucOverRandomPair = randompair %>% apply(MARGIN=1, function(row){
  Ohnolog1 = row[2]
  Ohnolog2 = row[3]
  ExpectedCV = sqrt((Vod[Ohnolog1] + Vod[Ohnolog2]))/(Eod[Ohnolog1] + Eod[Ohnolog2])
  OhnologJoitData = c(as.numeric(raw_count_sum[Ohnolog1,]), as.numeric(raw_count_sum[Ohnolog2,]))
  RealCV = sqrt(var(OhnologJoitData))/mean(OhnologJoitData)
  data.frame(id = row[1], O1 = Ohnolog1, O2 = Ohnolog2, ExpectedCV = ExpectedCV, RealCV = RealCV) %>% 
    return()
}) %>% do.call('rbind',.)
ScerFlucOverRandomPair %>%  na.omit() %>% ggplot(aes(x=ExpectedCV, y=RealCV)) + geom_point()

#VanderSluis.B et al., Mol Syst Biol. 2010;6: 429 10.1038/msb.2010.82 
#These data reports synthetic lethal gene which undergo duplicate. 
ScerSLOhnologs <- read.csv(file.path(getwd(),"msb201082-s2.csv"), sep="\t") %>% select(c(1:10)) %>% 
  filter(interactionInSGA == "TRUE") %>% select(-CommonNameA, -CommonNameB, -SGASimilarity) %>% 
  filter(WGD == TRUE)
ScerFlucOverSLOhnolog = ScerSLOhnologs %>% apply(MARGIN=1, function(row){
  Ohnolog1 = row[1]
  Ohnolog2 = row[2]
  ExpectedCV = sqrt((Vod[Ohnolog1] + Vod[Ohnolog2]))/(Eod[Ohnolog1] + Eod[Ohnolog2])
  OhnologJoitData = c(as.numeric(raw_count_sum[Ohnolog1,]), as.numeric(raw_count_sum[Ohnolog2,]))
  RealCV = sqrt(var(OhnologJoitData))/mean(OhnologJoitData)
  data.frame(O1 = Ohnolog1, O2 = Ohnolog2, ExpectedCV = ExpectedCV, RealCV = RealCV) %>% 
    return()
}) %>% do.call('rbind',.)
ScerFlucOverSLOhnolog %>%  na.omit() %>% ggplot(aes(x=ExpectedCV, y=RealCV)) + geom_point()


cor.test(ScerFlucOverSLOhnolog$ExpectedCV, ScerFlucOverSLOhnolog$RealCV)
cor.test(ScerFlucOverDuplicate$ExpectedCV, ScerFlucOverDuplicate$RealCV)
cor.test(ScerFlucOverRandomPair$ExpectedCV, ScerFlucOverRandomPair$RealCV)




