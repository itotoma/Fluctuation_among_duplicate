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

#downloaded from Kevin B et al., Genome Res. 2005 Oct; 15(10): 1456–1461. doi: 10.1101/gr.3672305
ScerOhnologs = read.csv(file.path(getwd(),"genores_gr.3672305_Byrne_Supp_Table2.csv")) %>% select(-No.)

#VanderSluis.B et al., Mol Syst Biol. 2010;6: 429 10.1038/msb.2010.82 
#These data reports synthetic lethal gene which undergo duplicate. 
ScerSLOhnologs <- read.csv(file.path(getwd(),"msb201082-s2.csv"), sep="\t") %>% select(c(1:10)) %>% 
  filter(interactionInSGA == "TRUE") %>% select(-CommonNameA, -CommonNameB, -SGASimilarity) %>% 
  filter(WGD == TRUE)

FlucAmongDuplicateFrom <- function(countdata, ohnologlist){
  #raw_count_mean_over_deal
  Eod = rowMeans(as.matrix(countdata))
  names(Eod) = rownames(countdata)
  #raw_count_mean_over_deal
  Vod = rowVars(as.matrix(countdata))
  names(Vod) = rownames(countdata)
  ScerFlucOverOhnolog = ohnologlist %>% apply(MARGIN=1, function(row){
    Ohnolog1 = row[1]
    Ohnolog2 = row[2]
    if(is.element(Ohnolog1, rownames(countdata))){
      if(is.element(Ohnolog2, rownames(countdata))){
        ExpectedCV = sqrt((Vod[Ohnolog1] + Vod[Ohnolog2]))/(Eod[Ohnolog1] + Eod[Ohnolog2])
        O1CV = sqrt(Vod[Ohnolog1])/Eod[Ohnolog1]
        O2CV = sqrt(Vod[Ohnolog2])/Eod[Ohnolog2]
        maxCV = max(O1CV, O2CV)
        ExpectedCV = sqrt(Vod[Ohnolog1] + Vod[Ohnolog2])/(Eod[Ohnolog1] + Eod[Ohnolog2])
        OhnologJoitData = c(as.numeric(countdata[Ohnolog1,]) + as.numeric(countdata[Ohnolog2,]))
        ArithmeticCV = sqrt(var(OhnologJoitData))/mean(OhnologJoitData)
        DiffCV = ArithmeticCV -  maxCV
        data.frame(O1 = Ohnolog1, O2 = Ohnolog2, O1CV = O1CV, O2CV = O2CV, ArithmeticCV = ArithmeticCV, ExpectedCV = ExpectedCV,  maxCV = maxCV, DiffCV = DiffCV) %>% 
          return()
      }
    }
  }) %>% do.call('rbind',.) %>% 
    mutate(chi = (ArithmeticCV-ExpectedCV)^2/ExpectedCV) %>% 
    return()
}
#WTCountFlucAmongDuplicate = FlucAmongDuplicateFrom(t(WTCount), ScerOhnologs)
#ScerFlucOverRandomPair = FlucAmongDuplicateFrom(t(WTCount), randompair)

#Comparision of random pair and ohnolog pair
WTCountFlucAmongDuplicate = FlucAmongDuplicateFrom(raw_count_sum, ScerOhnologs)
randompair = data.frame(O1 = sample(rownames(raw_count_sum), 554, replace = F), O2 = sample(rownames(raw_count_sum), 554, replace = F))
ScerFlucOverRandomPair = FlucAmongDuplicateFrom(raw_count_sum, randompair)
CompOhnoRand = rbind(mutate(ScerFlucOverRandomPair, OhnologPair = FALSE), mutate(WTCountFlucAmongDuplicate,  OhnologPair=TRUE))
CompOhnoRand %>% 
  ggplot(aes(x=ArithmeticCV, y=maxCV)) + geom_point(aes(colour= OhnologPair)) + xlab("Arithmetic CV")+ geom_abline(intercept=0,slope=1, alpha=0.5) +
  ggtitle('ArthmeticCV vs Max CV')
CompOhnoRand %>% 
  ggplot(aes(x=DiffCV, fill=OhnologPair)) + geom_histogram(position="identity", alpha=0.8)


#Comparision of synthetic lethal ohnolog and ohnolog
CompSLOhno = WTCountFlucAmongDuplicate %>% 
  mutate(SyntheticLethality = ifelse(is.element(O1, c(ScerSLOhnologs$ORF.A, ScerSLOhnologs$ORF.B)), TRUE, FALSE)) %>% 
  distinct(O1, .keep_all = TRUE)
CompSLOhno %>% 
  ggplot(aes(x=ArithmeticCV, y=maxCV)) + geom_point(aes(colour=SyntheticLethality)) + geom_abline(intercept=0,slope=1, alpha=0.5)  
CompSLOhno %>% 
  ggplot(aes(x=chi, fill=SyntheticLethality)) + geom_histogram(position = 'identity')




#合成致死遺伝子 | 少なくともどちらかが量的近郊遺伝子 | ペアの間に物理的相互作用がない|
tmpA = ScerSLDuplicate %>% filter(is.element(ORF.A,ScerDBG$LocusName)) %>% select(ORF.A, ORF.B)
tmpB = ScerSLDuplicate %>% filter(is.element(ORF.B,ScerDBG$LocusName)) %>% select(ORF.A, ORF.B)
#SL and At Least One DBG
SL_ALODBG = rbind(tmpA, tmpB)

#Comparision Ohnolog and (Synthetic Lethal & At least one of them is DBG)
CompSLDBGOhno =  CompSLOhno %>% 
  mutate(SLDBG = ifelse(is.element(O1, c(SL_ALODBG$ORF.A, SL_ALODBG$ORF.B)), TRUE, FALSE)) %>% 
  distinct(O1, .keep_all = TRUE) %>% filter(SyntheticLethality == TRUE)
CompSLDBGOhno %>% 
  ggplot(aes(x=ArithmeticCV, y=maxCV)) + geom_point(aes(colour=SLDBG)) + geom_abline(intercept=0,slope=1, alpha=0.5)  



SL_ALSDBGFluc = FlucAmongDuplicateFrom(raw_count_sum, SL_ALODBG) %>% 
  mutate(PhyInt = ifelse(O1 %in% c(RedunduntCandidateDBG$ORF.A, RedunduntCandidateDBG$ORF.B), FALSE, TRUE)) %>% 
  distinct(O1, .keep_all=TRUE)
SL_ALSDBGFluc %>% 
  #ggplot(aes(x=ArithmeticCV, y=maxCV)) + geom_point(aes(colour=PhyInt)) + xlab("Arithmetic CV")+ geom_abline(intercept=0,slope=1, alpha=0.5) 
  ggplot(aes(x = CV, fill = PhyInt)) + geom_histogram(position = "identity", alpha = 0.8) + 
  ggtitle("Ohnolog & Synthetic lethal pair & at least one is DBG")



ScerFlucOverOhnolog = FlucAmongDuplicateFrom(raw_count_sum, ScerOhnologs)
ScerFlucOverOhnolog %>%  na.omit() %>% ggplot(aes(x=(O1CV+O2CV)/2, y=RealCV)) + geom_point() + xlim(0,5) + ylim(0,5)




ScerFlucOverSLOhnolog = FlucAmongDuplicateFrom(raw_count_sum, ScerSLOhnologs)
ScerFlucOverSLOhnolog = FlucAmongDuplicateFrom(t(WTCount), ScerSLOhnologs)
ScerFlucOverSLOhnolog %>%  na.omit() %>% ggplot(aes(x=ExpectedCV, y=RealCV)) + geom_point()








cor.test(ScerFlucOverSLOhnolog$ExpectedCV, ScerFlucOverSLOhnolog$RealCV)
cor.test(ScerFlucOverDuplicate$ExpectedCV, ScerFlucOverDuplicate$RealCV)
cor.test(ScerFlucOverRandomPair$ExpectedCV, ScerFlucOverRandomPair$RealCV)

ScerFlucOverOhnolog$chi %>% na.omit() %>% hist()
ScerFlucOverRandomPair$chi %>% na.omit() %>% hist()
ScerFlucOverSLOhnolog$chi %>% na.omit() %>% hist()


ScerFlucOverOhnolog %>% 
  mutate(SL = ifelse(is.element(ScerFlucOverOhnolog$O1,c(ScerSLOhnologs$ORF.A, ScerSLOhnologs$ORF.B)),TRUE, FALSE)) %>%
  ggplot(aes(x = RealCV, fill = SL)) + geom_histogram(position = "identity", alpha = 0.8)


ScerFlucOverOhnolog %>% 
  mutate(Redundunt = ifelse(is.element(ScerFlucOverOhnolog$O1,c(RedunduntCandidate$ORF.A, RedunduntCandidate$ORF.B)),TRUE, FALSE)) %>%
  ggplot(aes(x = chi, fill = Redundunt)) + geom_histogram(position = "identity", alpha = 0.8)



test =  WTCount %>% colSums() %>% as.data.frame() %>% rownames_to_column() %>% 
  setNames(c("ORF", "Sum"))

test[is.element(test$ORF, c(NonRedunduntCandidate$ORF.A, NonRedunduntCandidate$ORF.B)),] %>% 
  arrange(desc(Sum))

ScerSLDuplicate %>% filter(ORF.B == 'YLL045C')
WTCount['YHL033C']$YHL033C %>% hist()
WTCount['YLL045C']$YLL045C %>% hist()
c(WTCount['YHL033C']$YHL033C + WTCount['YLL045C']$YLL045C) %>% hist()

sd(WTCount['YHL033C']$YHL033C)/mean(WTCount['YHL033C']$YHL033C)
sd(WTCount['YLL045C']$YLL045C)/mean(WTCount['YLL045C']$YLL045C)

sd(WTCount['YHL033C']$YHL033C + WTCount['YLL045C']$YLL045C)/mean(WTCount['YHL033C']$YHL033C + WTCount['YLL045C']$YLL045C)
sqrt( var(WTCount['YHL033C']$YHL033C)+var(WTCount['YLL045C']$YLL045C) )/(mean(WTCount['YHL033C']$YHL033C)+mean(WTCount['YLL045C']$YLL045C))


ScerFlucOverOhnolog %>% filter(O1 == 'YHL033C')


