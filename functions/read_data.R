library(seqinr)
library(dplyr)
library(biogram)
library(cvTools)
library(ranger)
library(hmeasure)

calc_perf <- function(true.class, scores, what, threshold = 0.1) {
  HMeasure(true.class, scores, threshold = threshold)[["metrics"]][[what]]
}

calc_mcc <- function(TP, TN, FP, FN)
  (TP*TN - FP*FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))

raw_dat <- read.csv("./data/old_db.csv", skip = 1)
#filter
octamers <- filter(raw_dat, Date < 2013) %>%
  mutate(Description = as.character(Description)) %>%
  mutate(seq_length = nchar(Description)) %>%
  filter(seq_length == 8) 

only_sure <- octamers %>%
  mutate(BindingB = grepl("Positive", Qualitative.Measure)) %>%
  group_by(Description) %>%
  summarise(conc = mean(BindingB)) %>%
  filter(conc %in% c(0, 1)) %>%
  select(Description) %>% 
  unlist

# only_sure <- octamers %>%
#   group_by(Description) %>%
#   summarise(count_uni = length(unique(Qualitative.Measure))) %>%
#   filter(count_uni == 1) %>%
#   select(Description) %>% 
#   unlist


filtered <- octamers %>%
  filter(Description %in% only_sure) %>%
  mutate(BindingB = grepl("Positive", Qualitative.Measure)) %>%
  group_by(Description) %>%
  summarise(target = as.logical(mean(BindingB))) 




seqs <- select(filtered, Description) %>% 
  unlist %>% 
  strsplit("") %>% 
  do.call(rbind, .)
rownames(seqs) <- NULL

targets <- select(filtered, target) %>% unlist %>% as.numeric

bitrigrams <- as.matrix(count_multigrams(ns = c(1, rep(2, 7), rep(3, 4)), 
                                         ds = c(list(0), 0L:6, list(c(0, 0), c(0, 1), c(1, 0), c(1, 1))),
                                         seq = seqs,
                                         pos = FALSE,
                                         u = a()[-1]))

bitrigrams <- bitrigrams > 0
storage.mode(bitrigrams) <- "integer"


# set.sed(1)
# fold_list <- lapply(list(pos = which(targets == 1), neg = which(targets == 0)), function(single_n) {
#   folded <- cvFolds(length(single_n), K = 5)
#   data.frame(id = single_n[folded[["subsets"]]], which = folded[["which"]])
# })
# save(fold_list, file = "./functions/fold_list.RData")
load("./functions/fold_list.RData")

#test_features(targets, bitrigrams[, c("K_0", "L_0")])

preds <- do.call(rbind, lapply(1L:5, function(fold) {
  train_dat <- rbind(
    data.frame(bitrigrams[fold_list[["pos"]][fold_list[["pos"]][, "which"] != fold, "id"], ], tar = 1),
    data.frame(bitrigrams[fold_list[["neg"]][fold_list[["neg"]][, "which"] != fold, "id"], ], tar = 0)
  )
  
  test_dat <- rbind(
    data.frame(bitrigrams[fold_list[["pos"]][fold_list[["pos"]][, "which"] == fold, "id"], ], tar = 1),
    data.frame(bitrigrams[fold_list[["neg"]][fold_list[["neg"]][, "which"] == fold, "id"], ], tar = 0)
  )
  
  all_feats <- test_features(train_dat[, ncol(train_dat)], train_dat[, -ncol(train_dat)], adjust = NULL)
  imp_feats <- cut(all_feats, breaks = c(0, 0.01, 1))[[1]]
  
  train_dat[["tar"]] <- factor(train_dat[["tar"]], label = c("neg", "pos"))
  model <- ranger(tar ~ ., train_dat[, c(na.omit(imp_feats), "tar")], write.forest = TRUE, probability = TRUE)
  data.frame(fold = fold, pred = predict(model, test_dat)[["predictions"]][, "pos"], tar = test_dat[["tar"]])
}))

HMeasure(preds[["tar"]], preds[["pred"]])[["metrics"]]



group_by(preds, fold) %>%
  summarise(AUC = calc_perf(tar, pred, "AUC"),
            Sens = calc_perf(tar, pred, "Sens"),
            Spec = calc_perf(tar, pred, "Spec"),
            TP = calc_perf(tar, pred, "TP"),
            FP = calc_perf(tar, pred, "FP"),
            TN = calc_perf(tar, pred, "TN"),
            FN = calc_perf(tar, pred, "FN")) %>%
  mutate(MCC = calc_mcc(TP, FP, TN, FN))



  
  
