source("./functions/simulate_sequence.R")

# define alphabet
alph <- as.character(1L:4)

# n_seq number of positive and negative sequences
# l_seq length of sequence
sim_quipt <- pblapply(1L:20, function(replication) {
  lapply(c(250, 500, 750, 1000), function(n_seq) {
    lapply(c(8, 12, 16, 20, 24), function(l_seq) {
      motifs <- generate_motif(alph, 5)
      dat <- generate_seqs(n_seq, l_seq, motifs)
      
      lapply(c("ig", "kl", "cs"), function(single_criterion) {
        res <- test_quipt(dat, n_seq, criterion, motifs)
        
        lapply(filter(res, motif)[["ngram"]], function(single_motif) {
          # n-gram of interest
          noi <- grepl(decode_ngrams(single_motif), decode_ngrams(res[["ngram"]]))
          # remove from the data the exact motif
          only_others <- res[single_motif != res[["ngram"]], ]
          
          data.frame(replication = replication,
                     l_seq = l_seq,
                     n_seq = n_seq,
                     motif = single_motif, 
                     p.value = filter(res, ngram == single_motif)[["p.value"]],
                     p.value.adj = filter(res, ngram == single_motif)[["p.value.adj"]],
                     criterion = single_criterion,
                     
                     # number of the n-grams of interest
                     nnoi = sum(noi),
                     
                     # n-grams with the motif
                     n.pos05 = sum(only_others[noi, "p.value"] < 0.05, na.rm = TRUE),
                     n.pos001 = sum(only_others[noi, "p.value"] < 0.001, na.rm = TRUE),
                     
                     # n-grams without the motif
                     n.neg05 = sum(only_others[!noi, "p.value"] < 0.05, na.rm = TRUE),
                     n.neg001 = sum(only_others[!noi, "p.value"] < 0.001, na.rm = TRUE),
                     
                     # n-grams with the motif (adjusted p-value)
                     n.pos.adj05 = sum(only_others[noi, "p.value.adj"] < 0.05, na.rm = TRUE),
                     n.pos.adj001 = sum(only_others[noi, "p.value.adj"] < 0.001, na.rm = TRUE),
                     
                     # n-grams without the motif (adjusted p-value)
                     n.neg.adj05 = sum(only_others[!noi, "p.value.adj"] < 0.05, na.rm = TRUE),
                     n.neg.adj001 = sum(only_others[!noi, "p.value.adj"] < 0.001, na.rm = TRUE)
          ) 
        }) %>% 
          do.call(rbind, .)
      }) %>% 
        do.call(rbind, .)
    }) %>% 
      do.call(rbind, .)
  }) %>% 
    do.call(rbind, .)
}) %>% 
  do.call(rbind, .)

write.csv(sim_quipt, file = "./results/sim2.csv")

cat("", file = "/home/michal/Dropbox/done.txt")
