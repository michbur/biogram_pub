library(biogram)
library(dplyr)

sim_single_seq <- function(len, u)
  sample(u, size = len, replace = TRUE)

add_single_motif <- function(motif, motif_len, seq, len) {
  motif_pos <- sample(1L:(len - motif_len + 1), 1)
  ids <- 0L:(motif_len - 1)
  # exclude gaps
  ids <- ids[motif != "_"] + motif_pos
  seq[ids] <- motif[motif != "_"]
  seq
}

# tmp <- sim_single_seq(16, alph)
# add_single_motif(c("1", "_", "1"), 3, tmp, 16)

simulate_sequences <- function(n_seq, len, u, motif_l, fraction = 0.5) {
  n_pos <- round(fraction*n_seq, 0)
  
  t(cbind(sapply(1L:n_pos, function(dummy) {
    motif <- sample(motif_l, 1)[[1]]
    add_single_motif(motif, length(motif), sim_single_seq(len, u), len)
  }),
  sapply(1L:(n_seq - n_pos), function(dummy)
    sim_single_seq(len, u)
  )
  ))
}

generate_single_motif <- function(u) {
  ns <- c(rep(2, 6), rep(3, 9))
  ds <- c(as.list(0L:5), 
          expand.grid(0L:2, 0L:2)  %>% split(1L:nrow(.)))
  
  ngram_id <- sample(1L:length(ns), 1)
  
  n <- ns[ngram_id]
  d <- ds[[ngram_id]]
  
  c(unlist(lapply(1L:length(d), function(d_id) {
    c(sample(u, 1), rep("_", d[d_id]))
  })), sample(u, 1))
}

generate_motif <- function(u, n_motif) {
  lapply(1L:n_motif, function(dummy) generate_single_motif(u))
}



# define alphabet
alph <- as.character(1L:4)


generate_seqs <- function(n_seq, l_seq, motifs) {
  # generate sequence data
  test_dat <- simulate_sequences(n_seq*2, l_seq, alph, motif_l = motifs)
  
  # perform QuiPT
  test_res <- count_multigrams(test_dat, 
                               ns = c(1, rep(2, 6), rep(3, 9)), 
                               ds = c(0, as.list(0L:5), 
                                      expand.grid(0L:2, 0L:2)  %>% split(1L:nrow(.))),
                               u = alph) %>% 
    binarize() 
}


test_quipt <- function(simulated_seqs, n_seq, criterion, motifs) {
  test_features(target = c(rep(1, n_seq), rep(0, n_seq)), features = simulated_seqs)
  res_df <- data.frame(test_res)
  
  data.frame(#n_seq = n_seq, 
    #l_seq = l_seq, 
    res_df, 
    motif = res_df[["ngram"]] %in% code_ngrams(sapply(motifs, paste0, collapse = "")),
    p.value.adj = p.adjust(res_df[["p.value"]], "BH"))
}

# number of positive and negative sequences
n_seq <- 500

# length of sequence
l_seq <- 8


# lapply(c(8, 12, 16, 20), function(l_seq) {
lapply(c(8, 12), function(l_seq) {
  motifs <- generate_motif(alph, 5)
  dat <- generate_seqs(n_seq, l_seq, motifs)
  
  lapply(c("ig", "kl", "cs"), function(single_criterion) {
    res <- test_quipt(dat, n_seq, criterion, motifs)
    
    lapply(filter(res, motif)[["ngram"]], function(single_motif) {
      # n-gram of interest
      noi <- grepl(decode_ngrams(single_motif), decode_ngrams(res[["ngram"]]))
      # remove from the data the exact motif
      only_others <- res[single_motif != res[["ngram"]], ]
      
      data.frame(l_seq = l_seq,
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
    })
  })
}) 