library(biogram)
library(dplyr)

sim__single_seq <- function(len, u)
  sample(u, size = len, replace = TRUE)

add_single_motif <- function(motif, motif_len, seq, len) {
  motif_pos <- sample(1L:(len - motif_len + 1), 1)
  ids <- 0L:(motif_len - 1)
  # exclude gaps
  ids <- ids[motif != "_"] + motif_pos
  seq[ids] <- motif[motif != "_"]
  seq
}

# tmp <- sim__single_seq(16, alph)
# add_single_motif(c("1", "_", "1"), 3, tmp, 16)

simulate_sequences <- function(n_seq, len, u, motif_l, fraction = 0.5) {
  n_pos <- round(fraction*n_seq, 0)
  
  t(cbind(sapply(1L:n_pos, function(dummy) {
    motif <- sample(motif_l, 1)[[1]]
    add_single_motif(motif, length(motif), sim__single_seq(len, u), len)
  }),
  sapply(1L:(n_seq - n_pos), function(dummy)
    sim__single_seq(len, u)
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


test_quipt <- function(n_seq, l_seq) {
  # randomly generate motifs
  motifs <- generate_motif(alph, 5)
  
  # generate sequence data
  test_dat <- simulate_sequences(n_seq*2, l_seq, alph, motif_l = motifs)
  
  # perform QuiPT
  test_res <- count_multigrams(test_dat, 
                               ns = c(1, rep(2, 6), rep(3, 9)), 
                               ds = c(0, as.list(0L:5), 
                                      expand.grid(0L:2, 0L:2)  %>% split(1L:nrow(.))),
                               u = alph) %>% 
    binarize() %>% 
    test_features(target = c(rep(1, n_seq), rep(0, n_seq)), features = .)
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

res <- test_quipt(n_seq, l_seq)

# analyse quipt must consider p-value adjustment
