library(biogram)

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
  ns <- c(rep(2, 4), rep(3, 3))
  ds <- list(0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0))
  
  ngram_id <- sample(1L:7, 1)
  
  n <- ns[ngram_id]
  d <- ds[[ngram_id]]
  
  c(unlist(lapply(1L:length(d), function(d_id) {
    c(sample(u, 1), rep("_", d[d_id]))
  })), sample(u, 1))
}

generate_motif <- function(u, n_motif) {
  lapply(1L:n_motif, function(dummy) generate_single_motif(u))
}

# number of positive and negative sequences
n_seq <- 500
  
# length of sequence
l_seq <- 6

# define alphabet
alph <- as.character(1L:4)

# randomly generate motifs
motifs <- generate_motif(alph, 5)

# generate sequence data
test_dat <- simulate_sequences(n_seq*2, l_seq, alph, motif_l = motifs)

# perform QuiPT
test_res <- test_features(binarize(count_multigrams(test_dat,
                                                    ns = c(1, rep(2, 4), rep(3, 3)), 
                                                    ds = list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0)),
                                                              u = alph)), 
                          target = c(rep(1, n_seq), rep(0, n_seq)))
res_df <- data.frame(test_res)

res_df[["motif"]] <- res_df[["ngram"]] %in% code_ngrams(sapply(motifs, paste0, collapse = ""))
res_df
