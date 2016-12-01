library(biogram)
library(dplyr)
library(infotheo)
library(pbapply)


generate_properties <- function(n_prop, u) {
  matrix(rnorm(20*n_prop), ncol = n_prop, dimnames = list(u))
}

sim_single_seq <- function(len, u)
  sample(u, size = len, replace = TRUE)

add_single_region <- function(seq, agg_alph, reg_len) {
  region_pos <- sample(1L:(length(seq) - reg_len + 1), 1)
  ids <- 0L:(reg_len - 1)
  # exclude gaps
  ids <- ids + region_pos
  probs <- sample(1L:length(agg_alph), reg_len, replace = TRUE, 
                  prob = c(0.7, rep(0.3/length(agg_alph), length(agg_alph) - 1)))
  region <- sapply(agg_alph[probs], sample, size = 1)
  seq[ids] <- region
  #list(seq = seq, ids = ids)
  seq
}



props1 <- generate_properties(3, letters[1L:20])
props2 <- generate_properties(7, letters[1L:20])

alph <- cutree(hclust(dist(props1, method = "euclidean"), method = "ward.D2"), 3)
agg_alph <- lapply(unique(alph), function(single_group) names(alph[alph == single_group]))

sample_seqs <- lapply(1L:10, function(i) {
  sim_single_seq(50, letters[1L:20]) %>% 
    add_single_region(agg_alph = agg_alph, reg_len = round(runif(1, 6, 15), 0))
})

raw_val <- list2matrix(sample_seqs) %>% 
  apply(1, function(ith_seq) {
    sapply(1L:(length(ith_seq) - 2), function(id1) 
      sapply((id1 + 1):(length(ith_seq) - 1), function(id2)
        mutinformation(ith_seq[c(id1, id1 + 1)], ith_seq[c(id2, id2 + 1)])
      )
    )
  }) %>% 
  unlist %>% 
  mean


red_val <- list2matrix(sample_seqs) %>% 
  degenerate(agg_alph) %>% 
  pbapply(1, function(ith_seq) {
    sapply(1L:(length(ith_seq) - 2), function(id1) 
      sapply((id1 + 1):(length(ith_seq) - 1), function(id2)
        mutinformation(ith_seq[c(id1, id1 + 1)], ith_seq[c(id2, id2 + 1)])
      )
    )
  }) %>% 
  unlist %>% 
  mean