library(biogram)

# alphabet
alph <- as.character(1L:4)

sim_seq <- function(len, u)
  sample(u, size = len, replace = TRUE)

add_motif <- function(motif, motif_len, seq, len) {
  motif_pos <- sample(1L:(len - motif_len), 1)
  seq[motif_pos:(motif_pos + motif_len - 1)] <- motif
  seq
}

tmp <- sim_seq(16, alph)
add_motif(c("1", "1", "1"), 3, tmp, 16)
