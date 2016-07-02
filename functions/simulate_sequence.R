library(biogram)

# alphabet
alph <- as.character(2L:4)

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

tmp <- sim__single_seq(16, alph)
add_single_motif(c("1", "_", "1"), 3, tmp, 16)

simulate_sequences <- function(n_seq, len, u, motif, fraction = 0.5) {
  motif_len <- length(motif)
  n_pos <- round(fraction*n_seq, 0)
  
  t(cbind(sapply(1L:n_pos, function(dummy)
    add_single_motif(motif, motif_len, sim__single_seq(len, u), len)
    ),
    sapply(1L:(n_seq - n_pos), function(dummy)
      sim__single_seq(len, u)
      )
  ))
}

simulate_sequences(50, 6, alph, c("1", "_", "1"))
