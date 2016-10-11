library(dplyr)
library(magrittr)
library(slam)
library(knitr)

slam_size <- lapply(10^(1L:4), function(m_size) {
  m <- matrix(0, nrow = m_size, ncol = m_size)
  c(matrix = format(object.size(m), units = "B"), 
    slam = format(object.size(as.simple_triplet_matrix(m)), units = "B"))
})


unlist(slam_size) %>%
  strsplit(" ") %>% 
  unlist %>%
  matrix(ncol = 4, byrow = TRUE) %>%
  extract(, c(1, 3)) %>%
  as.numeric %>%
  matrix(ncol = 2, byrow = FALSE, dimnames = list(c(), c("matrix", "slam"))) %>%
  cbind(size = 10^(1L:4)^2, .) %>%
  data.frame() %>%
  kable
  
# real simulated example

source("./functions/simulate_sequence.R")

posssible_ngrams <- list(ns = c(1, 
                                rep(2, 6), 
                                rep(3, 9)), 
     ds = c(0, as.list(0L:5), expand.grid(0L:2, 0L:2)  %>% split(1L:nrow(.)))
)

slam_size_dat <- lapply(1L:10000, function(dummy) sim_single_seq(len = 7, u = LETTERS[1L:20])) %>% 
  do.call(rbind, .) %>% 
  count_multigrams(ns = posssible_ngrams[["ns"]], ds = posssible_ngrams[["ds"]], seq = ., u = LETTERS[1L:20])

save(slam_size_dat, file = "./results/slam_size_dat.RData")

# lapply(c(1, 7, 16), function(n_ngram)
# lapply(c(100, 500, 1000, 5000, 10000), function(n_seq) 
#   data.frame(n_gram = n_gram,
#              n_seq = n_seq,
#              normal = format(object.size(dat[1L:nseq, 1L:n_ngram]), units = "B"),
#              slam = format(object.size(as.matrix(dat[1L:nseq, 1L:n_ngram])), units = "B")
#   )
# )
# )

