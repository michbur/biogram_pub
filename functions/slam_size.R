library(dplyr)
library(magrittr)
library(slam)
library(knitr)
library(reshape2)

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

# slam_size_dat <- lapply(1L:10000, function(dummy) sim_single_seq(len = 7, u = LETTERS[1L:20])) %>% 
#   do.call(rbind, .) %>% 
#   count_multigrams(ns = posssible_ngrams[["ns"]], ds = posssible_ngrams[["ds"]], seq = ., u = LETTERS[1L:20])
#save(slam_size_dat, file = "./results/slam_size_dat.RData")

load("./results/slam_size_dat.RData")

slam_summary <- lapply(c(20, 400*6 + 20), function(n_ngram)
  lapply(1L:10*1000, function(n_seq)
    data.frame(n_ngram = n_ngram,
               n_seq = n_seq,
               slam = format(object.size(slam_size_dat[1L:n_seq, 1L:n_ngram]), units = "B"),
               normal = format(object.size(as.matrix(slam_size_dat[1L:n_seq, 1L:n_ngram])), units = "B")
    )
  ) %>% do.call(rbind, .)
) %>% do.call(rbind, .) %>% 
  mutate(ngrams = factor(n_ngram, labels = c("1-grams", "1- to 2-grams"))) %>% 
  select(-n_ngram) %>% 
  melt(id.vars = c("n_seq", "ngrams"), 
       variable.name = "method", value.name = "size") %>% 
  mutate(size = as.numeric(as.character(sub(" bytes", "", size))))

format(object.size(slam_size_dat[, 1L:(400*6 + 20)]), units = "MB")
# 4.2 Mb
format(object.size(as.matrix(slam_size_dat[, 1L:(400*6 + 20)])), units = "MB")
# 184.8 Mb

save(slam_summary, file = "./results/slam_summary.RData")

# library(ggplot2)
# ggplot(slam_summary, aes(x = n_seq, y = size, color = method)) +
#   geom_point() +
#   scale_x_continuous("Number of sequences") +
#   scale_y_continuous("Size [bytes]") +
#   scale_color_discrete("Matrix type") +
#   facet_wrap(~ ngrams, ncol = 1, scales = "free_y")
# 
