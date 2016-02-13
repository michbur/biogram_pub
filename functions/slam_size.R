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
  