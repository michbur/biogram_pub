library(biogram)
library(dplyr)
library(pbapply)

# how it works? 
# 1. Generate unigrams and associate properties with them.
#  --each unigram may have specified level of the property 
#    (e.g., unigrams U1 has property P1 below 0.2 and P2 between 0.4 and 0.6). 
# 2. Generate rules for regions
#  --rules decide which unigrams are inside regions based on their properties
#    (e.g., region A should have aminoacids with high P1 property and low P2 property)
#  --regions must have specified length

#' @param unigram_ranges list of ranges which contain the property. If named, names 
#' are preserved.
generate_single_unigram <- function(unigram_ranges) {
  unigram_props <- lapply(unigram_ranges, function(i) {
    runif(1, min = i[1], max = i[2])
  })
  names(unigram_props) <- names(unigram_ranges)
  unlist(unigram_props)
}

generate_single_unigram(list(P1 = c(0, 0.5), 
                             P2 = c(0.2, 0.4),
                             P3 = c(0.5, 1),
                             P4 = c(0, 0)))

#' @param unigram_list named list. Names correspond to unigrams creating the 
#' alphabet.
generate_unigrams <- function(unigram_list) {
  unigrams <- do.call(rbind, lapply(unigram_list, function(single_unigram) {
    generate_single_unigram(single_unigram)
  }))
  
  unigrams
}

props1 <- list(P1 = c(0, 0.5), 
               P2 = c(0.2, 0.4),
               P3 = c(0.5, 1),
               P4 = c(0, 0))

props2 <- list(P1 = c(0.5, 1), 
               P2 = c(0.4, 1),
               P3 = c(0, 0.5),
               P4 = c(1, 1))

alph <- list(a = props1, 
             b = props1, 
             c = props1, 
             d = props2,
             e = props2)
generate_unigrams(alph)
