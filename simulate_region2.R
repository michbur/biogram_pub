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
    round(runif(1, min = i[1], max = i[2]), 2)
  })
  names(unigram_props) <- names(unigram_ranges)
  unlist(unigram_props)
}

generate_single_unigram(list(P1 = c(0, 0.5), 
                             P2 = c(0.2, 0.4),
                             P3 = c(0.5, 1),
                             P4 = c(0, 0)))

#' @param unigram_list a list of unigrams' parameters. See Details.
#' @param unigram_names names of unigrams. If not \code{NULL}, will
#' overwrite any existing unigram names.
#' @param prop_names names of properties. If not \code{NULL}, will 
#' overwrite any existing names.
#' @details Unigram parameters are represented as a list of ranges which contain 
#' the property. All list of ranges should have the same length, which is an 
#' equivalent of describing each unigram using the same properties. 
generate_unigrams <- function(unigram_list,
                              unigram_names = NULL,
                              prop_names = NULL) {
  if(length(unique(lengths(unigram_list))) != 1)
    stop("All unigrams must be defined by the same number of properties (each element of \
         unigram_list must have the same length).")
  
  unigrams <- do.call(cbind, lapply(unigram_list, function(single_unigram) {
    generate_single_unigram(single_unigram)
  }))

  if(!is.null(prop_names)) 
    rownames(unigrams) <- prop_names
    
  if(!is.null(unigram_names)) 
    colnames(unigrams) <- unigram_names
  
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


alph <- generate_unigrams(c(replicate(8, props1, simplify = FALSE),
                    replicate(12, props2, simplify = FALSE)),
                    unigram_names = letters[1L:20])


#' @param alphabet the unigram alphabet. Columns are equivalent to unigrams 
#' and rows to particular properties.
#' @param reg_len the number of unigrams inside the region.
#' @param prop_rules required intervals of properties of unigrams in the region. 
#' See Detailes.
generate_single_region <- function(alphabet, reg_len, prop_ranges) {
  
  unigrams <- colnames(alphabet)
  min_range <- sapply(prop_ranges, function(i) i[1])
  max_range <- sapply(prop_ranges, function(i) i[2])
  
  unigrams <- names(which(apply(alphabet >= min_range & alphabet <= max_range, 2, all)))
  
  sample(unigrams, reg_len, replace = TRUE)
}

rules1 <- list(
    P1 = c(0.5, 1), 
    P2 = c(0.4, 1),
    P3 = c(0, 0.5),
    P4 = c(1, 1))

generate_single_region(alph, 10, rules1)

# simulation ------------------------------

rules_set <- list(c(0, 0.5),
                  c(0.5, 1),
                  c(0.25, 0.75))


# three groups of AAs
set.seed(15390)
aa_props1 <- sample(rules_set, 6, replace = TRUE)
aa_props2 <- sample(rules_set, 6, replace = TRUE)
aa_props3 <- sample(rules_set, 6, replace = TRUE)

alph <- generate_unigrams(c(replicate(7, aa_props1, simplify = FALSE),
                            replicate(7, aa_props2, simplify = FALSE),
                            replicate(6, aa_props3, simplify = FALSE)),
                          unigram_names = letters[1L:20])
rownames(alph) <- paste0("P", 1L:6)

library(dplyr)
library(reshape2)
library(ggplot2)

data.frame(alph) %>% 
  mutate(prop = rownames(.)) %>% 
  melt %>% 
  ggplot(aes(x = variable, y = value)) +
  geom_point() +
  facet_wrap(~ prop)

# two regions

reg1_rules <- list(
  P1 = c(0.5, 1), 
  P2 = c(0, 0.25),
  P3 = c(0.5, 1),
  P4 = c(0, 1),
  P5 = c(0, 1),
  P6 = c(0, 1))

reg2_rulesA <- list(
  P1 = c(0, 0.5), 
  P2 = c(0.25, 1),
  P3 = c(0, 0.5),
  P4 = c(0, 1),
  P5 = c(0, 1),
  P6 = c(0, 1))

reg2_rulesB <- list(
  P1 = c(0.25, 0.75), 
  P2 = c(0.25, 1),
  P3 = c(0, 0.75),
  P4 = c(0.5, 1),
  P5 = c(0, 1),
  P6 = c(0, 1))

generate_single_region(alph, 10, reg1_rules)
generate_single_region(alph, 10, reg2_rulesA)
generate_single_region(alph, 10, reg2_rulesA)
