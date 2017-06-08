library(dplyr)
library(biogram)
library(ggplot2)

sim_quipt <- read.csv("./results/sim_evaluation50.csv") %>%
  mutate(motif_len = nchar(gsub(".", "", ngrams2df(as.character(motif))[["ngram"]], fixed = TRUE))) 

# p-values of true motifs
tm_pvals <- select(sim_quipt, l_seq, n_seq, criterion, p.value, motif_len)
ggplot(tm_pvals, aes(x = criterion, y = p.value, color = as.factor(motif_len))) +
  geom_boxplot() +
  facet_grid(l_seq ~ n_seq, labeller = label_both)
# its harder to detect the true motif the less sequences we have and the longer they are


filter(sim_quipt, motif_len == 2) %>% 
  mutate(frac05 = n.pos05/nnoi) %>% 
  ggplot(aes(x = criterion, y = frac05)) +
  geom_boxplot() +
  facet_grid(l_seq ~ n_seq, labeller = label_both)

# its harder to find motifs containing the true motif the less sequences we have and the longer they are

# no differences between citerions! how to create a scenario when difference will be visible?

filter(sim_quipt, criterion == "ig") %>% 
  select(l_seq, n_seq, criterion, p.value, motif_len) %>% 
  group_by(l_seq, n_seq, motif_len) %>% 
  summarise(acc = mean(p.value < 0.05)) %>% 
  ggplot(aes(x = factor(l_seq), y = factor(n_seq), fill = acc)) +
  geom_tile(color = "black") +
  facet_wrap(~ motif_len)
