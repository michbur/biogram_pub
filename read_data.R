library(dplyr)

raw_dat <- read.csv("./data/old_db.csv", skip = 1)
#filter
octamers <- filter(raw_dat, Date < 2013) %>%
  mutate(Description = as.character(Description)) %>%
  mutate(seq_length = nchar(Description)) %>%
  filter(seq_length == 8) 

filtered <- octamers %>%
  mutate(BindingB = grepl("Positive", Qualitative.Measure)) %>%
  group_by(Description) %>%
  summarise(conc = mean(BindingB)) %>%
  filter(conc %in% c(0, 1)) %>%
  select(Description) %>% unlist

dat <- octamers %>%
  filter(Description %in% filtered) %>%
  mutate(BindingB = grepl("Positive", Qualitative.Measure)) %>%
  group_by(Description) %>%
  summarise(target = as.logical(mean(BindingB))) 


