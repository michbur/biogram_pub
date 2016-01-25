library(dplyr)

raw_dat <- read.csv("./data/old_db.csv", skip = 1)
#filter
dat <- filter(raw_dat, Date < 2013) %>%
  mutate(Description = as.character(Description)) %>%
  mutate(seq_length = nchar(Description)) %>%
  filter(seq_length == 8) %>%
  mutate(BindingB = grepl("Positive", dat[["Qualitative.Measure"]])) %>%
  group_by(Description)

# if concordance is larger tha 0 and smaller than 1, octamer was classified as both positive and negative
concordance <- sapply(unique(dat[["Description"]]), function(single_octamer) {
  bindings <- dat[dat[["Description"]] == single_octamer, "BindingB"]
  sum(bindings)/length(bindings)
})

sum(dat[["Description"]] %in% names(which(concordance > 0 & concordance < 1)))

processed_dat <- filter(dat, !(Description %in% names(which(concordance > 0 & concordance < 1))))


dat[["Description"]]
length(dat[["Qualitative.Measure"]])
dat[["Date"]]

table(dat[["Date"]])
