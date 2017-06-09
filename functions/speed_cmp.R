library(biogram)


dat1 <- create_feature_target(25, 25, 25, 25) 

times <- lapply(c(50, 100, 150, 200, 250, 300), function(sample_size) {
  dat2 <- sapply(1L:10, function(dummy) 
    sample(0L:1, size = 100, replace = TRUE))
  lapply(1L:10, function(dummy) {
    list(slow = system.time(test_features(dat1[, "tar"], dat2, quick = FALSE, times = 1e+05)),
         quick = system.time(test_features(dat1[, "tar"], dat2, quick = TRUE, times = 1e+05)))
  })
})

save(times, file = "times.RData")
