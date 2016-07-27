library(biogram)
test_features(target = c(rep(1, 500), rep(0, 500)), features = read.csv("problematic_data.csv"))
