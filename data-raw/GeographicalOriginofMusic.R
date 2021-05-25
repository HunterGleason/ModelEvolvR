## code to prepare `GeographicalOriginofMusic` dataset goes here

default_plus_chromatic_features_1059_tracks <- read.csv("data-raw/default_plus_chromatic_features_1059_tracks.txt", header=FALSE)

colnames(default_plus_chromatic_features_1059_tracks)[117]<-"latitude"
colnames(default_plus_chromatic_features_1059_tracks)[118]<-"longitude"

usethis::use_data(default_plus_chromatic_features_1059_tracks, overwrite = TRUE)
