# Code and data for gaze fingerprinting paper


This repository has all the code and data for the analyses in Crockford et al., Detection of idiosyncratic gaze fingerprint signatures in humans.

The code directory has all of the code for running the analyses. The ```gaze_fingerprint_analysis.m``` script does analysis on the raw gaze data, computing fixation heatmaps and running all of the primary fingerprinting analysis. The ```*analysis.Rmd``` and ```*analysis.html``` files are the downstream analysis implemented in R, which reads the output of the primary analysis and does all the end analyses, statistics and figures seen in the paper.

