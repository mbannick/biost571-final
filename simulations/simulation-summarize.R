##############################
# Summarize Simulation Results
##############################

library(data.table)
library(reshape2)
library(magrittr)

args <- commandArgs(trailingOnly=TRUE)
in.dir <- args[1]
out.dir <- args[2]

# READ IN THE SIMULATION RESULTS

f <- list.files(in.dir, full.names=T)
df <- lapply(f, fread) %>% rbindlist(fill=T)

# SUMMARIZE THE PERFORMANCE METRICS
# ...

# WRITE THE SUMMARY DATA

write.csv(df, paste0(out.dir, "summary.csv"), row.names=F)
