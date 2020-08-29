library(data.table)
library(tidyr)
library(splitstackshape)
library(pbapply)
library(stringr)
library(svMisc)

source("convart_helper.R")
source("data.R")
source("analysis_function.R")
source("analysis_main_function.R")


# Analysis

analysis1<-analysis_main()

mousecelegans<-analysis(analysis1, "mouse", "celegans")
mutagencelegans<-analysis(analysis1, "mutagen", "celegans")

