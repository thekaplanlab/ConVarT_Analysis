
source("analysis_main_function.R")
source("analysis_function.R")

# Indexing MSA

main_analysis<-analysis_main()

# All orthologous combinations
# Some of these function calls would take up to 2 days to finish, so it is recommended to run each in parallel, at the same time

mouse_celegans<-orthologous_analysis(main_analysis, mouse, celegans)
mutagen_celegans<-orthologous_analysis(main_analysis, mutagen, celegans)
clinvar_celegans<-orthologous_analysis(main_analysis, clinvar, celegans)
dbsnp_celegans<-orthologous_analysis(main_analysis, dbsnp, celegans)
cosmic_celegans<-orthologous_analysis(main_analysis, cosmic, celegans)
gnomad_celegans<-orthologous_analysis(main_analysis, gnomad, celegans)
mouse_clinvar<-orthologous_analysis(main_analysis, mouse, clinvar)
mouse_dbsnp<-orthologous_analysis(main_analysis, mouse, dbsnp)
mouse_cosmic<-orthologous_analysis(main_analysis, mouse, cosmic)
mouse_gnomad<-orthologous_analysis(main_analysis, mouse, gnomad)

mutagen_clinvar<-orthologous_analysis(main_analysis, mutagen, clinvar)
mutagen_dbsnp<-orthologous_analysis(main_analysis, mutagen, dbsnp)
mutagen_cosmic<-orthologous_analysis(main_analysis, mutagen, cosmic)
mutagen_gnomad<-orthologous_analysis(main_analysis, mutagen, gnomad)


# Indexing for double MSA

main_analysis_double<-analysis_main_double()

# Orthologous variants from protein only found in human-mouse or human- C. elegans

d_mousecelegans<-double_orthologous_analysis(main_analysis_double, mouse, celegans)
d_mutagencelegans<-double_orthologous_analysis(main_analysis_double, mutagen, celegans)
d_clinvarcelegans<-double_orthologous_analysis(main_analysis_double, clinvar, celegans)
d_dbsnpcelegans<-double_orthologous_analysis(main_analysis_double, dbsnp, celegans)
d_cosmiccelegans<-double_orthologous_analysis(main_analysis_double, cosmic, celegans)
d_gnomadcelegans<-double_orthologous_analysis(main_analysis_double, gnomad, celegans)
d_mouseclinvar<-double_orthologous_analysis(main_analysis_double, mouse, clinvar)
d_mousedbsnp<-double_orthologous_analysis(main_analysis_double, mouse, dbsnp)
d_mousecosmic<-double_orthologous_analysis(main_analysis_double, mouse, cosmic)
d_mousegnomad<-double_orthologous_analysis(main_analysis_double, mouse, gnomad)

d_mutagenclinvar<-double_orthologous_analysis(main_analysis_double, mutagen, clinvar)
d_mutagendbsnp<-double_orthologous_analysis(main_analysis_double, mutagen, dbsnp)
d_mutagencosmic<-double_orthologous_analysis(main_analysis_double, mutagen, cosmic)
d_mutagengnomad<-double_orthologous_analysis(main_analysis_double, mutagen, gnomad)

