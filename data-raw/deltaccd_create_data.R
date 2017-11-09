library('annotate')
library('tidyverse')
library('metapredict')

####################

cg = read_csv(file.path('data-raw', 'genes_clock.csv'), col_types=cols())
cg$entrez_mm = as.character(cg$entrez_mm)
cg$entrez_hs = as.character(cg$entrez_hs)
clockGenes = cg

refCorMouseEntrez = readRDS(file.path('data-raw', 'result_mouse_ref.rds'))
refCorMouseEntrez = refCorMouseEntrez[cg$entrez_mm, cg$entrez_mm]

devtools::use_data(clockGenes, refCorMouseEntrez, internal=TRUE, overwrite=TRUE)

####################

sampleMetadataGeo = read_csv(file.path('data-raw', 'sample_metadata_cancer.csv'), col_types=cols()) %>%
	filter(condition %in% c('tumor', 'nontumor')) %>%
	mutate(condition = ifelse(condition=='nontumor', 'non-tumor', condition))
ematListGeo = extractExpressionData(readRDS(file.path('data-raw', 'proc_cancer.rds')), sampleMetadataGeo)

studyNow = 'GSE19188'
groupVec = sampleMetadataGeo$condition[sampleMetadataGeo$study==studyNow]
emat = ematListGeo[[studyNow]][,sampleMetadataGeo$sample[sampleMetadataGeo$study==studyNow]]
GSE19188 = list(emat=emat, groupVec=groupVec)

devtools::use_data(GSE19188, overwrite=TRUE)
