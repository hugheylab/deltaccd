library('annotate')
library('data.table')
library('metapredict')

####################

cg = fread(file.path('data-raw', 'genes_clock.csv'))
cg[, entrez_mm := as.character(entrez_mm)]
cg[, entrez_hs := as.character(entrez_hs)]
clockGenes = cg

bg = fread(file.path('data-raw', 'genes_blood.csv'))
bg[, entrez_hs := as.character(entrez_hs)]
bloodGenes = bg

refCorMouseEntrez = readRDS(file.path('data-raw', 'result_mouse_ref.rds'))
refCorMouseEntrez = refCorMouseEntrez[cg$entrez_mm, cg$entrez_mm]

refCorHumanBlood = readRDS(file.path('data-raw', 'result_blood_ref.rds'))

usethis::use_data(clockGenes, bloodGenes, refCorMouseEntrez, refCorHumanBlood,
                   internal=TRUE, overwrite=TRUE)

####################

sampleMetadataGeo = fread(file.path('data-raw', 'sample_metadata_cancer.csv'))
sampleMetadataGeo = sampleMetadataGeo[condition %in% c('tumor', 'nontumor')]
sampleMetadataGeo[condition == 'nontumor', condition := 'non-tumor']
ematListGeo = extractExpressionData(readRDS(file.path('data-raw', 'proc_cancer.rds')), sampleMetadataGeo)

studyNow = 'GSE19188'
groupVec = sampleMetadataGeo$condition[sampleMetadataGeo$study==studyNow]
emat = ematListGeo[[studyNow]][,sampleMetadataGeo$sample[sampleMetadataGeo$study==studyNow]]
GSE19188 = list(emat=emat, groupVec=groupVec)

usethis::use_data(GSE19188, overwrite=TRUE)
