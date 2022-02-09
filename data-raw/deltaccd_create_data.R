library('data.table')

parentDir = 'data-raw'

####################

cg = fread(file.path(parentDir, 'genes_clock.csv'))
cg[, entrez_mm := as.character(entrez_mm)]
cg[, entrez_hs := as.character(entrez_hs)]
clockGenes = cg

bg = fread(file.path(parentDir, 'genes_blood.csv'))
bg[, entrez_hs := as.character(entrez_hs)]
bloodGenes = bg

refCorMouseEntrez = readRDS(file.path(parentDir, 'result_mouse_ref.rds'))
refCorMouseEntrez = refCorMouseEntrez[cg$entrez_mm, cg$entrez_mm]

refCorHumanBlood = readRDS(file.path(parentDir, 'result_blood_ref.rds'))

usethis::use_data(clockGenes, bloodGenes, refCorMouseEntrez, refCorHumanBlood,
                  internal = TRUE, overwrite = TRUE)

####################

study = 'GSE19188'
seeker::seekerArray(list(study = study, geneIdType = 'entrez'), parentDir)

metadata = fread(file.path(parentDir, study, 'sample_metadata.csv'))
groupVec = metadata[['tissue type:ch1']]

emat = qs::qread(file.path(parentDir, study, 'gene_expression_matrix.qs'))
idx = c(which(rownames(emat) %in% clockGenes$entrez_hs),
        seq(1, nrow(emat), 30))
emat = emat[unique(idx), ]

GSE19188 = list(emat = emat, groupVec = groupVec)

usethis::use_data(GSE19188, overwrite = TRUE)

unlink(file.path(parentDir, study), recursive = TRUE)
