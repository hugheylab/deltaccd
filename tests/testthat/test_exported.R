library('ggplot2')
library('vdiffr')

test_that('getRefCor', {
  expect_error(getRefCor(species = 'mouse', tissue = 'blood'),
               'Blood reference is only available for species = \'human\'.',
               fixed = TRUE)
  
  expect_equal(getRefCor(species = 'mouse', tissue = 'pan'), refCorMouseEntrez)
  
  expect_equal(getRefCor(species = 'human', tissue = 'blood'), refCorHumanBlood)
})

test_that('calcCCD', {
  ref = cbind(c(1, 1), c(1, 1))
  rownames(ref) = paste0('gene_', 1:2)
  colnames(ref) = rownames(ref)
  
  emat = rbind(c(1, 0, 1), c(0, 1, 0))
  rownames(emat) = paste0('gene_', 1:2)
  
  ccd = calcCCD(ref, emat)
  expect_equal(ccd, data.table(group = 'all', CCD = 2, Pvalue = 0.5))
  
  ccd = calcCCD(ref, emat, nPerm = 0)
  expect_equal(ccd, data.table(group = 'all', CCD = 2, Pvalue = NA))
  
  groupVec = c('a', 'b', 'c', 'd')
  expect_error(calcCCD(ref, emat, groupVec), 
               'Length of groupVec does not match the number of columns in emat.',
               fixed = TRUE)
  
  groupVec = c('a', 'b', 'c')
  expect_error(calcCCD(ref, emat, groupVec), 
               'Each unique group in groupVec must have at least three samples.',
               fixed = TRUE)
})

test_that('calcDeltaCCD', {
  refCor = rbind(c(1, 1), c(1, 1))
  rownames(refCor) = paste0('gene_', 1:2)
  colnames(refCor) = rownames(refCor)
  
  emat = rbind(c(1, 0.9, 0.8, 1, 0.9, 0.8), c(1, 0.9, 0.8, 1, 0.9, 0.8))
  rownames(emat) = paste0('gene_', 1:2)
  groupVec = c(rep('a', 3), rep('b', 3))
  
  dccd = calcDeltaCCD(refCor, emat, groupVec, groupNormal = 'a')
  expect_equal(dccd$DeltaCCD, 0)
  
  expect_error(calcDeltaCCD(refCor, emat, groupVec, groupNormal = 'j'),
               'The supplied value for groupNormal is not present in groupVec.',
               fixed = TRUE)
  
  groupVec = rep('a', 6)
  expect_error(calcDeltaCCD(refCor, emat, groupVec, groupNormal = 'a'),
               'groupVec contains only one unique group.', fixed = TRUE)
  
  groupVec = c(rep('a', 4), 'b', 'b')
  expect_error(calcDeltaCCD(refCor, emat, groupVec, groupNormal = 'a'),
               'Each unique group in groupVec must have at least three samples.',
               fixed = TRUE)
  
  groupVec = c(groupVec, 'b')
  expect_error(calcDeltaCCD(refCor, emat, groupVec, groupNormal = 'a'),
               'Length of groupVec does not match the number of columns in emat.',
               fixed = TRUE)
})

test_that('plotHeatmap', {
  ematNow = rbind(c(1, 0.9, 0.8, 1, 0.9, 0.8), c(1, 0.9, 0.8, -1, -0.9, -0.8))
  rownames(ematNow) = paste0('gene_', 1:2)
  groupVec = c(rep('a', 3), rep('b', 3))
  geneNames = paste0('gene_', 1:2)
  p = plotHeatmap(geneNames, ematNow, groupVec)
  
  expect_doppelganger('basic heatmap', p)
  
  groupVec = c(rep('a', 3), rep('b', 4))
  expect_error(plotHeatmap(geneNames, ematNow, groupVec),
               'Length of groupVec does not match the number of columns in emat.',
               fixed = TRUE)  
  
  groupVec = c('a', 'a', 'b', 'b', 'c', 'c')
  expect_error(plotHeatmap(geneNames, ematNow, groupVec),
               'Each unique group in groupVec must have at least three samples.',
               fixed = TRUE)

  groupVec = c(rep('a', 3), rep('b', 3))
  geneNames = paste0('gene_', c(1, 3))
  expect_error(plotHeatmap(geneNames, ematNow, groupVec),
               'Fewer than two genes in the supplied vector are in the expression matrix.',
               fixed = TRUE)

  geneNames = paste0('gene_', 1:3)
  expect_warning(plotHeatmap(geneNames, ematNow, groupVec),
               '1 gene(s) in the supplied vector is/are not in the expression matrix.',
               fixed = TRUE)
})

test_that('plotRefHeatmap', {
  ref = cbind(c(1, 0.5, -1), c(0.5, 1, -0.5), c(-1, -0.5, 1)) 
  rownames(ref) = paste0('gene_', 1:3)
  colnames(ref) = rownames(ref)
  p = plotRefHeatmap(ref)
  expect_doppelganger('ref heatmap', p)
  
  ref = cbind(ref, c(1, 1, 1))
  colnames(ref) = paste0('gene_', 1:4)
  expect_error(suppressWarnings(plotRefHeatmap(ref)), 
               'refCor must be a correlation matrix, with identical rownames and colnames.',
               fixed = TRUE)
})