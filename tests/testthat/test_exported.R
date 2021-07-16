test_that('getRefCor', {
  expect_error(getRefCor(species = 'mouse', tissue = 'blood'),
               'Blood reference is only available for species = \'human\'.')
  
  expect_equal(getRefCor(species = 'mouse', tissue = 'pan'), refCorMouseEntrez)
  
  expect_equal(getRefCor(species = 'human', tissue = 'blood'), refCorHumanBlood)
})

test_that('calcCCD', {
  ref = diag(1, 2)
  rownames(ref) = paste0('gene_', 1:2)
  colnames(ref) = rownames(ref)
  
  emat = rbind(c(1, 0, 1), c(0, 1, 0))
  rownames(emat) = paste0('gene_', 1:2)
  
  ccd = calcCCD(ref, emat)
  expect_equal(ccd, data.table(group = 'all', CCD = 1, Pvalue = 0.5))
  
  groupVec = c('a', 'b', 'c', 'd')
  expect_error(calcCCD(ref, emat, groupVec), 
               'Length of groupVec does not match the number of columns in emat.')
  
  groupVec = c('a', 'b', 'c')
  expect_error(calcCCD(ref, emat, groupVec), 
               'Each unique group in groupVec must have at least three samples.')
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
               'The supplied value for groupNormal is not present in groupVec.')
  
  groupVec = rep('a', 6)
  expect_error(calcDeltaCCD(refCor, emat, groupVec, groupNormal = 'a'),
               'groupVec contains only one unique group.')
  
  groupVec = c(rep('a', 4), 'b', 'b')
  expect_error(calcDeltaCCD(refCor, emat, groupVec, groupNormal = 'a'),
               'Each unique group in groupVec must have at least three samples.')
  
  groupVec = c(groupVec, 'b')
  expect_error(calcDeltaCCD(refCor, emat, groupVec, groupNormal = 'a'),
               'Length of groupVec does not match the number of columns in emat.')
})

test_that('plotHeatmap', {
  ematNow = rbind(c(1, 0.9, 0.8, 1, 0.9, 0.8), c(1, 0.9, 0.8, -1, -0.9, -0.8))
  rownames(ematNow) = paste0('gene_', 1:2)
  groupVec = c(rep('a', 3), rep('b', 3))
  geneNames = paste0('gene_', 1:2)
  p = plotHeatmap(geneNames, ematNow, groupVec)
  
  vdiffr::expect_doppelganger('basic heatmap', p)
  
  groupVec = c(rep('a', 3), rep('b', 4))
  expect_error(plotHeatmap(geneNames, ematNow, groupVec),
               'Length of groupVec does not match the number of columns in emat.')  
  
  groupVec = c('a', 'a', 'b', 'b', 'c', 'c')
  expect_error(plotHeatmap(geneNames, ematNow, groupVec),
               'Each unique group in groupVec must have at least three samples.')

  groupVec = c(rep('a', 3), rep('b', 3))
  geneNames = paste0('gene_', c(1, 3))
  expect_error(plotHeatmap(geneNames, ematNow, groupVec),
               'Fewer than two genes in the supplied vector are in the expression matrix.')

  # geneNames = paste0('gene_', 1:3)
  # expect_warning(plotHeatmap(geneNames, ematNow, groupVec),
  #              '1 gene(s) in the supplied vector is/are not in the expression matrix.')
})

test_that('plotRefHeatmap', {
  ref = diag(1, 2)
  rownames(ref) = paste0('gene_', 1:2)
  colnames(ref) = rownames(ref)
  p = plotRefHeatmap(ref)
  vdiffr::expect_doppelganger('ref heatmap', p)
  
  ref = cbind(ref, c(1, 1))
  colnames(ref) = paste0('gene_', 1:3)
  expect_error(plotRefHeatmap(ref), 
               'refCor must be a correlation matrix, with identical rownames and colnames.')
})