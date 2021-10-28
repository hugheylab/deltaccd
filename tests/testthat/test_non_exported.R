library('ggplot2')
library('vdiffr')

geneNames = paste0('gene_', 1:2)

genCors = function() {
  ematNow = rbind(c(1, 0.9, 0.8, 1, 0.9, 0.8), c(1, 0.9, 0.8, -1, -0.9, -0.8))
  rownames(ematNow) = geneNames
  groupVec = c('a' ,'a', 'a', 'b', 'b', 'b')
  cors = calcCorr(ematNow, groupVec)
  return(cors)}


test_that('calcDist', {
  d = calcDist(c(1, 0), c(1, 0))
  expect_equal(d, 0)
  
  d = calcDist(c(1, 0), c(0, 1))
  expect_equal(d, sqrt(2))
})

test_that('calcCCDSimple', {
  ccd = calcCCDSimple(diag(1, 2), diag(1, 2))
  expect_equal(ccd, 1)
  
  ref = cbind(c(1, 1), c(1, 1))
  emat = rbind(c(1, 0, 1), c(0, 1, 0))
  ccd = calcCCDSimple(ref, emat)
  expect_equal(ccd, 2)
  
  ref = cbind(c(1, 0, -1), c(0, 1, -1), c(-1, -1, 1)) 
  emat = rbind(c(1, 0, 0), c(0, 1, -1), c(0, 1, 1))
  ccd = calcCCDSimple(ref, emat, scale = TRUE)
  expect_equal(ccd, 1/3)
})

test_that('checkVar', {
  groupVec = c('a', 'a', 'b', 'b')
  emat = cbind(diag(1, 2), diag(1, 2))
  expect_invisible(checkVar(emat, groupVec))
  
  emat = cbind(diag(1, 2), matrix(1, 2, 2))
  rownames(emat) = paste0('gene_', 1:2)
  expect_error(
    checkVar(emat, c('a', 'a', 'b', 'b')), 
    paste0('Zero variance in 2 gene-group pairs. ',
           'See _deltaccd_zero_var_genes.csv for full list.'), fixed = TRUE)
})

test_that('checkGenes', {
  emat = diag(1, 2)
  rownames(emat) = geneNames 
  refCor = emat 
  rownames(refCor) = geneNames
  expect_equal(checkGenes(emat, refCor), geneNames) 
  
  refCor = diag(1, 4)
  rownames(refCor) = paste0('gene_', 1:4)
  expect_error(
    checkGenes(emat, refCor), 
    paste0('2 gene(s) is/are not in the expression matrix. ', 
           'See _deltaccd_missing_genes.csv for full list.'), fixed = TRUE)
  
  emat = diag(1, 6)
  refCor = diag(1, 4)
  geneNames2 = paste0('gene_', 1:6)
  rownames(emat) = geneNames2
  rownames(refCor) = geneNames2[1:4]
  expect_equal(checkGenes(emat, refCor, geneNames2[1:5]), geneNames2[1:4])
})

test_that('checkRefCor', {
  refCor = diag(1, 2)
  rownames(refCor) = geneNames
  colnames(refCor) = rownames(refCor)
  expect_equal(checkRefCor(refCor), refCor)
  
  expect_error(checkRefCor(), 'Either refCor or refEmat must be supplied.',fixed = TRUE)
  
  rownames(refCor) = paste0('gene_', 3:4)
  expect_error(checkRefCor(refCor),
               'refCor must be a correlation matrix, with identical rownames and colnames.',
               fixed = TRUE)
  
  refCor = rbind(c(1, 0, 0), c(1, 0, 0))
  expect_error(checkRefCor(refCor),
               'refCor must be a correlation matrix, with identical rownames and colnames.',
               fixed = TRUE)
})

test_that('calcDeltaCCDSimple', {
  ref = cbind(c(1, -1), c(-1, 1)) 
  emat = ref
  expect_equal(calcDeltaCCDSimple(ref, emat, 1:2), 0)
})

test_that('makePerms', {
  perms = makePerms(1:10)
  expect_equal(dim(perms), c(1000, 10))
})

test_that('calcCorr', {
  cors = genCors()
  
  expect_equal(cors, data.table(group = c('a', 'b', 'a', 'b'), 
                                gene1 = factor(paste0('gene_', c(2, 2, 1, 1))), 
                                gene2 = factor(paste0('gene_', c(1, 1, 2, 2)),
                                               levels = paste0('gene_', 2:1)),
                                rho = c(1, -1, 1, -1)))
})

test_that('calcColorLimits', {
  cors = genCors()
  cLims = calcColorLimits(cors$rho)

  expect_equal(cLims, c('#E66101', '#5E3C99'))
})

test_that('plotHeatmapSimple', {
  cors = genCors()
  cLims = calcColorLimits(cors$rho)

  p = plotHeatmapSimple(
    ggplot(cors) + facet_wrap(vars(group)), cLims)
  expect_doppelganger('basic heatmap', p)
})
