test_that('getRefCor', {
  expect_error(getRefCor(species = 'mouse', tissue = 'blood'),
               'Blood reference is only available for species = \'human\'.')
  
  expect_equal(getRefCor(species = 'mouse', tissue = 'pan'), refCorMouseEntrez)
  
  expect_equal(getRefCor(species = 'human', tissue = 'pan'), refCorHumanBlood)
})

# test_that('calcCCD', {
#   
# })
# 
# test_that('calcDeltaCCD', {
#   
# })
# 
# test_that('plotHeatmap', {
#   
# })
# 
# test_that('plotRefHeatmap', {
#   
# })