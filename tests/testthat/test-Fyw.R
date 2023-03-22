test_that("multiple function works", {
  # Test regression coefficient
  testValue <- IRLS(Y = c(1:10),X = Y+ runif(10)*2)
  expect_true(testValue$coefficients[2] > 0)
  testValue <- IRLS(Y = c(1:10),X = -Y  - runif(10)*2)
  expect_true(testValue$coefficients[2] < 0)

  # Test find Fyw
  testValue <- findFyw(AP = 3, phiP = pi, AS = 2, phiS = pi*2)
  expect_equal(testValue$Fyw, 0.667)
  expect_equal(testValue$alpha, 12.305)
  expect_equal(testValue$beta, 0.042)
  expect_equal(testValue$meanTT_year, 0.511)
})
