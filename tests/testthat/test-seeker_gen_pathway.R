test_that("seeker gen pathway format", {
  expect_equal(class(seeker_gen_pathway("MAPT")), "data.frame")
  df <- data.frame(gen=c("MAPT", "APOE", "MMP12"))
  expect_equal(class(seeker_gen_pathway(df)), "data.frame")
})


