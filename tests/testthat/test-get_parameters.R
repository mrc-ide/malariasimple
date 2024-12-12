test_that("n_days must be > 0", {
  expect_error(get_parameters(n_days = -1))
})

##age_rendering
test_that("Default values give correct age_rendering outputs",{
  params <- list(age_vector = c(0,0.25,0.5,1,1.5,2,2.5,3,3.5,4,5,6,7,8.5,10,20,30,40,60,80)*365)
  age_render <- age_rendering(params = params,
                              prevalence_rendering_max_ages = NULL,
                              prevalence_rendering_min_ages = NULL,
                              clin_inc_rendering_max_ages = NULL,
                              clin_inc_rendering_min_ages = NULL)
  expect_equal(age_render$min_age_prev, 6)
  expect_equal(age_render$max_age_prev, 14)
  expect_equal(age_render$min_age_inc, 1)
  expect_equal(age_render$max_age_inc, 20)
})

