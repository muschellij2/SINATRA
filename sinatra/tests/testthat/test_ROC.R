context("tests ROC generation")

test_that("ROC Generation Works", {
expect_error(generate_ROC_with_coned_directions(nsim = 10, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1,
                                   causal_points = 5,shared_points = 5,num_cones = 5, eta = 0.1,
                                   truncated = 250, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
                                   min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT', mode = 'sphere',
                                   subdivision = 3,num_causal_region = 3, num_shared_region = 6), NA)
})
