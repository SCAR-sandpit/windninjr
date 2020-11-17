context("Config checking")

test_that("config works in general", {
    demfile <- wn_demo_file("missoula_valley_elevation")
    my_config <- my_config <- wn_config_domain_average(elevation = demfile, input_speed = 10, input_direction = 45)
})
