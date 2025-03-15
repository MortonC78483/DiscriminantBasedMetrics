n = 10

# simluate draws from normal distribution
conf_data = rnorm(n)
resampled = sample(conf_data, 2*n, replace = TRUE)
assigned_conf_status = c(rep(0, 10), rep(1, 10))

# for each data point, how many times was it sampled
s
