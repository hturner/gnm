gnmStart <- function(n, scale = 0.1) {
    theta <- runif(n, -1, 1) * scale
    theta + (-1)^(theta < 0) * scale
}
