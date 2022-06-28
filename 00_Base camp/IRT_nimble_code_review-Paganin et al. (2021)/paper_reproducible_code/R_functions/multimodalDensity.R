dMultiModal <- function(x, weights = c(1,1,1), means = c(-2, 0, 3)){
    require(sn)

    prop <- weights/sum(weights)

    prop[1]*dnorm(x, mean = means[1], sd = sqrt(1)) + 
    prop[2]*dnorm(x, mean = means[2], sd = sqrt(0.5)) + 
    prop[3]*dsn(x, xi = means[3], omega = 1, alpha = -3)
}

## curve(dMultiModal(x, weights = c(4,4,2)), from = -5, to = 5)

rMultiModal <- function(n, weights = c(4,4,2), means = c(-2, 0, 3)){
    prop <- weights/sum(weights)

    out <- c(rnorm(n*prop[1], mean = means[1], sd = sqrt(1)),
    rnorm(n*prop[2], mean = means[2], sd = sqrt(0.5)),
    rsn(n*prop[3], xi = means[3], omega = 1, alpha = -3))

    out
}
