y <- 1:300000
elim <- factor(rep(1:1000, 300))

missing <- sample(1:300000, 300)
elim <- elim[-missing]
y <- y[-missing]

#1
system.time(tapply(y, elim, mean))

#3 save quite a bit adding mean.default
system.time(sapply(split(y, elim), mean.default))

#4 - a bit of saving over sapply - now 1/2 time of tapply
# - also about 1/4 time of llbt!

system.time(as.vector(lapply(split(y, elim), mean.default)))

#5 use sum instead? much faster
size <- table(elim)
system.time(sapply(split(y, elim), sum)/size)
system.time(unlist(lapply(split(y, elim), sum))/si)

#6 use cumsum? - superfast!
# more comparable to #5 if many categories with large size
# still faster with 1000 cat, size approx 300
y.ord <-as.numeric(y[order(elim)])
size <- table(elim)
end <- cumsum(size)

system.time({
    cy <- cumsum(y.ord)[end]
    res <- diff(c(0,cy))/size
})

grp.mean <- function(x, grp.end, grp.size){
    diff(c(0,cumsum(x)[grp.end]))/size
}

########### best answer ###########
grp.sum <- function(x, grp.end){
    diff(c(0,cumsum(x)[grp.end]))
}

system.time(res2 <- grp.sum(y.ord, end)/size)
###################################

## slower - equiv to #4
elim.ord <- sort.int(elim)

system.time(res3 <- drop(rowsum(y.ord, elim.ord))/size)

#7 colSums instead - test idea
# about the same as cumsum, but more overhead dealing with NAs etc
num <- sequence(size)
y.ord <- c(xtabs(y ~ elim + num)) #could be sped up

system.time({res2 <- rowSums(matrix(y.ord, nc = 300))/size})



#5 - save by only doing split once!
ind <- split(seq_along(elim), elim)
lev <- levels(elim)

res <- numeric(length(lev))
names(res) <- lev

## takes ages for many levels
system.time(for (k in lev) res[k] <- mean(y[ind[[k]]]))

## ah! if x is classless, uses internal split, not this method:(
res <- vector("list", length(lev))
names(res) <- lev
system.time({for (k in lev) res[[k]] <- y[ind[[k]]]; unlist(lapply(res,mean))})

#############

ngrp <- 1000
grpsize <- 50
nx <- 50
x <- matrix(1:(ngrp*grpsize*nx),nc= nx)
keep.x <- x
eliminate <- factor(rep(1:ngrp, grpsize))

missing <- sample(1:(ngrp*grpsize), grpsize)
eliminate <- eliminate[-missing]
x <- x[-missing,]
keep.x <- x

x <- keep.x
system.time({
subtracted <- matrix(0, ncol(x), nlevels(eliminate))
for (j in 2:ncol(x)) {  ## sweeps needed to get the rank right
    colj <- x[, j]
    subtracted[j, ] <- sj <- sapply(split(colj, eliminate), mean)
    x[, j] <- colj - sj[eliminate]
}
})

x <- keep.x
## best idea from above
x.ord <-x[order(eliminate),]
elim.ord <- sort.int(eliminate)
size <- as.vector(table(eliminate))
end <- cumsum(size)

quickdiff <- function(x) {x[-1] - x[-length(x)]}

grp.sum <- function(x, grp.end){
    quickdiff(c(0, cumsum(as.numeric(x))[grp.end]))
}

system.time({
 subtracted <- matrix(0, ncol(x), nlevels(eliminate))
for (j in 2:ncol(x)) {  ## sweeps needed to get the rank right
    cj <- x.ord[,j]
    subtracted[j, ] <- sj <- grp.sum(cj, end)/size
    x.ord[, j] <- cj - sj[elim.ord]
}
})

#slower
system.time({
 subtracted <- matrix(0, ncol(x), nlevels(eliminate))
for (j in 2:ncol(x)) {  ## sweeps needed to get the rank right
    cj <- x.ord[,j]
    subtracted[j, ] <- sj <- unlist(lapply(split(cj, eliminate), sum))/size
    x.ord[, j] <- cj - sj[elim.ord]
}
})

x <- keep.x
elim.ord <- sort.int(eliminate)
x.ord <- x[order(eliminate),]
size <- as.vector(table(eliminate))
end <- cumsum(size)

nelim <- nlevels(eliminate)

endcols <- seq.int(nrow(x), length(x) - nrow(x), nrow(x))
endcols <- matrix(endcols, ncol(x) - 1, length(end))
endx <- c(t(endcols) + end)

#endx <- c(outer(end, endcols, "+"))

system.time({
    endcols <- nrow(x)*seq(ncol(x) - 1)
    endx <- unlist(lapply(endcols, "+", end))
    endx <- c(nrow(x), endx)
    subtracted <- matrix(c(numeric(nelim), quickdiff(cumsum(as.numeric(x.ord))[endx])/size),
                         nlevels(eliminate))
    x.ord <- x.ord - subtracted[elim.ord,]
})

############ best time #################
x <- keep.x
elim.ord <- sort.int(eliminate)
x.ord <- x[order(eliminate),]
size <- as.vector(table(eliminate))

system.time({
    subtracted <- rowsum(x.ord, elim.ord)/size
    subtracted[,1] <- 0
    x.ord <- x.ord - subtracted[elim.ord,]
})

## same
x <- keep.x
elim.ord <- sort.int(eliminate)
x.ord <- x[order(eliminate),]
size <- as.vector(table(eliminate))

system.time({
    subtracted <- rowsum(x.ord, elim.ord)/size
    subtracted[,1] <- numeric(nelim)
    x.ord <- x.ord - subtracted[elim.ord,]
})

## longer
x <- keep.x
elim.ord <- sort.int(eliminate)
x.ord <- x[order(eliminate),]
size <- as.vector(table(eliminate))

system.time({
    subtracted <- cbind(numeric(nelim), rowsum(x.ord[,-1], elim.ord)/size)
    x.ord <- x.ord - subtracted[elim.ord,]
})

## ages
x <- keep.x
elim.ord <- sort.int(eliminate)
x.ord <- x[order(eliminate),]
size <- as.vector(table(eliminate))

system.time({
    subtracted <- rowsum(x.ord[,-1], elim.ord)/size
    x.ord[,-1] <- x.ord[,-1] - subtracted[elim.ord,]
})



### good for low number of levels!
ind <- split(seq_along(eliminate), eliminate)
lev <- levels(eliminate)
system.time({
subtracted <- matrix(0, nlevels(eliminate), ncol(x) - 1)
rownames(subtracted) <- lev
for (k in lev) subtracted[k,] <- colMeans(x[ind[[k]], -1, drop =FALSE])
x[,-1] <- x[,-1] - subtracted[eliminate,]
})

## BUT expecting lev to be much higher than ncol(x) - still better to use colMeans approach

x <- keep.x

system.time({
subtracted <- matrix(0, ncol(x), nlevels(eliminate))
colnames(subtracted) <- lev
for (j in 2:ncol(x)) {  ## sweeps needed to get the rank right
    for (k in lev) subtracted[j, k] <- mean(x[ind[[k]], j])
    x[, j] <- x[, j] - subtracted[j, eliminate]
}
})
