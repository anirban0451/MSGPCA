test_that('neighbor indices and labels match upstream including tied grids', {
  grids <- list(as.matrix(expand.grid(x=1:5,y=1:5)),
    cbind(rep(1:5,5)+rep(c(0,0.5,0,0.5,0),each=5),rep(1:5,each=5)*sqrt(3)/2),
    matrix(sin(1:50),25,2))
  for(z in grids) for(shape in c('square','hexagon')) {
    q <- if(shape=='square') 4L else 6L
    d <- as.matrix(dist(z))
    nn <- t(vapply(1:25,function(i) order(d[i,])[seq_len(q+1)],integer(q+1)))
    expect_identical(MSGPCA:::.refine_neighbors(z,q),nn)
    labels <- rep(c('alpha','beta','gamma'),length.out=25)
    expect_identical(refine_cluster(labels,z,shape),upstream(labels,z,shape))
  }
})
test_that('label identities and invalid cases are explicit', {
  z <- as.matrix(expand.grid(x=1:3,y=1:3)); a <- c(10,10,10,10,30,10,10,10,10)
  expect_equal(refine_cluster(a,z),rep('10',9))
  expect_identical(refine_cluster(factor(a,levels=c(10,30,99)),z),refine_cluster(as.character(a),z))
  expect_error(refine_cluster(c(NA,rep(1,8)),z),'nonmissing')
  expect_error(refine_cluster(a,z,'triangle'),'arg')
  expect_error(refine_cluster(1:4,z[1:4,]),'At least')
  z[2,] <- z[1,]; expect_error(refine_cluster(a,z),'Duplicate')
})
test_that('strict rule differs from generic majority and is sequential', {
  z <- cbind(1:8,0); a <- c('a','b','b','a','a','b','b','a')
  expected <- upstream(a,z)
  expect_identical(refine_cluster(a,z),expected)
  # Exhaust binary labels to ensure fixtures discriminate update semantics.
  nn <- MSGPCA:::.refine_neighbors(z,4L); sequential_diff <- strict_diff <- FALSE
  for(k in 0:255) {
    a <- ifelse(as.integer(intToBits(k))[1:8]==1,'a','b')
    seq <- refine_cluster(a,z)
    simultaneous <- vapply(1:8,function(i) {
      tab <- table(a[nn[i,]])
      if(tab[a[i]]<2 && max(tab)>2) names(tab)[which.max(tab)] else a[i]
    },'')
    generic <- vapply(1:8,function(i) names(which.max(table(a[nn[i,]]))),'')
    sequential_diff <- sequential_diff || !identical(seq,simultaneous)
    strict_diff <- strict_diff || !identical(seq,generic)
  }
  expect_true(sequential_diff); expect_true(strict_diff)
})
test_that('accelerated and rowwise engines return exactly identical labels', {
  set.seed(72)
  for(z in list(as.matrix(expand.grid(x=1:8,y=1:8)),matrix(runif(128),64,2))) {
    for(shape in c('square','hexagon')) {
      lab <- sample(c('2','10','99'),nrow(z),TRUE)
      expect_identical(refine_cluster(lab,z,shape,engine='knn'),refine_cluster(lab,z,shape,engine='reference'))
    }
  }
})
