laplacian.eigenmap <- function(W, d = 2, normalized = TRUE) {
  n <- nrow(W)
  if (normalized) {
    L <- normalized.laplacian(W)
  } else {
    L <- graph.laplacian(W)
  }
  L.eigen <- eigen(L, symmetric = TRUE)
  sweep(L.eigen$vectors[, seq(n - 1, n - d)], 2, 
        sqrt(L.eigen$values[seq(n - 1, n - d)]), `/`)
}

graph.laplacian <- function(W) {
  return(diag(colSums(W)) - W)
}

spectral.clustering <- function(W, K = 2, d = K, normalized = TRUE) {
  Y <- laplacian.eigenmap(W, d, normalized)
  kmeans(Y, K)$cluster
}

estimate.one.t.bezier <- function(x, p, 
                                  min.t = 0, max.t = 2, 
                                  intercept = TRUE, 
                                  tol = 1e-9,
                                  method = 'polyroot') {
  polynomial <- function(t, p) {
    if (intercept) {
      r <- nrow(p) - 1
      T <- sapply(seq(0, r), function(s) {
        choose(r, s) * (1 - t) ^ (r - s) * t ^ s
      })
    } else {
      r <- nrow(p)
      T <- sapply(seq(r), function(s) {
        choose(r, s) * (1 - t) ^ (r - s) * t ^ s
      })
    }
    as.numeric(T %*% p)
  }
  l <- function(t) log(as.numeric(crossprod(polynomial(t, p) - x)))
  if (method == 'optimize') {
    minima <- optimize(l, c(min.t, max.t), tol = tol)$minimum
  } else if (method == 'optim') {
    minima <- optim((min.t + max.t) / 2, l, method = 'SANN')$par
  } else if (method == 'polyroot') {
    if (!intercept) {
      p <- rbind(0, p)
    }
    p0 <- p[1, ]
    p1 <- p[2, ]
    p2 <- p[3, ]
    b0 <- sum((x - p0) * (p1 - p0))
    b1 <- sum((x - p0) * (p0 - 2 * p1 + p2) - 2 * (p1 - p0) ^ 2)
    b2 <- -3 * sum((p1 - p0) * (p0 - 2 * p1 + p2))
    b3 <- -sum((p0 - 2 * p1 + p2) ^ 2)
    roots <- polyroot(c(b0, b1, b2, b3))
    roots <- Re(roots[abs(Im(roots)) < tol])
    if (length(roots) > 1) {
      minima <- roots[which.min(sapply(roots, l))]
    } else if (length(roots) == 1) {
      minima <- roots
    } else if (length(roots) == 0) {
      warning('failed to find a real root')
      minima <- optimize(l, c(min.t, max.t), tol = tol)$minimum
    }
  } else {
    stop('invalid method for finding t')
  }
  # if (abs(minima) == max.t) warning(minima)
  # if (length(minima) > 1) warning(minima)
  return(max(minima))
}

estimate.t.bezier <- function(X, p, t.prev, 
                              intercept = TRUE, 
                              min.t = 0, max.t = 2, 
                              tol = 1e-9,
                              parallel = FALSE) {
  if (missing(t.prev)) {
    t.hat <- plyr::laply(seq_len(nrow(X)), function(i) {
      estimate.one.t.bezier(X[i, ], p, 
                            min.t = min.t, max.t = max.t, 
                            tol = tol, intercept = intercept)
    }, .parallel = parallel)
  } else {
    t.hat <- plyr::laply(seq_along(t.prev), function(i) {
      estimate.one.t.bezier(X[i, ], p, 
                            t.prev[i] - 2, t.prev[i] + 2, 
                            tol = tol, intercept = intercept)
    }, .parallel = parallel)
    loss.new <- sapply(seq_along(t.hat),
                       function(i) loss.one.t(X[i, ], t.hat[i], p))
    loss.prev <- sapply(seq_along(t.prev),
                        function(i) loss.one.t(X[i, ], t.prev[i], p))
    prop.change <- (mean(loss.new < loss.prev))
    if (prop.change < 1) {
      warning(paste('Failed to find new t for', 1 - prop.change))
    }
    t.hat <- ifelse(loss.new <= loss.prev, t.hat, t.prev)
  }
  
  return(t.hat)
}

construct.bezier.model.matrix <- function(t.hat, r = 2, intercept = TRUE) {
  if (intercept) {
    sapply(seq(0, r), function(s) {
      choose(r, s) * (1 - t.hat) ^ (r - s) * t.hat ^ s
    })
  } else {
    sapply(seq(r), function(s) {
      choose(r, s) * (1 - t.hat) ^ (r - s) * t.hat ^ s
    })
  }
}

bezier.curve <- function(t, p, intercept = TRUE) {
  if (intercept) {
    r <- nrow(p) - 1
  } else {
    r <- nrow(p)
  }
  T <- construct.bezier.model.matrix(t, r, intercept = intercept)
  return(T %*% p)
}

normalize.umvue <- function(t.hat) {
  n <- length(t.hat)
  w <- n ^ 2 - 1
  min.t <- min(t.hat)
  max.t <- max(t.hat)
  a <- (n * (n + 1) * min.t - (n + 1) * max.t) / w
  b <- (-(n + 1) * min.t + n * (n + 1) * max.t) / w
  return((t.hat - a) / (b - a))
}

normalize.ecdf <- function(t.hat) {
  ecdf(t.hat)(t.hat)
}

estimate.bezier.curve.2 <- function(X, 
                                    init.params, 
                                    weights,
                                    eps = 1e-3, 
                                    maxit = 200,
                                    initialization = 'isomap', 
                                    normalize = TRUE, 
                                    normalize.method = 'umvue',
                                    k.isomap = as.integer(sqrt(nrow(X))),
                                    intercept = TRUE,
                                    parallel = FALSE) {
  n <- nrow(X)
  d <- ncol(X)
  r <- 2
  
  if (missing(weights)) {
    weights <- rep(1, n)
  }
  W <- diag(weights)
  
  if (missing(init.params)) {
    if (initialization == 'random') {
      if (intercept) {
        p <- matrix(rnorm(d * (r + 1)), nrow = r + 1, ncol = d)
      } else {
        p <- matrix(rnorm(d * r), nrow = r, ncol = d)
      }
      t.hat <- estimate.t.bezier(X, p, intercept = intercept, parallel = parallel)
      if (!intercept) {
        if (sum((c(0, 0) %*% p) ^ 2) > sum((c(1, 1) %*% p) ^ 2)) {
        # if (sum(X[which.min(t.hat), ] ^ 2) > sum(X[which.max(t.hat), ] ^ 2)) {
          t.hat <- 1 - t.hat
        }
      }
      
    } else if (initialization == 'isomap') {
      isomap.out <- estimate.bezier.curve.isomap(X, k.isomap, weights, 
                                                 intercept = intercept)
      p <- isomap.out$p
      t.hat <- isomap.out$t
    } else {
      stop('initialization must be random or isomap')
    }
  } else {
    p <- init.params
    t.hat <- estimate.t.bezier(X, p, intercept = intercept, parallel = parallel)
    # if (!intercept) {
    #   if (sum(X[which.min(t.hat), ] ^ 2) > sum(X[which.max(t.hat), ] ^ 2)) {
    #     t.hat <- 1 - t.hat
    #     T <- construct.bezier.model.matrix(t.hat, r, intercept = intercept)
    #     # p <- solve(t(T) %*% W %*% T, t(T) %*% W %*% X)
    #     # t.hat <- estimate.t.bezier(X, p, intercept = intercept, parallel = parallel)
    #   }
    # }
  }
  
  if (normalize) {
    if (normalize.method == 'umvue') {
      t.hat <- normalize.umvue(t.hat)
    } else if (normalize.method == 'ecdf') {
      t.hat <- normalize.ecdf(t.hat)
    } else {
      stop('normalize.method must be umvue or ecdf')
    }
  } else {
    # t.hat <- t.hat + min(t.hat)
  }
  
  mse <- bezier.mse(X, t.hat, p, intercept)
  
  mse.b <- c()
  
  niter <- 0
  while (TRUE) {
    mse.prev <- mse

    # if (!intercept) {
    #   if (sum((c(0, 0) %*% p) ^ 2) > sum((c(1, 1) %*% p) ^ 2)) {
    #   # if (sum(X[which.min(t.hat), ] ^ 2) > sum(X[which.max(t.hat), ] ^ 2)) {
    #     t.hat <- 1 - t.hat
    #   }
    # }
    
    T <- construct.bezier.model.matrix(t.hat, r, intercept = intercept)
    p <- solve(t(T) %*% W %*% T, t(T) %*% W %*% X)
    
    # mse <- bezier.mse(X, t.hat, p, intercept = intercept)
    mse <- norm(X - T %*% p, type = 'F') ^ 2
    mse.b <- c(mse.b, mse)
    
    t.hat <- estimate.t.bezier(X, p, t.hat, 
                               intercept = intercept, 
                               parallel = parallel)
    
    mse <- bezier.mse(X, t.hat, p, intercept = intercept)
    mse.b <- c(mse.b, mse)
    
    d.mse <- (mse.prev - mse) / mse.prev
    if (d.mse < eps) break
    
    if (normalize) {
      if (normalize.method == 'umvue') {
        t.hat <- normalize.umvue(t.hat)
      } else if (normalize.method == 'ecdf') {
        t.hat <- normalize.ecdf(t.hat)
      } else {
        stop('normalize.method must be umvue or ecdf')
      }
    }
    
    niter <- niter + 1
    if (niter >= maxit) {
      warning('failed to converge to a bezier curve')
      break
    }
  }
  
  X.hat <- T %*% p
  
  return(list(p = p, 
              X = X.hat,
              t = t.hat,
              mse = mse.b))
}

knn.graph <- function(X, k) {
  dist.matrix <- mds.edm1(X)
  knn.graph <- graph.knn(dist.matrix, k) * dist.matrix
  dist.matrix[dist.matrix == 0] <- Inf
  diag(dist.matrix) <- 0
  return(graph.short(dist.matrix))
}

estimate.bezier.curve.isomap <- function(X, k = as.integer(sqrt(nrow(X))), 
                                         weights, 
                                         intercept = TRUE,
                                         normalize = TRUE) {
  k <- as.integer(k)
  n <- nrow(X)
  d <- 2
  r <- 2
  
  if (missing(weights)) {
    W <- diag(n)
  } else {
    W <- diag(weights)
  }
  D <- knn.graph(X, k)
  
  t.hat <- as.vector(cmdscale(sqrt(D), 1))
  if (!intercept) {
    if (sum(X[which.min(t.hat), ] ^ 2) > sum(X[which.max(t.hat), ] ^ 2)) {
      t.hat <- 1 - t.hat
    }
  }
  if (normalize) t.hat <- ecdf(t.hat)(t.hat)
  
  T <- construct.bezier.model.matrix(t.hat, r, intercept = intercept)
  p <- solve(t(T) %*% W %*% T, t(T) %*% W %*% X)
  
  # X.hat <- bezier.curve(t.hat, p, intercept = intercept)
  X.hat <- T %*% p
  
  return(list(p = p, 
              X = X.hat,
              t = t.hat))
}

compute.distances.bezier <- function(X, p, 
                                     min.t = 0, max.t = 1, 
                                     parallel = FALSE, intercept = TRUE) {
  t.hat <- estimate.t.bezier(X, p, parallel = parallel, intercept = intercept,
                             min.t = 0, max.t = 1)
  X.hat <- bezier.curve(t.hat, p, intercept = intercept)
  apply(X - X.hat, 1, function(x) sum(x ^ 2))
}

manifold.clustering <- function(X, K = 2, 
                                initialization = 'random',
                                intercept = TRUE, 
                                maxit = 50,
                                k.isomap = as.integer(sqrt(nrow(X)) / K),
                                eps = 1e-3,
                                parallel = TRUE,
                                verbose = FALSE,
                                animate = FALSE,
                                animation.dir = '.',
                                animation.title = 'plot') {
  n <- nrow(X)
  
  if (length(initialization) == n) {
    z.hat <- initialization
  } else if (initialization == 'spectral') {
    z.hat <- spectral.clustering(A, K)
  } else if (initialization == 'random') {
    z.hat <- sample(seq(K), n, replace = TRUE)
  } else {
    stop('choose a valid initialization')
  }
  
  niter <- 0
  
  loss <- c()
  new.loss <- 1 / eps
  
  if (animate) {
    iter.df <- data.frame(x = numeric(0),
                          y = numeric(0),
                          z = integer(0),
                          iter = numeric(0),
                          type = character(0))
    curve1.df <- curve2.df <- data.frame(X = numeric(0),
                                         Y = numeric(0),
                                         iter = numeric(0),
                                         type = character(0))
  }
  
  while (TRUE) {
    z.hat.prev <- z.hat
    
    X1 <- X[z.hat == 1, ]
    X2 <- X[z.hat == 2, ]
    
    if (verbose) print('fitting curve 1')
    if (exists('curve1')) {
      curve1 <- estimate.bezier.curve.2(X1, 
                                        init.params = curve1$p,
                                        k.isomap = k.isomap, 
                                        intercept = intercept,
                                        parallel = parallel)
    } else {
      curve1 <- estimate.bezier.curve.2(X1, 
                                        k.isomap = k.isomap, 
                                        intercept = intercept,
                                        parallel = parallel)
    }
    if (verbose) print('fitting curve 2')
    if (exists('curve2')) {
      curve2 <- estimate.bezier.curve.2(X2, 
                                        init.params = curve2$p,
                                        k.isomap = k.isomap, 
                                        intercept = intercept,
                                        parallel = parallel)
    } else {
      curve2 <- estimate.bezier.curve.2(X2, 
                                        k.isomap = k.isomap, 
                                        intercept = intercept,
                                        parallel = parallel)
    }
    
    if (verbose) print('reassigning clusters')
    d1 <- compute.distances.bezier(X, curve1$p, 
                                   # min.t = min(curve1$t), 
                                   # max.t = max(curve1$t),
                                   min.t = 0, max.t = 1,
                                   parallel = parallel, 
                                   intercept = intercept)
    d2 <- compute.distances.bezier(X, curve2$p, 
                                   # min.t = min(curve2$t), 
                                   # max.t = max(curve2$t),
                                   min.t = 0, max.t = 1,
                                   parallel = parallel, 
                                   intercept = intercept)
    distances <- cbind(d1, d2)
    
    # ggplot() +
    #   # viridis::scale_colour_viridis() +
    #   geom_point(aes(x = X[, 1], y = X[, 2],
    #                  colour = factor(z.hat),
    #                  shape = factor(z.hat)),
    #              size = 5) +
    #   geom_point(aes(x = curve1$X[, 1], y = curve1$X[, 2]), colour = 'red') +
    #   geom_point(aes(x = curve2$X[, 1], y = curve2$X[, 2]), colour = 'blue') +
    #   coord_fixed()
    
    if (animate) {
      update.curve1.df <- dplyr::tibble(X = curve1$X[, 1], 
                                        Y = curve1$X[, 2],
                                        iter = niter,
                                        type = 'curvefit')
      update.curve2.df <- dplyr::tibble(X = curve2$X[, 1], 
                                        Y = curve2$X[, 2],
                                        iter = niter,
                                        type = 'curvefit')
      curve1.df %<>% dplyr::bind_rows(update.curve1.df)
      curve2.df %<>% dplyr::bind_rows(update.curve2.df)
      
      update.iter.df <- dplyr::tibble(x = X[, 1], 
                                      y = X[, 2],
                                      z = z.hat,
                                      iter = niter,
                                      type = 'curvefit')
      iter.df %<>% dplyr::bind_rows(update.iter.df)
    }
    
    prev.loss <- new.loss
    curves <- list(curve1, curve2)
    new.loss <- sapply(seq(K), function(k) {
      norm(X[z.hat == k, ] - curves[[k]]$X, 'F') ^ 2
    }) %>% 
      sum()
    d.loss <- abs(prev.loss - new.loss) / new.loss
    loss <- c(loss, new.loss)
    
    z.hat <- apply(distances, 1, which.min)
    # plot(X, col = z.hat)
    # points(curve1$X, pch = '*')
    # points(curve2$X, pch = '*', col = 2)
    # ggplot() +
    #   # viridis::scale_colour_viridis() +
    #   geom_point(aes(x = X[, 1], y = X[, 2],
    #                  colour = factor(z.hat),
    #                  shape = factor(z.hat)),
    #              size = 5) +
    #   geom_point(aes(x = curve1$X[, 1], y = curve1$X[, 2]), colour = 'red') +
    #   geom_point(aes(x = curve2$X[, 1], y = curve2$X[, 2]), colour = 'blue') +
    #   geom_text(aes(x = X[, 1], y = X[, 2],
    #                 # colour = factor(z.hat),
    #                 label = seq(n)),
    #             size = 2) + 
    #   coord_fixed() + 
    #   labs(x = NULL, y = NULL, colour = NULL, shape = NULL)
    # print(loss)
    # print(curve1$mse %>% diff())
    # print(curve2$mse %>% diff())
    # print(table(z.hat, z.hat.prev))
    # print(d.loss)
    # 
    # sapply(seq(K), function(k) {
    #   X.k <- X[z.hat == k, ]
    #   p.k <- curves[[k]]$p
    #   t.k <- estimate.t.bezier(X.k, p.k,
    #                            min.t = min(curves[[k]]$t),
    #                            max.t = max(curves[[k]]$t),
    #                            parallel = parallel, intercept = intercept)
    #   T.k <- construct.bezier.model.matrix(t.k, r = 2, intercept = intercept)
    #   norm(X.k - T.k %*% p.k, 'F') ^ 2
    # }) %>%
    #   sum()
    
    if (animate) {
      update.curve1.df <- dplyr::tibble(X = curve1$X[, 1], 
                                        Y = curve1$X[, 2],
                                        iter = niter + .5,
                                        type = 'reassign')
      update.curve2.df <- dplyr::tibble(X = curve2$X[, 1], 
                                        Y = curve2$X[, 2],
                                        iter = niter + .5,
                                        type = 'reassign')
      curve1.df %<>% dplyr::bind_rows(update.curve1.df)
      curve2.df %<>% dplyr::bind_rows(update.curve2.df)
      
      update.iter.df <- dplyr::tibble(x = X[, 1], 
                                      y = X[, 2],
                                      z = z.hat,
                                      iter = niter + .5,
                                      type = 'reassign')
      iter.df %<>% dplyr::bind_rows(update.iter.df)
    }
    
    if (d.loss < eps) break
    
    niter <- niter + 1
    if (niter >= maxit) {
      warning('failed to converge to a clustering')
      break
    }
  }
  
  p.list <- lapply(curves, function(curve) curve$p)
  if (!intercept) {
    p.list <- lapply(p.list, function(p) rbind(0, p))
  }
  
  if (animate) {
    anim <- ggplot(iter.df) +
      geom_point(aes(x = x, y = y,
                     colour = factor(z)),
                 size = 1) +
      geom_point(data = curve1.df,
                 aes(x = X, y = Y), colour = 'red') +
      geom_point(data = curve2.df,
                 aes(x = X, y = Y), colour = 'blue') +
      coord_fixed() +
      labs(x = NULL, y = NULL, colour = NULL, shape = NULL) + 
      transition_states(iter,
                        wrap = FALSE,
                        transition_length = 0,
                        state_length = 2)
    anim_save(paste0(animation.title, '.gif'),
              anim,
              path = animation.dir)
  }
  
  return(list(z = z.hat, 
              X = list(curve1$X, curve2$X),
              t = list(curve1$t, curve2$t),
              p = p.list,
              niter = niter,
              loss = loss))
}

bezier.mse <- function(X, t, p, intercept = TRUE) {
  if (intercept) {
    r <- nrow(p) - 1L
  } else {
    r <- nrow(p)
  }
  T <- construct.bezier.model.matrix(t, r, intercept = intercept)
  # norm(X - T %*% p, 'F') ^ 2 / prod(dim(X))
  norm(X - T %*% p, 'F') ^ 2
}

loss.one.t <- function(x, t, p, intercept = TRUE) {
  if (intercept) {
    r <- nrow(p) - 1
    T <- sapply(seq(0, r), function(s) {
      choose(r, s) * (1 - t) ^ (r - s) * t ^ s
    })
  } else {
    r <- nrow(p)
    T <- sapply(seq(r), function(s) {
      choose(r, s) * (1 - t) ^ (r - s) * t ^ s
    })
  }
  g <- as.numeric(T %*% p)
  as.numeric(crossprod(x - g))
}

clustering.loss <- function(X.list, p.list, t.list, intercept = TRUE) {
  K <- length(X.list)
  sapply(seq(K), function(k) {
    X.k <- X.list[[k]]
    p.k <- p.list[[k]]
    t.k <- t.list[[k]]
    bezier.mse(X.k, t.k, p.k, intercept = intercept)
  }) %>% 
    sum() %>% 
    return()
}
