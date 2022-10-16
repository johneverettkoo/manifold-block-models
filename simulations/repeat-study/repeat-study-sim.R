import::from(magrittr, `%>%`)
import::from(foreach, foreach, `%do%`, `%dopar%`)
library(ggplot2)

setwd('~/dev/manifold-block-models')
source('functions.R')

n <- 512
nsamp <- 0
iter.per.core <- 32
maxit <- 32

ncores <- parallel::detectCores()
doMC::registerDoMC(ncores)

setwd('simulations/repeat-study')

reorder.mat <- gtools::permutations(2, 2)

clustering.df <- foreach(i = seq(ncores), .combine = dplyr::bind_rows) %dopar% {
  if (paste0('sim_', i, '.csv') %in% dir()) {
    out <- readr::read_csv(paste0('sim_', i, '.csv'))
  } else {
    out <- foreach(j = seq(iter.per.core), .combine = dplyr::bind_rows, .errorhandling = 'remove') %do% {
      clustering <- simulate.and.compute.error(
        n, 
        initialization = 'ground truth',
        ground.truth.sample = nsamp,
        maxit = maxit
      )
      error.rate <- clustering$error
      loss <- clustering$loss
      p <- clustering$p
      error.count = error.rate * n
      dplyr::tibble(n = n,
                    nsamp = nsamp,
                    error.rate = error.rate,
                    error.count = error.count,
                    loss = loss,
                    p121 = p[[1]][2, 1],
                    p122 = p[[1]][2, 2],
                    p131 = p[[1]][3, 1],
                    p132 = p[[1]][3, 2],
                    p221 = p[[2]][2, 1],
                    p222 = p[[2]][2, 2],
                    p231 = p[[2]][3, 1],
                    p232 = p[[2]][3, 2])
    }
    readr::write_csv(out, paste0('sim_', i, '.csv'))
  }
  return(out)
}

readr::write_csv(clustering.df, 'repeat-study.csv')

plot.repeat.study <- function(repeat.df) {
  repeat.df <- na.omit(repeat.df)
  min.t <- 0
  max.t <- 1
  length.out <- 1e2
  t. <- seq(min.t, max.t, length.out = length.out)
  niter <- nrow(repeat.df)
  p.df <- repeat.df[c('p121', 'p122', 'p131', 'p132', 
                      'p221', 'p222', 'p231', 'p232')]
  p.clust <- kmeans(p.df, 4)$cluster
  loss.clust <- kmeans(repeat.df$loss, 2)$cluster
  curves.df <- foreach::foreach(i = seq(niter), .combine = dplyr::bind_rows) %dopar% {
    p1 <- matrix(c(repeat.df$p121[i], repeat.df$p122[i], 
                   repeat.df$p131[i], repeat.df$p132[i]),
                 nrow = 2, byrow = TRUE)
    p1 <- rbind(0, p1)
    p2 <- matrix(c(repeat.df$p221[i], repeat.df$p222[i],
                   repeat.df$p231[i], repeat.df$p232[i]),
                 nrow = 2, byrow = TRUE)
    p2 <- rbind(0, p2)
    p <- list(p1, p2)
    plyr::ldply(seq_along(p), function(k) {
      curve <- bezier.curve(t., p[[k]])
      if (mean(curve[, 1]) < 0) {
        curve[, 1] <- -curve[, 1]
      }
      dplyr::tibble(t = t.,
                    x = curve[, 1],
                    y = curve[, 2],
                    z = k,
                    iter = i,
                    p.clust = p.clust[i],
                    loss.clust = loss.clust[i])
    })
  }
  ggplot(curves.df) + 
    geom_path(aes(x = x, y = y, group = interaction(z, iter),
                  colour = factor(loss.clust)),
              alpha = 1 / log(niter)) + 
    theme_bw() + 
    facet_wrap(~loss.clust) + 
    coord_fixed() + 
    labs(colour = NULL) + 
    scale_colour_brewer(palette = 'Set1') + 
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
}

plot.repeat.study(clustering.df)
