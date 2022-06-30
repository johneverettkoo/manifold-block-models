source('~/dev/manifold-block-models/functions.R')
source('~/dev/pabm-grdpg/functions.R')
import::from(magrittr, `%>%`)
import::from(foreach, foreach, `%do%`, `%dopar%`)
library(ggplot2)

n.vec <- 2 ^ c(7, 8, 9, 10)
n.vec <- sort(n.vec, decreasing = TRUE)
nsamp.vec <- c(6, 8, 10)
iter <- 50

doMC::registerDoMC(parallel::detectCores())

clustering.df <- foreach(nsamp = nsamp.vec, .combine = dplyr::bind_rows) %do% {
  reorder.mat <- gtools::permutations(2, 2)
  out.df <- foreach(n = n.vec, .combine = dplyr::bind_rows) %do% {
    print(paste('nsamp =', nsamp, ', n = ', n))
    foreach(i = seq(iter), .combine = dplyr::bind_rows, .errorhandling = 'remove') %dopar% {
      error.rate <- simulate.and.compute.error(n, 
                                               initialization = 'ground truth',
                                               ground.truth.sample = nsamp)
      error.count = error.rate * n
      out <- dplyr::tibble(n = n,
                           nsamp = nsamp,
                           error.rate = error.rate,
                           error.count = error.count)
      gc()
      return(out)
    } %>% 
      return()
  }
  gc()
  return(out.df)
}

gc()

readr::write_csv(clustering.df, '~/dev/manifold-block-models/simulations/balanced.csv')

clust.summary.df <- clustering.df %>% 
  dplyr::group_by(n, nsamp) %>% 
  dplyr::summarise(
    med.count = median(error.count),
    first.q.count = quantile(error.count, .25),
    third.q.count = quantile(error.count, .75),
    med.rate = median(error.rate),
    first.q.rate = quantile(error.rate, .25),
    third.q.rate = quantile(error.rate, .75)
  ) %>% 
  dplyr::ungroup()

ggplot(clust.summary.df) + 
  theme_bw() + 
  theme(text = element_text(size = 10)) + 
  scale_y_log10() + 
  scale_x_log10(breaks = n.vec) + 
  labs(y = 'community detection error rate',
       colour = NULL, shape = NULL) + 
  geom_line(aes(x = n, y = med.rate, colour = factor(nsamp))) + 
  geom_point(aes(x = n, y = med.rate, colour = factor(nsamp))) + 
  geom_errorbar(aes(x = n, ymin = first.q.rate, ymax = third.q.rate, colour = factor(nsamp))) + 
  scale_colour_brewer(palette = 'Set1')
