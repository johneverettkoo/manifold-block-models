import::from(magrittr, `%>%`)
import::from(foreach, foreach, `%do%`, `%dopar%`)
library(ggplot2)

setwd('~/dev/manifold-block-models')
source('functions.R')

n.vec <- 2 ^ c(4, 5, 6, 7, 8, 9, 10, 11)
# n.vec <- sort(n.vec, decreasing = TRUE)
nsamp.vec <- c(0, 4, 8)
iter <- 50

doMC::registerDoMC(parallel::detectCores() - 0)

setwd('simulations/nonintersect')

reorder.mat <- gtools::permutations(2, 2)

clustering.df <- foreach(nsamp = nsamp.vec, .combine = dplyr::bind_rows) %do% {
  out.df <- foreach(n = n.vec, .combine = dplyr::bind_rows) %do% {
    print(paste('nsamp =', nsamp, ', n = ', n))
    if (paste0('nonintersect_n_', n, '_nsamp_', nsamp, '.csv') %in% dir()) {
      out <- readr::read_csv(paste0('nonintersect_n_', n, '_nsamp_', nsamp, '.csv'))
    } else {
      out <- foreach(i = seq(iter), .combine = dplyr::bind_rows, .errorhandling = 'remove') %dopar% {
        clustering <- simulate.and.compute.error.nonintersect(
          n, 
          initialization = 'ground truth',
          ground.truth.sample = nsamp,
          maxit = 32
        )
        loss <- clustering$loss
        error.curves <- clustering$error.curves
        error.singlelink <- clustering$error.singlelink
        dplyr::tibble(n = n,
                      nsamp = nsamp,
                      loss = loss,
                      error.curves = error.curves,
                      error.singlelink = error.singlelink)
      }
      readr::write_csv(out, paste0('nonintersect_n_', n, '_nsamp_', nsamp, '.csv'))
    }
    return(out)
  }
  gc()
  return(out.df)
}

gc()

readr::write_csv(
  clustering.df, '~/dev/manifold-block-models/simulations/nonintersect/nonintersect.csv')

clust.summary.df <- clustering.df %>% 
  dplyr::group_by(n, nsamp) %>% 
  dplyr::summarise(
    med.curves = median(error.curves),
    first.q.curves = quantile(error.curves, .25),
    third.q.curves = quantile(error.curves, .75),
    med.singlelink = median(error.singlelink),
    first.q.singlelink = quantile(error.singlelink, .25),
    third.q.singlelink = quantile(error.singlelink, .75)
  ) %>% 
  dplyr::ungroup()

ggplot(clust.summary.df) + 
  theme_bw() + 
  theme(text = element_text(size = 10)) + 
  scale_y_log10() +
  scale_x_log10(breaks = n.vec) +
  labs(y = 'community detection error rate',
       colour = NULL, shape = NULL) + 
  geom_line(aes(x = n, y = med.curves, colour = factor(nsamp))) + 
  geom_point(aes(x = n, y = med.curves, colour = factor(nsamp))) + 
  geom_errorbar(aes(x = n, 
                    ymin = first.q.curves, ymax = third.q.curves, 
                    colour = factor(nsamp))) + 
  scale_colour_brewer(palette = 'Set1')

ggplot(clust.summary.df) + 
  theme_bw() + 
  theme(text = element_text(size = 10)) + 
  # scale_y_log10() +
  scale_x_log10(breaks = n.vec) +
  labs(y = 'community detection error rate',
       colour = NULL, shape = NULL) + 
  geom_line(aes(x = n, y = med.singlelink, colour = factor(nsamp))) + 
  geom_point(aes(x = n, y = med.singlelink, colour = factor(nsamp))) + 
  geom_errorbar(aes(x = n, 
                    ymin = first.q.singlelink, ymax = third.q.singlelink, 
                    colour = factor(nsamp))) + 
  scale_colour_brewer(palette = 'Set1')
