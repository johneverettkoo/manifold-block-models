import::from(magrittr, `%>%`, `%<>%`)
import::from(foreach, foreach, `%do%`, `%dopar%`)
library(ggplot2)

setwd('~/dev/manifold-block-models')
source('functions.R')

n.vec <- 2 ^ c(4, 5, 6, 7, 8, 9, 10, 11)
# n.vec <- sort(n.vec, decreasing = TRUE)
nsamp.vec <- c(0, 4, 8)
iter <- 50

cores.to.reserve <- 0

doMC::registerDoMC(parallel::detectCores() - cores.to.reserve)

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

error.singlelink.df <- clustering.df %>% 
  dplyr::select(n, nsamp, error.singlelink) %>% 
  dplyr::mutate(nsamp = 'single linkage') %>% 
  dplyr::rename(error = error.singlelink)
clustering.df %<>% 
  dplyr::select(n, nsamp, error = error.curves) %>% 
  dplyr::mutate(nsamp = as.character(nsamp)) %>% 
  dplyr::bind_rows(error.singlelink.df)

readr::write_csv(
  clustering.df, '~/dev/manifold-block-models/simulations/nonintersect/nonintersect.csv')

clust.summary.df <- clustering.df %>% 
  dplyr::group_by(n, nsamp) %>% 
  dplyr::summarise(
    med.error = median(error),
    first.q.error = quantile(error, .25),
    third.q.error = quantile(error, .75)
  ) %>% 
  dplyr::ungroup()

clust.summary.df %>% 
  dplyr::filter(n <= 512) %>% 
  ggplot() + 
  theme_bw() + 
  theme(text = element_text(size = 10)) + 
  scale_y_log10() +
  scale_x_log10(breaks = n.vec) +
  labs(y = 'community detection error rate',
       colour = NULL, shape = NULL) + 
  geom_line(aes(x = n, y = med.error, colour = factor(nsamp))) + 
  geom_point(aes(x = n, y = med.error, colour = factor(nsamp))) + 
  geom_errorbar(aes(x = n, 
                    ymin = first.q.error, ymax = third.q.error, 
                    colour = factor(nsamp))) + 
  scale_colour_brewer(palette = 'Set1')
