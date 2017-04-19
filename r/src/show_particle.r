#' show_particle visualizes particle trajectories on a 3d plot
#' @author Ben Fasoli
#'
#' @param particle_rds location of PARTICLE.rds file
#' @param outfile location to save html for interactive figure
#'
#' @import plotly
#' @export

show_particle <- function(particle_rds, outfile = 'plot_particle.html') {
  require(plotly)

  particle <- readRDS(particle_rds)$data

  f <- plot_ly(
    type = 'scatter3d', mode = 'markers',
    x = particle$long,
    y = particle$lati,
    z = particle$zagl,
    color = particle$foot,
    text = paste('Index:', particle$indx, '\n',
                 'Footprint:', particle$foot),
    marker = list(
      opacity = 0.2,
      size = 5
    )
  ) %>%
  layout(
    scene = list(
      camera = list(
        eye = list(
          x = 0,
          y = -2,
          z = 2
        )
      )
    )
  )

  htmlwidgets::saveWidget(f, file = outfile)
  print(f)
}
