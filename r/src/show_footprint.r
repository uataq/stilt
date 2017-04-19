#' show_footprint visualizes a footprint .rds file on an interactive map
#' @author Ben Fasoli
#' 
#' @param footprint_rds location of FOOTPRINT.rds file
#' @param outfile location to save html for interactive figure
#' 
#' @import dplyr, htmlwidgets, leaflet, raster
#' @export

show_footprint <- function(footprint_rds, footprint_html = 'plot_footprint.html') {
  require(dplyr)
  require(htmlwidgets)
  require(leaflet)
  require(raster)
  
  footprint <- readRDS(footprint_rds)
  
  r <- footprint$data %>%
    .[c('long', 'lati', 'foot')] %>%
    filter(foot >= 0) %>%
    raster::rasterFromXYZ(crs = '+init=epsg:4326') %>%
    raster::trim()
  
  lvls <- pretty(footprint$data$foot)
  cpal <- colorNumeric(RColorBrewer::brewer.pal(11, 'Spectral') %>% rev,
                       lvls, na.color = '#00000000')
  
  l <- leaflet() %>%
    addProviderTiles('CartoDB.Positron') %>%
    fitBounds(r@extent@xmin, r@extent@ymin,
              r@extent@xmax, r@extent@ymax) %>%
    addRasterImage(r, colors = cpal, opacity = 0.3) %>%
    addCircleMarkers(footprint$receptors$r_long, footprint$receptors$r_lati,
                     color = 'black', radius = 8, weight = 2, fillOpacity = 0.1)
  
  saveWidget(l, file = footprint_html)
  print(l)
}