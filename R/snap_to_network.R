#' Snap a GPS coordinate to the nearest stream network segment
#'
#' @param lat Numeric. Latitude in decimal degrees (WGS84).
#' @param lon Numeric. Longitude in decimal degrees (WGS84).
#' @param network An \code{sf} object of stream network line segments.
#'   If \code{NULL}, NHD flowlines are fetched automatically.
#' @param buffer_km Numeric. Buffer radius (km) for NHD download. Default 50.
#' @param min_stream_order Integer. Minimum NHD stream order to snap to when
#'   using automatic NHD data. Prevents snapping to small headwater ditches.
#'   Default is 3. Ignored when \code{network} is user-supplied.
#'
#' @return A named list with \code{point} (sf), \code{network} (sf), and
#'   \code{segment_id} (row index of nearest qualifying segment).
#' @export
snap_to_network <- function(lat, lon, network = NULL, buffer_km = 50,
                            min_stream_order = 3) {
  point <- sf::st_sfc(sf::st_point(c(lon, lat)), crs = 4326)
  point_sf <- sf::st_sf(geometry = point)

  if (is.null(network)) {
    message("No network provided. Fetching NHD flowlines via nhdplusTools...")
    aoi <- sf::st_buffer(
      sf::st_transform(point_sf, crs = 5070),
      dist = buffer_km * 1000
    )
    aoi <- sf::st_transform(aoi, crs = 4326)
    network <- nhdplusTools::get_nhdplus(AOI = aoi, realization = "flowline")
    network <- sf::st_transform(network, crs = 4326)

    candidates <- network[network$streamorde >= min_stream_order, ]
    if (nrow(candidates) == 0) {
      stop("No streams of order >= ", min_stream_order,
           " found within ", buffer_km, " km. Try reducing min_stream_order or increasing buffer_km.")
    }
    nearest_row <- sf::st_nearest_feature(point_sf, candidates)
    segment_id <- which(network$comid == candidates$comid[nearest_row])[1]
  } else {
    network <- sf::st_transform(network, crs = 4326)
    segment_id <- sf::st_nearest_feature(point_sf, network)
  }

  list(point = point_sf, network = network, segment_id = segment_id)
}
