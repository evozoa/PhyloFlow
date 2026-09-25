#' Survey the NHD network upstream of a stream segment
#'
#' Summarises the drainage basin upstream of an NHDPlus segment and estimates,
#' for each possible resolution (minimum stream order), how many segments the
#' extracted network would contain and how long downloading it would take.
#' Used by \code{extract_upstream} to let the user choose a resolution before
#' committing to a large download.
#'
#' The NHDPlus web service can only be queried by bounding box, so the
#' survey asks it how many segments of each order fall in the basin's
#' bounding box (without downloading them). Download time is estimated from
#' that count, since the whole box must be fetched; the number of segments
#' in the basin itself is estimated by scaling it by the fraction of the box
#' the basin covers.
#'
#' @param comid Integer. NHDPlus COMID of the outlet segment.
#' @param max_area_km2 Numeric. Largest drainage area (km^2) to survey. Very
#'   large basins (e.g. major rivers) can contain hundreds of thousands of
#'   segments; choose a point further upstream or supply your own network.
#'   Default 100000.
#'
#' @return A list of class \code{"phyloflow_survey"} with elements
#'   \code{comid}, \code{outlet_name}, \code{outlet_order}, \code{area_km2},
#'   \code{basin} (\code{sf} polygon), and \code{options}, a data frame with
#'   one row per resolution giving \code{min_order}, \code{est_segments} and
#'   \code{est_seconds}.
#' @export
survey_upstream <- function(comid, max_area_km2 = 1e5) {
  basin <- nhdplusTools::get_nldi_basin(
    list(featureSource = "comid", featureID = as.character(comid))
  )
  area_km2 <- as.numeric(sf::st_area(sf::st_transform(basin, 5070))) / 1e6
  if (area_km2 > max_area_km2) {
    stop("Drainage area upstream of COMID ", comid, " is ",
         format(round(area_km2), big.mark = ","), " km^2, larger than ",
         "max_area_km2 (", format(max_area_km2, big.mark = ","), "). ",
         "Choose a point further upstream, supply your own network, or ",
         "raise max_area_km2.")
  }

  outlet <- nhdplusTools::get_nhdplus(
    comid = as.integer(comid), realization = "flowline",
    properties = c("comid", "streamorde", "gnis_name"), skip_geometry = TRUE
  )
  outlet_order <- outlet$streamorde[1]

  # The outlet's own order yields a single unbranched main stem, so the
  # coarsest useful resolution is one order below it
  min_orders <- seq_len(max(outlet_order - 1, 1))

  bbox <- sf::st_bbox(sf::st_transform(basin, 4326))
  bbox_km2 <- as.numeric(sf::st_area(
    sf::st_transform(sf::st_as_sfc(bbox), 5070))) / 1e6
  in_bbox <- count_flowlines(bbox, min_orders)

  # Download time fitted to NHDPlus OGC API timings (seconds vs. segments
  # in the bounding box, which is what get_nhdplus() pages through)
  structure(
    list(
      comid        = as.integer(comid),
      outlet_name  = outlet$gnis_name[1],
      outlet_order = outlet_order,
      area_km2     = area_km2,
      basin        = basin,
      options      = data.frame(
        min_order    = min_orders,
        est_segments = round(in_bbox * area_km2 / bbox_km2),
        est_seconds  = 3 + 0.0085 * in_bbox
      )
    ),
    class = "phyloflow_survey"
  )
}

# Number of NHDPlus flowlines of at least each order within a bounding box,
# using the OGC API's numberMatched so no features are downloaded
count_flowlines <- function(bbox, min_orders) {
  base <- paste0("https://api.water.usgs.gov/fabric/pygeoapi/",
                 "collections/nhdflowline_network/items")
  urls <- paste0(base, "?f=json&limit=1&skipGeometry=true",
                 "&bbox=", paste(round(bbox, 5), collapse = ","),
                 "&filter=", utils::URLencode(paste0("streamorde >= ",
                                                     min_orders),
                                              reserved = TRUE))

  counts <- rep(NA_real_, length(urls))
  pool <- curl::new_pool()
  for (i in seq_along(urls)) {
    local({
      j <- i
      curl::curl_fetch_multi(urls[j], pool = pool, done = function(res) {
        if (res$status_code == 200) {
          body <- jsonlite::fromJSON(rawToChar(res$content))
          if (!is.null(body$numberMatched)) counts[j] <<- body$numberMatched
        }
      })
    })
  }
  curl::multi_run(pool = pool)

  if (anyNA(counts))
    stop("Could not count NHDPlus flowlines; the USGS web service may be ",
         "unavailable. Pass `resolution` to skip the survey.")
  counts
}

#' @export
print.phyloflow_survey <- function(x, ...) {
  print_survey_header(x)
  cat("\n")
  print(data.frame(
    resolution = survey_labels(x),
    segments   = format_count(x$options$est_segments),
    download   = format_seconds(x$options$est_seconds)
  ), row.names = FALSE, right = FALSE)
  invisible(x)
}

print_survey_header <- function(x) {
  name <- if (is.na(x$outlet_name) || x$outlet_name == "") "unnamed stream"
          else x$outlet_name
  cat("Upstream of COMID ", x$comid, " (", name, ", stream order ",
      x$outlet_order, ")\n", sep = "")
  cat("  Drainage area ", format(round(x$area_km2), big.mark = ","),
      " km^2, ~", format(signif(x$options$est_segments[1], 2),
                          big.mark = ","),
      " stream segments\n", sep = "")
}

survey_labels <- function(x) {
  k <- x$options$min_order
  label <- paste0("Order >= ", k)
  label[k == 1] <- paste0(label[k == 1], " (every mapped stream)")
  coarsest <- k == max(k) & max(k) > 1
  label[coarsest] <- paste0(label[coarsest], " (major tributaries only)")
  label
}

format_count <- function(n) {
  paste0("~", vapply(signif(n, 2), format, character(1),
                     big.mark = ",", scientific = FALSE))
}

format_seconds <- function(s) {
  ifelse(s < 90,
         paste0("~", pmax(5, round(s / 5) * 5), " s"),
         paste0("~", round(s / 60, 1), " min"))
}
