#' Extract the upstream hydrologic network from a GPS coordinate
#'
#' For automatic NHD data, the point is snapped to the nearest qualifying
#' stream, the full drainage basin upstream of it is surveyed, and only
#' segments of at least the chosen stream order are downloaded. For
#' user-supplied networks, traverses the graph using igraph.
#'
#' Whole basins can contain many thousands of segments, so for NHD data the
#' resolution (minimum stream order) must be chosen. If \code{resolution} is
#' \code{NULL} in an interactive session, the basin is surveyed and a menu of
#' resolutions is shown with estimated segment counts and download times. In
#' a non-interactive session the survey is printed and an error asks for
#' \code{resolution} to be set.
#'
#' @param lat Numeric. Latitude in decimal degrees (WGS84).
#' @param lon Numeric. Longitude in decimal degrees (WGS84).
#' @param network An \code{sf} object of stream network line segments.
#'   If \code{NULL}, NHD data is fetched and traversed automatically.
#' @param resolution Integer or \code{NULL}. Minimum NHD stream order to
#'   include in the extracted network: 1 keeps every mapped stream, higher
#'   values keep only larger streams. If \code{NULL} (default), the user is
#'   asked. Ignored for user-supplied networks.
#' @param from_col Character. Column for upstream node ID in user-supplied
#'   networks. Default \code{"fromnode"}.
#' @param to_col Character. Column for downstream node ID in user-supplied
#'   networks. Default \code{"tonode"}.
#' @param buffer_km Numeric. Search radius (km) for the nearest qualifying
#'   stream when snapping the point. Default 10.
#' @param min_stream_order Integer. Minimum NHD stream order to snap to.
#'   Default 3. Ignored for user-supplied networks.
#' @param max_area_km2 Numeric. Largest drainage area (km^2) to extract.
#'   Default 100000. See \code{survey_upstream}.
#'
#' @return An \code{sf} object of all stream segments upstream of the input
#'   point (of at least order \code{resolution} for NHD data).
#' @export
extract_upstream <- function(lat, lon, network = NULL, resolution = NULL,
                             from_col = "fromnode", to_col = "tonode",
                             buffer_km = 10, min_stream_order = 3,
                             max_area_km2 = 1e5) {
  snapped <- snap_to_network(lat, lon, network = network,
                             buffer_km = buffer_km,
                             min_stream_order = min_stream_order)
  net <- snapped$network
  start_idx <- snapped$segment_id

  if (is.null(network)) {
    start_comid <- net$comid[start_idx]
    message("Snapped to COMID ", start_comid, " (",
            round(as.numeric(sf::st_distance(snapped$point,
                                             net[start_idx, ]))),
            " m away).")

    if (is.null(resolution)) {
      message("Surveying upstream basin...")
      survey <- survey_upstream(start_comid, max_area_km2 = max_area_km2)
      if (!interactive()) {
        print(survey)
        stop("Choose a resolution from the table above and pass it as ",
             "`resolution` (minimum stream order).", call. = FALSE)
      }
      print_survey_header(survey)
      choice <- utils::menu(
        paste0(survey_labels(survey), ": ",
               format_count(survey$options$est_segments), " segments, ",
               format_seconds(survey$options$est_seconds), " download"),
        title = "\nChoose a resolution (0 to cancel):"
      )
      if (choice == 0) stop("Cancelled.", call. = FALSE)
      resolution <- survey$options$min_order[choice]
      basin <- survey$basin
    } else {
      basin <- nhdplusTools::get_nldi_basin(
        list(featureSource = "comid", featureID = as.character(start_comid))
      )
    }

    message("Downloading streams of order >= ", resolution, "...")
    up <- nhdplusTools::get_nhdplus(
      AOI = basin, realization = "flowline",
      streamorder = if (resolution > 1) resolution else NULL
    )
    up <- sf::st_transform(up, crs = 4326)

    # The basin polygon can clip in neighbouring streams along its edge, so
    # keep only segments connected upstream of the start segment
    up <- upstream_by_hydroseq(up, start_comid)
    message("Extracted ", nrow(up), " segments (",
            round(sum(up$lengthkm)), " km of stream).")
    return(up)
  }

  # User-provided network: igraph traversal
  if (!from_col %in% names(net) || !to_col %in% names(net)) {
    stop(
      "Columns '", from_col, "' and '", to_col, "' not found in network.\n",
      "Available columns: ", paste(names(net), collapse = ", "), "\n",
      "Use from_col and to_col to specify the correct column names."
    )
  }

  edges <- data.frame(
    from = as.character(net[[from_col]]),
    to   = as.character(net[[to_col]]),
    stringsAsFactors = FALSE
  )

  g <- igraph::graph_from_data_frame(edges, directed = TRUE)
  outlet_node <- as.character(net[[to_col]][start_idx])

  if (!outlet_node %in% igraph::V(g)$name) {
    stop("Outlet node '", outlet_node, "' not found in network graph.")
  }

  upstream_nodes <- names(igraph::subcomponent(g, outlet_node, mode = "in"))
  mask <- net[[from_col]] %in% upstream_nodes | net[[to_col]] %in% upstream_nodes
  net[mask, ]
}

# Keep segments that drain to start_comid, following dnhydroseq upstream.
# Minor divergences (braided side channels) are dropped so the result is
# strictly dendritic; otherwise each would appear as a spurious tip.
upstream_by_hydroseq <- function(net, start_comid) {
  if ("divergence" %in% names(net))
    net <- net[net$divergence != 2 | net$comid == start_comid, ]

  start <- which(net$comid == start_comid)
  if (length(start) == 0)
    stop("Start segment COMID ", start_comid, " missing from download.")

  hydroseq <- as.character(net$hydroseq)
  children <- split(seq_len(nrow(net)), as.character(net$dnhydroseq))
  keep     <- logical(nrow(net))
  frontier <- start
  while (length(frontier) > 0) {
    keep[frontier] <- TRUE
    nxt <- unlist(children[hydroseq[frontier]], use.names = FALSE)
    frontier <- nxt[!keep[nxt]]
  }
  net[keep, ]
}
