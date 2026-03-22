#' Extract the upstream hydrologic network from a GPS coordinate
#'
#' For automatic NHD data, uses local hydrosequence topology for accurate
#' upstream routing. For user-supplied networks, traverses the graph using igraph.
#'
#' @param lat Numeric. Latitude in decimal degrees (WGS84).
#' @param lon Numeric. Longitude in decimal degrees (WGS84).
#' @param network An \code{sf} object of stream network line segments.
#'   If \code{NULL}, NHD data is fetched and traversed automatically.
#' @param from_col Character. Column for upstream node ID in user-supplied
#'   networks. Default \code{"fromnode"}.
#' @param to_col Character. Column for downstream node ID in user-supplied
#'   networks. Default \code{"tonode"}.
#' @param buffer_km Numeric. Buffer radius (km) for NHD download. Default 50.
#' @param min_stream_order Integer. Minimum NHD stream order to snap to.
#'   Default 3. Ignored for user-supplied networks.
#'
#' @return An \code{sf} object of all stream segments upstream of the input point.
#' @export
extract_upstream <- function(lat, lon, network = NULL,
                             from_col = "fromnode", to_col = "tonode",
                             buffer_km = 50, min_stream_order = 3) {
  snapped <- snap_to_network(lat, lon, network = network,
                             buffer_km = buffer_km,
                             min_stream_order = min_stream_order)
  net <- snapped$network
  start_idx <- snapped$segment_id

  if (is.null(network)) {
    # NHD path: derive tocomid from hydroseq (consistent across VPU boundaries)
    hydroseq_to_comid <- setNames(net$comid, as.character(net$hydroseq))
    net$tocomid <- hydroseq_to_comid[as.character(net$dnhydroseq)]
    net$tocomid[is.na(net$tocomid)] <- 0

    start_comid <- net$comid[start_idx]
    upstream_comids <- nhdplusTools::get_UT(
      network = sf::st_drop_geometry(net),
      comid = start_comid
    )
    return(net[net$comid %in% upstream_comids, ])
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
