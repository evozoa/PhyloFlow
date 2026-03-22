#' Compute pairwise stream network distances between sample sites
#'
#' Computes a matrix of pairwise distances between sample collection sites
#' measured along the stream network. Suitable for use in Mantel tests and
#' other isolation-by-distance analyses with nuclear markers.
#'
#' Optionally, segment lengths can be weighted by any NHD attribute (e.g.
#' \code{"slope"} or \code{"qa_ma"} for mean annual discharge) to produce
#' resistance-based distances.
#'
#' @param sites A data frame with one row per sample site.
#' @param site_id_col Character. Column name for sample IDs. Default
#'   \code{"sample_id"}.
#' @param lat_col Character. Column name for latitude (decimal degrees, WGS84).
#'   Default \code{"lat"}.
#' @param lon_col Character. Column name for longitude (decimal degrees,
#'   WGS84). Default \code{"lon"}.
#' @param network An \code{sf} object of stream network line segments. If
#'   \code{NULL}, NHD data is fetched automatically.
#' @param buffer_km Numeric. Buffer (km) around the bounding box of all sites
#'   for NHD download. Default 50.
#' @param min_stream_order Integer. Minimum NHD stream order to snap samples
#'   to. Default 3.
#' @param length_col Character. Column name for segment lengths. Default
#'   \code{"lengthkm"}. Set to \code{NULL} to compute from geometry in metres.
#' @param weight_col Character or \code{NULL}. Optional NHD column to use as a
#'   resistance weight. Segment lengths are multiplied by this value before
#'   summing path distances. For example, \code{"slope"} weights distances by
#'   channel gradient. Default \code{NULL} (unweighted stream distance).
#' @param format Character. \code{"matrix"} (default) returns a symmetric
#'   matrix; \code{"dist"} returns an object of class \code{dist}.
#'
#' @return A symmetric distance matrix (or \code{dist} object) of pairwise
#'   stream network distances, with rows and columns named by sample ID.
#' @export
sites_to_distances <- function(sites,
                               site_id_col      = "sample_id",
                               lat_col          = "lat",
                               lon_col          = "lon",
                               network          = NULL,
                               buffer_km        = 50,
                               min_stream_order = 3,
                               length_col       = "lengthkm",
                               weight_col       = NULL,
                               format           = "matrix") {

  # --- 1. Validate input ---
  req     <- c(site_id_col, lat_col, lon_col)
  missing <- req[!req %in% names(sites)]
  if (length(missing) > 0)
    stop("Missing columns in sites: ", paste(missing, collapse = ", "))

  sample_ids <- as.character(sites[[site_id_col]])
  lats       <- sites[[lat_col]]
  lons       <- sites[[lon_col]]
  n_sites    <- nrow(sites)

  if (n_sites < 2)
    stop("At least 2 sample sites are required.")

  # --- 2. Download NHD for bounding box + buffer ---
  if (is.null(network)) {
    message("Fetching NHD flowlines for study area...")
    all_pts  <- sf::st_sfc(
      lapply(seq_len(n_sites), function(i) sf::st_point(c(lons[i], lats[i]))),
      crs = 4326
    )
    bbox_sf  <- sf::st_as_sfc(sf::st_bbox(all_pts))
    bbox_sf  <- sf::st_sf(geometry = bbox_sf, crs = 4326)
    bbox_buf <- sf::st_buffer(sf::st_transform(bbox_sf, crs = 5070),
                              dist = buffer_km * 1000)
    aoi      <- sf::st_transform(bbox_buf, crs = 4326)
    network  <- nhdplusTools::get_nhdplus(AOI = aoi, realization = "flowline")
    network  <- sf::st_transform(network, crs = 4326)
  }

  seg_ids <- as.character(network$comid)

  # --- 3. Snap each sample to nearest qualifying stream ---
  net_qual <- network[network$streamorde >= min_stream_order, ]
  if (nrow(net_qual) == 0)
    stop("No streams of order >= ", min_stream_order,
         " found. Try reducing min_stream_order or increasing buffer_km.")

  site_comids <- character(n_sites)
  for (i in seq_len(n_sites)) {
    pt             <- sf::st_sf(geometry = sf::st_sfc(
      sf::st_point(c(lons[i], lats[i])), crs = 4326))
    nearest        <- sf::st_nearest_feature(pt, net_qual)
    site_comids[i] <- as.character(net_qual$comid[nearest])
  }
  names(site_comids) <- sample_ids
  message("All samples snapped to stream network.")

  # --- 4. Build hydroseq topology ---
  hydroseq_to_seg   <- stats::setNames(seg_ids, as.character(network$hydroseq))
  parent_seg        <- hydroseq_to_seg[as.character(network$dnhydroseq)]
  names(parent_seg) <- seg_ids

  # --- 5. Find lowest common ancestor (LCA) of all sample COMIDs ---
  get_path_to_outlet <- function(start_comid) {
    path    <- character(0)
    current <- start_comid
    visited <- character(0)
    while (!is.na(current) && current %in% seg_ids && !current %in% visited) {
      visited <- c(visited, current)
      path    <- c(path, current)
      current <- parent_seg[current]
    }
    path
  }

  unique_sample_comids <- unique(site_comids)
  paths  <- lapply(unique_sample_comids, get_path_to_outlet)
  common <- Reduce(intersect, paths)

  if (length(common) == 0)
    stop("No common ancestor found for all sample sites. ",
         "Try increasing buffer_km.")

  hydroseq_vals <- network$hydroseq[match(common, seg_ids)]
  lca_comid     <- common[which.max(hydroseq_vals)]
  message("Common outlet found: COMID ", lca_comid)

  # --- 6. Extract subnetwork upstream of LCA ---
  net_df         <- sf::st_drop_geometry(network)
  net_df$tocomid <- as.integer(hydroseq_to_seg[as.character(network$dnhydroseq)])
  net_df$tocomid[is.na(net_df$tocomid)] <- 0L

  upstream_ids <- nhdplusTools::get_UT(network = net_df,
                                       comid   = as.integer(lca_comid))
  upstream_net <- network[network$comid %in% upstream_ids, ]
  up_seg_ids   <- as.character(upstream_net$comid)
  up_parent    <- stats::setNames(
    hydroseq_to_seg[as.character(upstream_net$dnhydroseq)],
    up_seg_ids
  )

  # --- 7. Build minimal connecting network ---
  keep_set <- character(0)
  for (sc in unique_sample_comids) {
    current <- sc
    while (!is.na(current) && current %in% up_seg_ids) {
      keep_set <- c(keep_set, current)
      current  <- up_parent[current]
    }
  }
  keep_set <- unique(keep_set)

  pruned_net <- upstream_net[upstream_net$comid %in% keep_set, ]
  p_seg_ids  <- as.character(pruned_net$comid)

  # --- 8. Compute segment lengths (with optional resistance weighting) ---
  if (!is.null(length_col) && length_col %in% names(pruned_net)) {
    lengths <- stats::setNames(pruned_net[[length_col]], p_seg_ids)
  } else {
    net_m   <- sf::st_transform(pruned_net, crs = 3857)
    lengths <- stats::setNames(as.numeric(sf::st_length(net_m)), p_seg_ids)
  }

  if (!is.null(weight_col)) {
    if (!weight_col %in% names(pruned_net))
      stop("weight_col '", weight_col, "' not found in network. ",
           "Available columns: ", paste(names(pruned_net), collapse = ", "))
    weights <- stats::setNames(pruned_net[[weight_col]], p_seg_ids)
    if (any(is.na(weights))) {
      warning("NA values in weight_col '", weight_col, "'; replacing with 1.")
      weights[is.na(weights)] <- 1
    }
    lengths <- lengths * weights
    message("Distances weighted by '", weight_col, "'.")
  }

  # --- 9. Compute cumulative distance from root ---
  p_parent     <- stats::setNames(
    hydroseq_to_seg[as.character(pruned_net$dnhydroseq)],
    p_seg_ids
  )
  p_has_parent <- !is.na(p_parent) & p_parent %in% p_seg_ids

  root_candidates <- p_seg_ids[!p_has_parent]
  root_seg <- root_candidates[
    which.min(pruned_net$hydroseq[!p_has_parent])
  ]

  # BFS from root to get traversal order
  build_children <- function(seg_set, parent_vec) {
    par   <- parent_vec[seg_set]
    valid <- !is.na(par) & par %in% seg_set
    if (!any(valid)) return(list())
    split(seg_set[valid], par[valid])
  }

  p_children  <- build_children(p_seg_ids, p_parent)
  bfs_order   <- character(0)
  queue       <- root_seg
  while (length(queue) > 0) {
    curr      <- queue[1]
    queue     <- queue[-1]
    bfs_order <- c(bfs_order, curr)
    kids      <- p_children[[curr]]
    if (!is.null(kids)) queue <- c(queue, kids)
  }

  # Cumulative distance from root to each node (includes own length)
  cum_dist <- stats::setNames(numeric(length(bfs_order)), bfs_order)
  cum_dist[root_seg] <- lengths[root_seg]
  for (node in bfs_order[-1]) {
    cum_dist[node] <- cum_dist[p_parent[node]] + lengths[node]
  }

  # --- 10. Compute pairwise distances ---
  # dist(A, B) = cum_dist[A] + cum_dist[B] - 2*cum_dist[LCA] + lengths[LCA]
  get_path <- function(start) {
    path    <- character(0)
    current <- start
    while (!is.na(current) && current %in% p_seg_ids) {
      path    <- c(path, current)
      if (current == root_seg) break
      current <- p_parent[current]
    }
    path
  }

  dist_mat <- matrix(0, nrow = n_sites, ncol = n_sites,
                     dimnames = list(sample_ids, sample_ids))

  for (i in seq_len(n_sites)) {
    for (j in seq_len(i - 1)) {
      ci <- site_comids[i]
      cj <- site_comids[j]
      if (ci == cj) {
        dist_mat[i, j] <- dist_mat[j, i] <- 0
        next
      }
      path_i  <- get_path(ci)
      path_j  <- get_path(cj)
      common  <- intersect(path_i, path_j)
      if (length(common) == 0) {
        warning("No common path found for ", sample_ids[i], " and ",
                sample_ids[j], ". Setting distance to NA.")
        dist_mat[i, j] <- dist_mat[j, i] <- NA
        next
      }
      lca_ij  <- common[which.max(cum_dist[common])]
      d       <- cum_dist[ci] + cum_dist[cj] - 2 * cum_dist[lca_ij] + lengths[lca_ij]
      dist_mat[i, j] <- dist_mat[j, i] <- d
    }
  }

  if (format == "dist") return(stats::as.dist(dist_mat))
  dist_mat
}
