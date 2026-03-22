#' Build a phylogenetic tree from sample collection sites
#'
#' The primary user-facing function for phylogeographic workflows. Accepts a
#' data frame of sample collection coordinates, downloads the connecting NHD
#' stream network, and returns a rooted phylogenetic tree where each tip is a
#' sample site. Branch lengths reflect stream segment lengths.
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
#'   used for NHD download. Increase if sites span a large area or if the
#'   common outlet is not found. Default 50.
#' @param min_stream_order Integer. Minimum NHD stream order to snap samples
#'   to. Default 3.
#' @param length_col Character. Column name for segment lengths. Default
#'   \code{"lengthkm"}. Set to \code{NULL} to compute from geometry in metres.
#' @param format Character. \code{"newick"} (default) or \code{"nexus"}.
#' @param file Character or \code{NULL}. Output file path. If \code{NULL} the
#'   tree is returned without writing to disk.
#'
#' @return An \code{ape} \code{phylo} object with tips labelled by sample ID.
#' @export
sites_to_tree <- function(sites,
                          site_id_col      = "sample_id",
                          lat_col          = "lat",
                          lon_col          = "lon",
                          network          = NULL,
                          buffer_km        = 50,
                          min_stream_order = 3,
                          length_col       = "lengthkm",
                          format           = "newick",
                          file             = NULL) {

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
  # Keep only COMIDs on the path from each sample downstream to the LCA.
  # This is the minimal Steiner tree connecting all samples via the stream network.
  keep_set <- character(0)
  for (sc in unique_sample_comids) {
    current <- sc
    while (!is.na(current) && current %in% up_seg_ids) {
      keep_set <- c(keep_set, current)
      current  <- up_parent[current]
    }
  }
  keep_set <- unique(keep_set)

  # --- 8. Build tree from minimal network ---
  build_children <- function(seg_set, parent_vec) {
    par   <- parent_vec[seg_set]
    valid <- !is.na(par) & par %in% seg_set
    if (!any(valid)) return(list())
    split(seg_set[valid], par[valid])
  }

  pruned_net <- upstream_net[upstream_net$comid %in% keep_set, ]
  p_seg_ids  <- as.character(pruned_net$comid)

  if (!is.null(length_col) && length_col %in% names(pruned_net)) {
    lengths <- stats::setNames(pruned_net[[length_col]], p_seg_ids)
  } else {
    net_m   <- sf::st_transform(pruned_net, crs = 3857)
    lengths <- stats::setNames(as.numeric(sf::st_length(net_m)), p_seg_ids)
  }

  p_parent     <- stats::setNames(
    hydroseq_to_seg[as.character(pruned_net$dnhydroseq)],
    p_seg_ids
  )
  p_has_parent <- !is.na(p_parent) & p_parent %in% p_seg_ids

  root_candidates <- p_seg_ids[!p_has_parent]
  if (length(root_candidates) == 0)
    stop("No root found in pruned network.")
  root_seg <- root_candidates[which.min(pruned_net$hydroseq[!p_has_parent])]

  # BFS from root to ensure connectivity
  p_children <- build_children(p_seg_ids, p_parent)
  reachable  <- character(0)
  queue      <- root_seg
  while (length(queue) > 0) {
    curr      <- queue[1]
    queue     <- queue[-1]
    reachable <- c(reachable, curr)
    kids      <- p_children[[curr]]
    if (!is.null(kids))
      queue <- c(queue, kids[!kids %in% reachable])
  }

  pruned_net   <- pruned_net[pruned_net$comid %in% reachable, ]
  p_seg_ids    <- as.character(pruned_net$comid)
  lengths      <- lengths[p_seg_ids]
  p_parent     <- p_parent[p_seg_ids]
  p_has_parent <- !is.na(p_parent) & p_parent %in% p_seg_ids
  p_children   <- build_children(p_seg_ids, p_parent)

  has_children  <- p_seg_ids %in% names(p_children)
  tip_segs      <- p_seg_ids[!has_children]
  internal_segs <- c(root_seg, setdiff(p_seg_ids[has_children], root_seg))
  n_tips        <- length(tip_segs)
  n_internal    <- length(internal_segs)

  # --- 9. Label tips with sample IDs ---
  comid_to_sample <- stats::setNames(sample_ids, site_comids)
  tip_labels <- ifelse(
    tip_segs %in% names(comid_to_sample),
    comid_to_sample[tip_segs],
    tip_segs
  )

  unmapped <- tip_segs[!tip_segs %in% names(comid_to_sample)]
  if (length(unmapped) > 0)
    warning(length(unmapped), " tip(s) labelled by COMID (not matched to a sample ID).")

  nested <- sample_ids[!site_comids %in% tip_segs]
  if (length(nested) > 0)
    warning("The following samples are internal nodes (downstream of other ",
            "samples): ", paste(nested, collapse = ", "),
            "\nThey represent confluence points in the tree rather than tips.")

  # --- 10. Assemble phylo object ---
  node_idx <- c(
    stats::setNames(seq_len(n_tips), tip_segs),
    stats::setNames(n_tips + seq_len(n_internal), internal_segs)
  )

  non_root      <- p_seg_ids[p_has_parent]
  parents       <- p_parent[non_root]
  valid         <- parents %in% names(node_idx) & non_root %in% names(node_idx)
  non_root      <- non_root[valid]
  parents       <- parents[valid]

  edge_mat        <- matrix(0L, nrow = length(non_root), ncol = 2)
  edge_mat[, 1]   <- node_idx[parents]
  edge_mat[, 2]   <- node_idx[non_root]
  edge_len        <- lengths[non_root]

  phylo_tree <- structure(
    list(edge        = edge_mat,
         edge.length = edge_len,
         tip.label   = tip_labels,
         Nnode       = as.integer(n_internal)),
    class = "phylo"
  )

  if (!is.null(file)) {
    if (format == "newick")     ape::write.tree(phylo_tree, file = file)
    else if (format == "nexus") ape::write.nexus(phylo_tree, file = file)
    else stop("format must be 'newick' or 'nexus'")
    message("Tree written to: ", file)
    invisible(phylo_tree)
  } else {
    phylo_tree
  }
}
