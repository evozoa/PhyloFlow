#' Build a phylogenetic tree from sample collection sites
#'
#' The primary user-facing function for phylogeographic workflows. Accepts a
#' data frame of sample collection coordinates, finds the stream network
#' connecting them, and returns a rooted phylogenetic tree where each tip is a
#' sample site. Branch lengths reflect stream segment lengths.
#'
#' With automatic NHD data, each site is snapped to the nearest qualifying
#' stream and traced downstream along the main stem (via the NLDI) until all
#' sites meet, however far away that confluence is. Only the stream paths
#' connecting the sites are downloaded, so run time scales with the number
#' of sites rather than the size of the basin.
#'
#' @param sites A data frame with one row per sample site.
#' @param site_id_col Character. Column name for sample IDs. Default
#'   \code{"sample_id"}.
#' @param lat_col Character. Column name for latitude (decimal degrees, WGS84).
#'   Default \code{"lat"}.
#' @param lon_col Character. Column name for longitude (decimal degrees,
#'   WGS84). Default \code{"lon"}.
#' @param network An \code{sf} object of stream network line segments with
#'   NHDPlus attributes (\code{comid}, \code{hydroseq}, \code{dnhydroseq},
#'   \code{streamorde}). If \code{NULL}, NHD data is fetched automatically.
#' @param buffer_km Numeric. Search radius (km) around each site for the
#'   nearest qualifying stream. Default 10. Ignored for user-supplied
#'   networks.
#' @param min_stream_order Integer. Minimum NHD stream order to snap samples
#'   to. Default 3.
#' @param length_col Character. Column name for segment lengths. Default
#'   \code{"lengthkm"}. Set to \code{NULL} to compute from geometry in metres.
#' @param format Character. \code{"newick"} (default) or \code{"nexus"}.
#' @param file Character or \code{NULL}. Output file path. If \code{NULL} the
#'   tree is returned without writing to disk.
#' @param collapse_singles Logical. If \code{TRUE} (default), nodes with a
#'   single descendant (stream segments joined without a confluence) are
#'   removed and their branch lengths summed.
#'
#' @return An \code{ape} \code{phylo} object with tips labelled by sample ID.
#' @export
sites_to_tree <- function(sites,
                          site_id_col      = "sample_id",
                          lat_col          = "lat",
                          lon_col          = "lon",
                          network          = NULL,
                          buffer_km        = 10,
                          min_stream_order = 3,
                          length_col       = "lengthkm",
                          format           = "newick",
                          file             = NULL,
                          collapse_singles = TRUE) {

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

  # --- 2-7. Minimal network connecting all samples down to their LCA ---
  connected <- if (is.null(network)) {
    connect_sites_nhd(lats, lons, buffer_km, min_stream_order, length_col)
  } else {
    connect_sites_local(network, lats, lons, min_stream_order)
  }
  pruned_net      <- connected$pruned_net
  site_comids     <- connected$site_comids
  hydroseq_to_seg <- connected$hydroseq_to_seg
  names(site_comids)   <- sample_ids

  # --- 8. Build tree from minimal network ---
  build_children <- function(seg_set, parent_vec) {
    par   <- parent_vec[seg_set]
    valid <- !is.na(par) & par %in% seg_set
    if (!any(valid)) return(list())
    split(seg_set[valid], par[valid])
  }

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

  if (collapse_singles && n_tips >= 2)
    phylo_tree <- ape::collapse.singles(phylo_tree)

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

# Snap sites to NHD, trace each downstream via the NLDI, and keep the path
# segments above the first segment shared by all sites
connect_sites_nhd <- function(lats, lons, buffer_km, min_stream_order,
                              length_col) {
  n_sites <- length(lats)
  message("Snapping ", n_sites, " sites to the NHD network (",
          format_seconds(1.5 * n_sites), ")...")
  site_comids <- vapply(seq_len(n_sites), function(i) {
    snapped <- suppressMessages(
      snap_to_network(lats[i], lons[i], buffer_km = buffer_km,
                      min_stream_order = min_stream_order)
    )
    as.character(snapped$network$comid[snapped$segment_id])
  }, character(1))

  unique_comids <- unique(site_comids)
  message("Tracing ", length(unique_comids), " site(s) downstream (",
          format_seconds(2.5 * length(unique_comids)), ")...")
  paths <- lapply(unique_comids, function(comid) {
    dm <- nhdplusTools::navigate_nldi(
      list(featureSource = "comid", featureID = comid),
      mode = "DM", data_source = "flowlines", distance_km = 9999
    )
    unique(c(comid, as.character(dm$DM_flowlines$nhdplus_comid)))
  })

  common <- Reduce(intersect, paths)
  if (length(common) == 0)
    stop("The sites do not share a downstream outlet: they drain to the ",
         "sea (or to closed basins) separately, so no single stream tree ",
         "connects them.")

  all_ids <- unique(unlist(paths))
  message("Downloading ", length(all_ids), " connecting segments (",
          format_seconds(3 + 0.0085 * length(all_ids)), ")...")
  net <- nhdplusTools::get_nhdplus(
    comid = as.integer(all_ids), realization = "flowline",
    properties = unique(c("comid", "hydroseq", "dnhydroseq", "streamorde",
                          length_col)),
    skip_geometry = !is.null(length_col)
  )
  if (inherits(net, "sf")) net <- sf::st_transform(net, crs = 4326)

  seg_ids   <- as.character(net$comid)
  lca_comid <- common[which.max(net$hydroseq[match(common, seg_ids)])]
  message("Common outlet found: COMID ", lca_comid)

  keep <- setdiff(all_ids, setdiff(common, lca_comid))
  list(pruned_net      = net[seg_ids %in% keep, ],
       site_comids     = site_comids,
       hydroseq_to_seg = stats::setNames(seg_ids, as.character(net$hydroseq)))
}

# Snap sites to a user-supplied NHD-attributed network and keep each site's
# path down to the most upstream segment shared by all sites
connect_sites_local <- function(network, lats, lons, min_stream_order) {
  network <- sf::st_transform(network, crs = 4326)
  seg_ids <- as.character(network$comid)

  net_qual <- network[network$streamorde >= min_stream_order, ]
  if (nrow(net_qual) == 0)
    stop("No streams of order >= ", min_stream_order,
         " found. Try reducing min_stream_order.")

  site_comids <- character(length(lats))
  for (i in seq_along(lats)) {
    pt             <- sf::st_sf(geometry = sf::st_sfc(
      sf::st_point(c(lons[i], lats[i])), crs = 4326))
    nearest        <- sf::st_nearest_feature(pt, net_qual)
    site_comids[i] <- as.character(net_qual$comid[nearest])
  }
  message("All samples snapped to stream network.")

  hydroseq_to_seg   <- stats::setNames(seg_ids, as.character(network$hydroseq))
  parent_seg        <- hydroseq_to_seg[as.character(network$dnhydroseq)]
  names(parent_seg) <- seg_ids

  get_path_to_outlet <- function(start_comid) {
    path    <- character(0)
    current <- start_comid
    while (!is.na(current) && current %in% seg_ids && !current %in% path) {
      path    <- c(path, current)
      current <- parent_seg[[current]]
    }
    path
  }

  paths  <- lapply(unique(site_comids), get_path_to_outlet)
  common <- Reduce(intersect, paths)
  if (length(common) == 0)
    stop("No common ancestor found for all sample sites within the ",
         "supplied network.")

  hydroseq_vals <- network$hydroseq[match(common, seg_ids)]
  lca_comid     <- common[which.max(hydroseq_vals)]
  message("Common outlet found: COMID ", lca_comid)

  keep <- setdiff(unique(unlist(paths)), setdiff(common, lca_comid))
  list(pruned_net      = network[seg_ids %in% keep, ],
       site_comids     = site_comids,
       hydroseq_to_seg = hydroseq_to_seg)
}
