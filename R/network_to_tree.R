#' Convert an upstream stream network to a phylogenetic tree
#'
#' Converts a directed stream network (as returned by \code{extract_upstream})
#' into a rooted phylogenetic tree. Each stream segment becomes a branch, with
#' length equal to segment length. For NHD data, topology is derived from
#' \code{hydroseq}/\code{dnhydroseq}, which are consistent across VPU
#' boundaries. For user-supplied networks, topology is derived from node ID
#' columns.
#'
#' @param upstream_network An \code{sf} object of upstream stream segments as
#'   returned by \code{extract_upstream}.
#' @param segment_col Character. Column uniquely identifying each segment.
#'   Default \code{"comid"} for NHD data.
#' @param length_col Character or \code{NULL}. Column for segment length. If
#'   \code{NULL}, lengths are computed from geometry in metres. For NHD data
#'   use \code{"lengthkm"} for lengths in kilometres.
#' @param from_col Character. Upstream node column for non-NHD networks.
#'   Default \code{"fromnode"}.
#' @param to_col Character. Downstream node column for non-NHD networks.
#'   Default \code{"tonode"}.
#' @param format Character. \code{"newick"} (default) or \code{"nexus"}.
#' @param file Character or \code{NULL}. Output file path. If \code{NULL} the
#'   tree is returned without writing to disk.
#'
#' @return An \code{ape} \code{phylo} object.
#' @export
network_to_tree <- function(upstream_network,
                            segment_col = "comid",
                            length_col = NULL,
                            from_col = "fromnode", to_col = "tonode",
                            format = "newick", file = NULL) {

  # --- Branch lengths ---
  if (is.null(length_col)) {
    net_m   <- sf::st_transform(upstream_network, crs = 3857)
    lengths <- as.numeric(sf::st_length(net_m))
  } else {
    lengths <- upstream_network[[length_col]]
  }

  use_nhd <- all(c("hydroseq", "dnhydroseq", segment_col) %in% names(upstream_network))

  if (use_nhd) {
    # --- NHD path: COMID-based topology via hydroseq ---
    seg_ids <- as.character(upstream_network[[segment_col]])
    names(lengths) <- seg_ids

    hydroseq_to_seg <- stats::setNames(seg_ids,
                                       as.character(upstream_network$hydroseq))
    parent_seg <- hydroseq_to_seg[as.character(upstream_network$dnhydroseq)]
    names(parent_seg) <- seg_ids

    # Root = most downstream segment (lowest hydroseq) with no parent in network
    has_parent   <- !is.na(parent_seg) & parent_seg %in% seg_ids
    root_candidates <- seg_ids[!has_parent]
    if (length(root_candidates) == 0)
      stop("No root segment found. The network may contain a cycle.")
    root_seg <- root_candidates[
      which.min(upstream_network$hydroseq[!has_parent])
    ]

    # BFS from root upstream — keep only segments reachable from root
    children_of <- split(seg_ids[has_parent],
                         parent_seg[has_parent])
    reachable <- character(0)
    queue     <- root_seg
    while (length(queue) > 0) {
      current   <- queue[1]
      queue     <- queue[-1]
      reachable <- c(reachable, current)
      kids      <- children_of[[current]]
      if (!is.null(kids))
        queue <- c(queue, kids[!kids %in% reachable])
    }

    # Filter to connected segments only
    keep       <- seg_ids %in% reachable
    seg_ids    <- seg_ids[keep]
    lengths    <- lengths[keep]
    parent_seg <- parent_seg[keep]
    names(parent_seg) <- seg_ids
    has_parent <- !is.na(parent_seg) & parent_seg %in% seg_ids

    tip_segs      <- seg_ids[!seg_ids %in% parent_seg[has_parent]]
    internal_segs <- c(root_seg,
                       setdiff(seg_ids[seg_ids %in% parent_seg[has_parent]],
                               root_seg))

    n_tips     <- length(tip_segs)
    n_internal <- length(internal_segs)

    node_idx <- c(
      stats::setNames(seq_len(n_tips), tip_segs),
      stats::setNames(n_tips + seq_len(n_internal), internal_segs)
    )

    non_root        <- seg_ids[has_parent]
    parents         <- parent_seg[non_root]
    edge_mat        <- matrix(0L, nrow = length(non_root), ncol = 2)
    edge_mat[, 1]   <- node_idx[parents]
    edge_mat[, 2]   <- node_idx[non_root]
    edge_len        <- lengths[non_root]

    tip_labels <- tip_segs

  } else {
    # --- User-supplied network: junction-node topology ---
    if (!from_col %in% names(upstream_network) ||
        !to_col   %in% names(upstream_network)) {
      stop(
        "Columns '", from_col, "' and '", to_col, "' not found.\n",
        "Available columns: ", paste(names(upstream_network), collapse = ", ")
      )
    }

    from_nodes <- as.character(upstream_network[[from_col]])
    to_nodes   <- as.character(upstream_network[[to_col]])

    edges_df <- data.frame(from = from_nodes, to = to_nodes,
                           length = lengths, stringsAsFactors = FALSE)
    edges_df <- edges_df[!is.na(edges_df$from) & !is.na(edges_df$to), ]

    g        <- igraph::graph_from_data_frame(edges_df, directed = TRUE)
    out_deg  <- igraph::degree(g, mode = "out")
    in_deg   <- igraph::degree(g, mode = "in")

    root_name <- names(which(out_deg == min(out_deg)))[1]
    tip_names <- names(which(in_deg == 0))

    all_nodes      <- igraph::V(g)$name
    internal_nodes <- c(root_name,
                        setdiff(all_nodes[!all_nodes %in% tip_names], root_name))

    n_tips     <- length(tip_names)
    n_internal <- length(internal_nodes)

    node_idx <- c(
      stats::setNames(seq_len(n_tips), tip_names),
      stats::setNames(n_tips + seq_len(n_internal), internal_nodes)
    )

    edf      <- igraph::as_data_frame(g, what = "edges")
    valid    <- edf$from %in% names(node_idx) & edf$to %in% names(node_idx)
    edf      <- edf[valid, ]

    edge_mat        <- matrix(0L, nrow = nrow(edf), ncol = 2)
    edge_mat[, 1]   <- node_idx[edf$to]     # parent = downstream junction
    edge_mat[, 2]   <- node_idx[edf$from]   # child  = upstream junction
    edge_len        <- edf$length

    tip_labels <- tip_names
  }

  # --- Validate ---
  if (any(is.na(edge_mat)) || any(edge_mat == 0L)) {
    stop("Invalid tree structure: edge matrix contains NA or zero values. ",
         "Check network topology.")
  }

  phylo_tree <- structure(
    list(edge        = edge_mat,
         edge.length = edge_len,
         tip.label   = tip_labels,
         Nnode       = as.integer(n_internal)),
    class = "phylo"
  )

  if (!is.null(file)) {
    if (format == "newick")      ape::write.tree(phylo_tree, file = file)
    else if (format == "nexus")  ape::write.nexus(phylo_tree, file = file)
    else stop("format must be 'newick' or 'nexus'")
    message("Tree written to: ", file)
    invisible(phylo_tree)
  } else {
    phylo_tree
  }
}
