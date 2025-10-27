#' Merge two FraserDataSet objects with different splice sites and junctions
#'
#' This function merges two FraserDataSet objects that have different
#' splice sites and junctions but no overlapping samples. It combines
#' only the rawCountsJ and rawCountsSS assays and recalculates startID
#' and endID for the merged rowRanges.
#'
#' @param fds1 First FraserDataSet object.
#' @param fds2 Second FraserDataSet object.
#' @param join_type Type of join for rowRanges: "outer" (default) or "inner".
#' @param fds_name Name for the merged FraserDataSet.
#' @return A merged FraserDataSet object containing combined samples and assays.
#' @examples
#' # merged_fds <- mergeFDS(fds1, fds2, join_type="outer")
#' @export
#' @importFrom data.table as.data.table setorder setnames fcoalesce
#' @importFrom GenomeInfoDb seqlevelsInUse
#' @importFrom SummarizedExperiment rowRanges rowData assays
mergeFDS <- function(fds1, fds2, join_type = c("outer", "inner"), fds_name="merged_fds") {
  join_type <- match.arg(join_type)
  
  create_junction_counts <- function(fds){
    chr_levels <- seqlevelsInUse(fds)
    
    row_ranges <- as.data.table(rowRanges(fds))
    
    junction_counts <- as.data.table(assays(fds)$rawCountsJ)
    junction_counts <- cbind(junction_counts, row_ranges[, c("seqnames", "start", "end", "width", "strand", "startID", "endID")])
    # Enforce order
    junction_counts[, seqnames := factor(as.character(seqnames), levels = chr_levels)]
    setorder(junction_counts, seqnames, start)
    return (junction_counts)
  }
  
  create_split_counts <- function(fds){
    
    chr_levels <- seqlevelsInUse(fds)
    
    row_ranges <- as.data.table(rowRanges(fds))
    
    splice_sites <- rowData(nonSplicedReads(fds))
    splice_site_counts <- as.data.table(assays(fds)$rawCountsSS)
    splice_site_counts <- as.data.table(cbind(splice_sites, splice_site_counts))
    
    splice_site_lookup <- unique(rbind(
      row_ranges[, .(spliceSiteID = startID, seqnames, position = start - 1)],  # donor positions
      row_ranges[, .(spliceSiteID = endID, seqnames, position = end)]      # acceptor positions
    )) 
    
    splice_site_lookup[, start := position]
    splice_site_lookup[, end := start + 1]
    splice_site_lookup[, width := 2]
    splice_site_lookup[, position := NULL]
    
    splice_site_counts <- merge(splice_site_counts, splice_site_lookup, by = "spliceSiteID", all.x = TRUE)        
    splice_site_counts[, seqnames := factor(as.character(seqnames), levels = chr_levels)]
    setorder(splice_site_counts, seqnames, start)
    
    return (splice_site_counts)
    
  }
  
  
  junction_counts_1 <- create_junction_counts(fds1)
  junction_counts_2 <- create_junction_counts(fds2)
  
  
  
  junction_counts_merged <- merge(junction_counts_1, junction_counts_2, by=c("start", "end", "seqnames", "width", "strand"), all=T, )
  junction_counts_merged[is.na(junction_counts_merged)] <- 0
  
  positions <- unique(rbind(
    junction_counts_merged[, .(seqnames, pos = start)],
    junction_counts_merged[, .(seqnames, pos = end)]
  ))
  setorder(positions, seqnames, pos)
  positions[, spliceSiteID := .I]
  
  
  junction_counts_merged <- merge(junction_counts_merged, positions, by.x = c("seqnames", "start"), by.y = c("seqnames", "pos"), all.x = TRUE)
  setnames(junction_counts_merged, "spliceSiteID", "startID")
  
  junction_counts_merged <- merge(junction_counts_merged, positions, by.x = c("seqnames", "end"), by.y = c("seqnames", "pos"), all.x = TRUE)
  setnames(junction_counts_merged, "spliceSiteID", "endID")
  
  junction_counts_merged[ ,c("startID.x","startID.y", "endID.x", "endID.y") := NULL]
  
  splice_site_counts_1 <- create_split_counts(fds1)
  splice_site_counts_2 <- create_split_counts(fds2)
  
  splice_site_counts_merged <- merge(splice_site_counts_1, splice_site_counts_2, by=c("start", "end", "seqnames", "width", "type"), all=T, )
  splice_site_counts_merged[ ,c("spliceSiteID.x","spliceSiteID.y") := NULL]
  splice_site_counts_merged[is.na(splice_site_counts_merged)] <- 0
  
  
  splice_site_counts_merged <- merge(splice_site_counts_merged, positions, by.x = c("seqnames", "start"), by.y = c("seqnames", "pos"), all.x = TRUE)
  positions[, pos := pos - 1] # to merge Donnor values 
  splice_site_counts_merged <- merge(splice_site_counts_merged, positions, by.x = c("seqnames", "start"), by.y = c("seqnames", "pos"), all.x = TRUE)
  
  splice_site_counts_merged[, spliceSiteID := fcoalesce(spliceSiteID.x, spliceSiteID.y)]
  splice_site_counts_merged[ ,c("spliceSiteID.x","spliceSiteID.y") := NULL]
  
  
  # Subset both colData
  # --- Combine colData ---
  # Get intersection of column names
  shared_cols <- intersect(colnames(colData(fds1)), colnames(colData(fds2)))
  
  cd1 <- colData(fds1)[, shared_cols, drop = FALSE]
  cd2 <- colData(fds2)[, shared_cols, drop = FALSE]
  
  merged_colData <- rbind(cd1, cd2)  
  
  # --- Create merged FraserDataSet ---
  merged_fds <- FraserDataSet(
    colData = as.data.table(merged_colData),
    junctions = junction_counts_merged,
    spliceSites = splice_site_counts_merged,
    name = paste0(fds_name)
  )
  
  merged_fds <- calculatePSIValues(merged_fds, types=fitMetrics(merged_fds))
  
  
  return(merged_fds)
}