

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Homo.sapiens")
# BiocManager::install("GenomicRanges")
 
library(Homo.sapiens)
library(dplyr)


# Code is from https://github.com/gwas-partitioning/bnmf-clustering, slightly adapted for the purposes of this work.
# snps is a vector of variants in the form "chr:pos" or "chr:pos:ref:alt"
# Returns a data frame with variant name, gene name and ENTREZID. 
# If a gene name was not found, defaults to variant name

query_locus_names <- function(snps){
  
  geneRanges <- function(db, column="ENTREZID"){
    g <- genes(db, columns=column)
    col <- mcols(g)[[column]]
    genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
    mcols(genes)[[column]] <- as.character(unlist(col))
    genes
  }
  
  gns = geneRanges(Homo.sapiens, column="SYMBOL")
  # need chr | start | end, with gene names for row names
  df_gns <- data.frame(gns)
  gtf.gene <- df_gns %>%
    mutate(chr = gsub("chr","",seqnames)) %>%
    dplyr::select(chr, start, end, SYMBOL)
  
  # need chr | start | end, with gene names for row names
  df_entrez <- geneRanges(Homo.sapiens, column="ENTREZID") %>%
    data.frame() %>%
    mutate(chr = gsub("chr","",seqnames)) %>%
    mutate(ChrPos = paste(chr,start,sep=":")) %>% 
    dplyr::select(ChrPos, ENTREZID)
  
  #' Convert from string to range
  #' 
  #' @param pos A vector of strings ex. chr1 2938302 2938329
  #' @param delim Delimiter for string splitting
  #' @param region Boolean of whether region or just one position
  #'
  #' @returns Dataframe of ranges
  #' 
  string2range <- function(pos, delim=' ', region=TRUE) {
    posp <- as.data.frame(do.call(rbind, strsplit(pos, delim)))
    posp[,1] <- posp[,1]
    posp[,2] <- as.numeric(as.character(posp[,2]))
    if(region) {
      posp[,3] <- as.numeric(as.character(posp[,3]))
    } else {
      posp[,3] <- posp[,2]
    }
    return(posp)
  }
  
  #' Convert from ranges to GRanges
  #' 
  #' @param df Dataframe with columns as sequence name, start, and end
  #' 
  #' @returns GRanges version 
  #' 
  range2GRanges <- function(df) {
    require(GenomicRanges)
    require(IRanges)
    gr <- GenomicRanges::GRanges(
      seqnames = df[,1],
      ranges=IRanges(start = df[,2], end = df[,3])
    )
    return(gr)
  }
  
  # convert SNPs to GRanges
  snps.ranges <- string2range(snps, delim=":", region=FALSE)
  snps.granges <- range2GRanges(snps.ranges)
  names(snps.granges) <- snps
  
  # convert genes to GRanges
  gtf.granges <- range2GRanges(gtf.gene)
  names(gtf.granges) <-  gtf.gene$SYMBOL  #gene.names
  
  hits <- GenomicRanges::nearest(snps.granges, gtf.granges)
  # make vector of SNPs to gene
  
  df_gns <- df_gns %>% 
    mutate(chr = gsub("chr","",seqnames)) %>%
    mutate(ChrPos = paste(chr,start,sep=":")) %>%
    dplyr::select(chr, start, end, ChrPos, SYMBOL) %>%
    merge(df_entrez, by="ChrPos")
  
  df_hits <- data.frame(gene=names(gtf.granges)[hits]) %>%
    mutate(variant = names(snps.granges)) %>%
    merge(df_gns[,c("SYMBOL","ENTREZID")], by.x="gene", by.y="SYMBOL") %>%
    filter(!duplicated(variant))
  
  # Make the duplicated genes have unique names.
  #If no gene is found for a variant, then the name of the variant is returned
  if (nrow(df_hits) > 0) {
    dup_genes <- df_hits$gene[!is.na(df_hits$gene)]
    duplicates <- data.frame(duplicated = table(dup_genes) > 1)
    duplicates$gene <- rownames(duplicates)
    
    for (r in 1:nrow(duplicates)) {
      row <- duplicates[r,]
      if (!is.na(row$duplicated) && row$duplicated) {
        idx <- which(df_hits$gene == row$gene)
        df_hits$gene[idx] <- paste0(row$gene, "_", seq_along(idx))
      }
    }
  }
  
  res <- data.frame(variant = snps) %>%
    left_join(df_hits, by="variant") %>%
    mutate(gene = dplyr::coalesce(gene, variant)) # coalese fills in empty genes w/ variant
  
  return(res)
}


## Function for making circle plots given a vector of weights and the variant/trait they correspond to.


cluster_circle_plot <- function(cluster_weights, nvariants, total_names, title, rotate = 0) {
  
  ### Inputs
  ## cluster_weights: a vector of concatenated variant and trait weights for a given component (cluster).
  ##                  The variant are assumed to be first, followed by the trait weights.
  ## nvariants: the number of variants in the data
  ## total_names: the corresponding names of the variant and trait weights,
  ## title: the plot's title
  ## rotate: used to rotate the circle if needed
  
  ### Output: the circle plot
  
  cluster_data = data.frame(
    names = total_names,
    group = c(ifelse(cluster_weights[1:nvariants] > 0, "Vpos", "Vneg"), ifelse(cluster_weights[(nvariants+1):length(cluster_weights)] > 0, "Tpos", "Tneg")),
    weights = cluster_weights*35
  )
  
  cluster_data = cluster_data[cluster_data$weights != 0, ]
  cluster_data = cluster_data %>% arrange(group, abs(weights))
  cluster_data = cluster_data %>% arrange(desc(abs(weights))) %>% head(30)
  
  empty_bar <- 4
  to_add <- data.frame(matrix(NA, empty_bar * nlevels(factor(cluster_data$group)), ncol(cluster_data)))
  colnames(to_add) <- colnames(cluster_data)
  to_add$group <- rep(levels(factor(cluster_data$group)), each=empty_bar)
  cluster_data <- rbind(cluster_data, to_add)
  cluster_data <- cluster_data %>% arrange(group)
  cluster_data$id <- seq(1, nrow(cluster_data))
  
  rot = 0.5
  p.rot = 14.8
  if (nrow(cluster_data) < 30) {
    rot = 0.2
    p.rot = 5.27
  }
  if (length(unique(cluster_data$group)) < 4 & nrow(cluster_data) >= 30) {
    rot = 0.2
    p.rot = 13.6
  }
  
  label_data <- cluster_data
  number_of_bar <- nrow(label_data)
  angle <- 90 - 360 * (label_data$id - rot) / number_of_bar - rotate
  label_data$hjust <- ifelse(angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle + 180, angle)
  
  group_colors <- c("Vpos" = "#228B22", "Vneg" = "#9370DB",
                    "Tpos" = "#F08080", "Tneg" = "#87CEEB")
  
  pl = ggplot(cluster_data, aes(x = as.factor(id), y = weights, fill = group)) +
    geom_bar(stat = "identity", alpha = 0.5) +
    ylim(-50, 50) +                  
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 40, face = "bold", vjust=-11),
      plot.margin = margin(t = -70, r = -60, b = -20, l = -100, unit = "pt")   # <-- increase top (t) if needed
    ) +
    coord_polar(start = rotate * pi/180) +
    scale_fill_manual(values = group_colors) +
    geom_text(data = label_data, 
              aes(x = id, 
                  y = ifelse(weights > 0, weights + 5, 5),
                  label = names, 
                  hjust = hjust), 
              color = "black", fontface = "bold", alpha = 0.6, 
              angle = label_data$angle, 
              inherit.aes = FALSE,
              size = 7) +
    annotate("path", x = seq(0, p.rot*pi, length.out = 200), y = rep(0, 200), 
             linewidth = 0.3, color = "black") +
    ggtitle(title)
  
  suppressWarnings(print(pl))
}


### Function to compute the jaccard matrix between two weight matrices


compute_jaccard_matrix <- function(A, B) {
  
  ### Inputs:
  ## A,B: the two weight matrices, where the weights of each component are assumed to be in the columns
  
  ### Output: a matrix with jaccard indices between the components of A and B, where the rows correspond to the
  ###         components of A and the columns correspond to the components of B
  
  k.1 = ncol(A)
  k.2 = ncol(B)
  
  jaccard.matrix = matrix(0, nrow = k.1, ncol = k.2)
  
  for (ii in 1:k.1) {
    for (jj in 1:k.2) {
      jaccard.matrix[ii, jj] <- sum(pmin(A[,ii], B[,jj])) / sum(pmax(A[,ii], B[,jj]))
    }
  }
  
  return(jaccard.matrix)
}




## Stacked barplot for contribution scores with base R (not used)


top_n_values_indices <- function(vec, n) {
  indices <- order(vec, decreasing = TRUE)[1:n]  # Get top n indices
  values <- vec[indices]  # Get corresponding values
  return(list(values = values, indices = indices))
}


stacked_barplot <- function(contr.scores, n, main = NULL, ylab = NULL, xlab = NULL) {
  
  k = ncol(contr.scores)
  
  data.matrix = matrix(0, nrow = n+1, ncol = k)
  
  for (ii in 1:k) {
    
    contr.lat = top_n_values_indices(contr.scores[,ii], n=n)
    data.matrix[,ii] = c(contr.lat$values, 1 - sum(contr.lat$values))
  }
  
  colors = c(colorRampPalette(c("darkblue", "blue", "cyan"))(n), "grey")
  
  default.par = par(no.readonly = TRUE)  
  par(mar = c(4, 4, 4, 2))  # Increase bottom margin to fit rotated labels
  
  bar_positions = barplot(data.matrix, beside = FALSE, col = colors, ylim = c(0, 1.05),
                          ylab = "", yaxt = "n", yaxs = "i", names.arg = rep("", k))
  
  axis(2, at = seq(0, 1, by = 0.2))
  mtext(ylab, side = 2, line = 2.5, adj = 0.475, cex = 1.4)
  title(main = main, cex.main = 2)
  
  labels <- paste0(xlab, 1:k)
  text(
    x = bar_positions,
    y = par("usr")[3] - 0.03,  # slightly below x-axis
    labels = labels,
    srt = 60,
    adj = 1,
    xpd = TRUE,
    cex = 0.9
  )
  
  par(default.par)
}




#### Modified pheatmap to support rownames from left 


library(grid)
library(gtable)

# Modified pheatmap:::heatmap_motor
heatmap_motor <- function (matrix, border_color, cellwidth, cellheight, tree_col, 
                           tree_row, treeheight_col, treeheight_row, filename, width, 
                           height, breaks, color, legend, annotation_row, annotation_col, 
                           annotation_colors, annotation_legend, annotation_names_row, 
                           annotation_names_col, main, fontsize, fontsize_row, fontsize_col, 
                           hjust_col, vjust_col, angle_col, fmat, fontsize_number, number_color, 
                           gaps_col, gaps_row, labels_row, labels_col, ...) 
{
  lo = pheatmap:::lo(coln = labels_col, rown = labels_row, nrow = nrow(matrix), 
                     ncol = ncol(matrix), cellwidth = cellwidth, cellheight = cellheight, 
                     treeheight_col = treeheight_col, treeheight_row = treeheight_row, 
                     legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, 
                     annotation_colors = annotation_colors, annotation_legend = annotation_legend, 
                     annotation_names_row = annotation_names_row, annotation_names_col = annotation_names_col, 
                     main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
                     fontsize_col = fontsize_col, angle_col = angle_col, gaps_row = gaps_row, 
                     gaps_col = gaps_col, ...)
  res = lo$gt
  mindim = lo$mindim
  if (!is.na(filename)) {
    if (is.na(height)) {
      height = convertHeight(gtable_height(res), "inches", valueOnly = T)
    }
    if (is.na(width)) {
      width = convertWidth(gtable_width(res), "inches", valueOnly = T)
    }
    r = regexpr("\\.[a-zA-Z]*$", filename)
    if (r == -1) 
      stop("Improper filename")
    ending = substr(filename, r + 1, r + attr(r, "match.length"))
    f = switch(ending, pdf = function(x, ...) pdf(x, ...), 
               png = function(x, ...) png(x, units = "in", res = 300, 
                                          ...), jpeg = function(x, ...) jpeg(x, units = "in", 
                                                                             res = 300, ...), jpg = function(x, ...) jpeg(x, 
                                                                                                                          units = "in", res = 300, ...), tiff = function(x, 
                                                                                                                                                                         ...) tiff(x, units = "in", res = 300, compression = "lzw", 
                                                                                                                                                                                   ...), bmp = function(x, ...) bmp(x, units = "in", 
                                                                                                                                                                                                                    res = 300, ...), stop("File type should be: pdf, png, bmp, jpg, tiff"))
    f(filename, height = height, width = width)
    gt = heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, 
                       border_color = border_color, tree_col = tree_col, 
                       tree_row = tree_row, treeheight_col = treeheight_col, 
                       treeheight_row = treeheight_row, breaks = breaks, 
                       color = color, legend = legend, annotation_col = annotation_col, 
                       annotation_row = annotation_row, annotation_colors = annotation_colors, 
                       annotation_legend = annotation_legend, annotation_names_row = annotation_names_row, 
                       annotation_names_col = annotation_names_col, filename = NA, 
                       main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
                       fontsize_col = fontsize_col, hjust_col = hjust_col, 
                       vjust_col = vjust_col, angle_col = angle_col, fmat = fmat, 
                       fontsize_number = fontsize_number, number_color = number_color, 
                       labels_row = labels_row, labels_col = labels_col, 
                       gaps_col = gaps_col, gaps_row = gaps_row, ...)
    grid.draw(gt)
    dev.off()
    return(gt)
  }
  if (mindim < 3) 
    border_color = NA
  if (!is.na(main)) {
    elem = pheatmap:::draw_main(main, fontsize = 1.3 * fontsize, ...)
    res = gtable_add_grob(res, elem, t = 1, l = 3, name = "main", 
                          clip = "off")
  }
  if (!pheatmap:::is.na2(tree_col) & treeheight_col != 0) {
    elem = pheatmap:::draw_dendrogram(tree_col, gaps_col, horizontal = T)
    res = gtable_add_grob(res, elem, t = 2, l = 3, name = "col_tree")
  }
  if (!pheatmap:::is.na2(tree_row) & treeheight_row != 0) {
    elem = pheatmap:::draw_dendrogram(tree_row, gaps_row, horizontal = F)
    res = gtable_add_grob(res, elem, t = 4, l = 1, name = "row_tree")
  }
  elem = pheatmap:::draw_matrix(matrix, border_color, gaps_row, gaps_col, 
                                fmat, fontsize_number, number_color)
  res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", 
                        name = "matrix")
  if (length(labels_col) != 0) {
    pars = list(labels_col, gaps = gaps_col, fontsize = fontsize_col, 
                hjust_col = hjust_col, vjust_col = vjust_col, angle_col = angle_col, 
                ...)
    elem = do.call(pheatmap:::draw_colnames, pars)
    res = gtable_add_grob(res, elem, t = 5, l = 3, clip = "off", 
                          name = "col_names")
  }
  if (length(labels_row) != 0) {
    pars = list(labels_row, gaps = gaps_row, fontsize = fontsize_row, 
                ...)
    elem = do.call(pheatmap:::draw_rownames, pars)
    res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", 
                          name = "row_names")
  }
  if (!pheatmap:::is.na2(annotation_col)) {
    converted_annotation = convert_annotations(annotation_col, 
                                               annotation_colors)
    elem = pheatmap:::draw_annotations(converted_annotation, border_color, 
                                       gaps_col, fontsize, horizontal = T)
    res = gtable_add_grob(res, elem, t = 3, l = 3, clip = "off", 
                          name = "col_annotation")
    if (annotation_names_col) {
      elem = pheatmap:::draw_annotation_names(annotation_col, fontsize, 
                                              horizontal = T)
      res = gtable_add_grob(res, elem, t = 3, l = 4, clip = "off", 
                            name = "col_annotation_names")
    }
  }
  if (!pheatmap:::is.na2(annotation_row)) {
    converted_annotation = convert_annotations(annotation_row, 
                                               annotation_colors)
    elem = pheatmap:::draw_annotations(converted_annotation, border_color, 
                                       gaps_row, fontsize, horizontal = F)
    res = gtable_add_grob(res, elem, t = 4, l = 2, clip = "off", 
                          name = "row_annotation")
    if (annotation_names_row) {
      elem = pheatmap:::draw_annotation_names(annotation_row, fontsize, 
                                              horizontal = F, hjust_col = hjust_col, vjust_col = vjust_col, 
                                              angle_col = angle_col)
      res = gtable_add_grob(res, elem, t = 5, l = 2, clip = "off", 
                            name = "row_annotation_names")
    }
  }
  annotation = c(annotation_col[length(annotation_col):1], 
                 annotation_row[length(annotation_row):1])
  annotation = annotation[unlist(lapply(annotation, function(x) !pheatmap:::is.na2(x)))]
  if (length(annotation) > 0 & annotation_legend) {
    elem = pheatmap:::draw_annotation_legend(annotation, annotation_colors, 
                                             border_color, fontsize = fontsize, ...)
    t = ifelse(is.null(labels_row), 4, 3)
    res = gtable_add_grob(res, elem, t = t, l = 6, b = 5, 
                          clip = "off", name = "annotation_legend")
  }
  if (!pheatmap:::is.na2(legend)) {
    elem = pheatmap:::draw_legend(color, breaks, legend, fontsize = fontsize, 
                                  ...)
    t = ifelse(is.null(labels_row), 4, 3)
    res = gtable_add_grob(res, elem, t = t, l = 5, b = 5, 
                          clip = "off", name = "legend")
  }
  return(res)
}

# Modified pheatmap:::lo    
lo <- function (rown, coln, nrow, ncol, cellheight = NA, cellwidth = NA, 
                treeheight_col, treeheight_row, legend, annotation_row, annotation_col, 
                annotation_colors, annotation_legend, annotation_names_row, 
                annotation_names_col, main, fontsize, fontsize_row, fontsize_col, 
                angle_col, gaps_row, gaps_col, ...) 
{
  if (!is.null(coln[1]) | (!pheatmap:::is.na2(annotation_row) & annotation_names_row)) {
    if (!is.null(coln[1])) {
      t = coln
    }
    else {
      t = ""
    }
    tw = strwidth(t, units = "in", cex = fontsize_col/fontsize)
    if (annotation_names_row) {
      t = c(t, colnames(annotation_row))
      tw = c(tw, strwidth(colnames(annotation_row), units = "in"))
    }
    longest_coln = which.max(tw)
    gp = list(fontsize = ifelse(longest_coln <= length(coln), 
                                fontsize_col, fontsize), ...)
    coln_height = unit(1, "grobheight", textGrob(t[longest_coln], 
                                                 rot = angle_col, gp = do.call(gpar, gp))) + unit(10, 
                                                                                                  "bigpts")
  }
  else {
    coln_height = unit(5, "bigpts")
  }
  if (!is.null(rown[1])) {
    t = rown
    tw = strwidth(t, units = "in", cex = fontsize_row/fontsize)
    if (annotation_names_col) {
      t = c(t, colnames(annotation_col))
      tw = c(tw, strwidth(colnames(annotation_col), units = "in"))
    }
    longest_rown = which.max(tw)
    gp = list(fontsize = ifelse(longest_rown <= length(rown), 
                                fontsize_row, fontsize), ...)
    rown_width = unit(1, "grobwidth", textGrob(t[longest_rown], 
                                               rot = 0, gp = do.call(gpar, gp))) + unit(10, "bigpts")
  }
  else {
    rown_width = unit(5, "bigpts")
  }
  gp = list(fontsize = fontsize, ...)
  if (!pheatmap:::is.na2(legend)) {
    longest_break = which.max(nchar(names(legend)))
    longest_break = unit(1.1, "grobwidth", 
                         textGrob(as.character(names(legend))[longest_break], 
                                  gp = do.call(gpar, gp)))
    title_length = unit(1.1, "grobwidth", textGrob("Scale", 
                                                   gp = gpar(fontface = "bold", ...)))
    legend_width = unit(12, "bigpts") + longest_break * 1.2
    legend_width = max(title_length, legend_width)
  }
  else {
    legend_width = unit(0, "bigpts")
  }
  if (is.na(main)) {
    main_height = unit(0, "npc")
  }
  else {
    main_height = unit(1.5, "grobheight", textGrob(main, 
                                                   gp = gpar(fontsize = 1.3 * fontsize, ...)))
  }
  textheight = unit(fontsize, "bigpts")
  if (!pheatmap:::is.na2(annotation_col)) {
    annot_col_height = ncol(annotation_col) * (textheight + 
                                                 unit(2, "bigpts")) + unit(2, "bigpts")
    t = c(as.vector(as.matrix(annotation_col)), colnames(annotation_col))
    annot_col_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], 
                                                             gp = gpar(...))) + unit(12, "bigpts")
    if (!annotation_legend) {
      annot_col_legend_width = unit(0, "npc")
    }
  }
  else {
    annot_col_height = unit(0, "bigpts")
    annot_col_legend_width = unit(0, "bigpts")
  }
  if (!pheatmap:::is.na2(annotation_row)) {
    annot_row_width = ncol(annotation_row) * (textheight + 
                                                unit(2, "bigpts")) + unit(2, "bigpts")
    t = c(as.vector(as.matrix(annotation_row)), colnames(annotation_row))
    annot_row_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], 
                                                             gp = gpar(...))) + unit(12, "bigpts")
    if (!annotation_legend) {
      annot_row_legend_width = unit(0, "npc")
    }
  }
  else {
    annot_row_width = unit(0, "bigpts")
    annot_row_legend_width = unit(0, "bigpts")
  }
  annot_legend_width = max(annot_row_legend_width, annot_col_legend_width)
  treeheight_col = unit(treeheight_col, "bigpts") + unit(5, 
                                                         "bigpts")
  treeheight_row = unit(treeheight_row, "bigpts") + unit(5, 
                                                         "bigpts")
  if (is.na(cellwidth)) {
    mat_width = unit(1, "npc") - rown_width - legend_width - 
      treeheight_row - annot_row_width - annot_legend_width
  }
  else {
    mat_width = unit(cellwidth * ncol, "bigpts") + length(gaps_col) * 
      unit(4, "bigpts")
  }
  if (is.na(cellheight)) {
    mat_height = unit(1, "npc") - main_height - coln_height - 
      treeheight_col - annot_col_height
  }
  else {
    mat_height = unit(cellheight * nrow, "bigpts") + length(gaps_row) * 
      unit(4, "bigpts")
  }
  gt = gtable(widths = unit.c(treeheight_row, rown_width,  
                              mat_width, treeheight_row, legend_width, annot_legend_width), 
              heights = unit.c(main_height, treeheight_col, annot_col_height, 
                               mat_height, coln_height), vp = viewport(gp = do.call(gpar, 
                                                                                    gp)))
  cw = convertWidth(mat_width - (length(gaps_col) * unit(4, 
                                                         "bigpts")), "bigpts", valueOnly = T)/ncol
  ch = convertHeight(mat_height - (length(gaps_row) * unit(4, 
                                                           "bigpts")), "bigpts", valueOnly = T)/nrow
  mindim = min(cw, ch)
  res = list(gt = gt, mindim = mindim)
  return(res)
}

# Modified pheatmap:::draw_rownames      
draw_rownames <- function (rown, gaps, ...) 
{
  coord = pheatmap:::find_coordinates(length(rown), gaps)
  y = unit(1, "npc") - (coord$coord - 0.5 * coord$size)
  res = textGrob(rown, x = unit(-3, "bigpts"), y = y, vjust = 0.5, 
                 hjust = 1, gp = gpar(...))
  return(res)
}

assignInNamespace(x="draw_rownames", value=draw_rownames, ns="pheatmap")
assignInNamespace(x="lo", value=lo, ns="pheatmap")
assignInNamespace(x="heatmap_motor", value=heatmap_motor, ns="pheatmap")






