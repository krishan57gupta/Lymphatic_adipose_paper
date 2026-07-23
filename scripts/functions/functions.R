library(Seurat)
library(ggplot2)
library(ggrastr)
library(ggpubr)
library(export)
library(Hmisc)
library(stringr)
library(MASS)
library(scales)
library(ggbreak) 
library(cowplot)

log_breaks = function(maj, radix=10) {
  function(x) {
    minx         = floor(min(logb(x,radix), na.rm=T)) - 1
    maxx         = ceiling(max(logb(x,radix), na.rm=T)) + 1
    n_major      = maxx - minx + 1
    major_breaks = seq(minx, maxx, by=1)
    if (maj) {
      breaks = major_breaks
    } else {
      steps = logb(1:(radix-1),radix)
      breaks = rep(steps, times=n_major) +
        rep(major_breaks, each=radix-1)
    }
    radix^breaks
  }
}

replace_all <- function(df, pattern, replacement, ignore.case = FALSE, fixed = FALSE) {
  repl_col <- function(x) {
    if (is.factor(x)) {
      lv <- levels(x)
      levels(x) <- gsub(pattern, replacement, lv,
                        ignore.case = ignore.case, fixed = fixed)
      x
    } else if (is.character(x)) {
      ifelse(is.na(x), x, gsub(pattern, replacement, x,
                               ignore.case = ignore.case, fixed = fixed))
    } else x
  }
  as.data.frame(lapply(df, repl_col), stringsAsFactors = FALSE, check.names = FALSE)
}

replace_all_many <- function(df, patterns, ignore.case = FALSE, fixed = FALSE) {
  for (p in names(patterns)) {
    df <- replace_all(df, pattern = p, replacement = patterns[[p]],
                      ignore.case = ignore.case, fixed = fixed)
  }
  df
}


left_until_nth <- function(x, n = 10, delim = "/") {
  sapply(x, function(s) {
    p <- gregexpr(delim, s, fixed = TRUE)[[1]]
    if (length(p) >= n) substr(s, 1, p[n] - 1) else s
  }, USE.NAMES = FALSE)
}

redix=10
log_breaks = function(maj, radix=10) {
  function(x) {
    minx         = floor(min(logb(x,radix), na.rm=T)) - 1
    maxx         = ceiling(max(logb(x,radix), na.rm=T)) + 1
    n_major      = maxx - minx + 1
    major_breaks = seq(minx, maxx, by=1)
    if (maj) {
      breaks = major_breaks
    } else {
      steps = logb(1:(radix-1),radix)
      breaks = rep(steps, times=n_major) +
        rep(major_breaks, each=radix-1)
    }
    radix^breaks
  }
}

# Median Normalization
medianNorm<-function(data){ # row as genes and columns as cells
  # cell filtering cells with atleast 10% expressed genes
  data<-data[,Matrix::colSums(data>0)>200]
  # gene filtering genes with positive expression in atleast 3 cells
  data<-data[which(Matrix::rowSums(data > 0) > 1),]
  # Median normalization
  data=as.matrix(data)
  cells_sum<-Matrix::rowSums(t(data))
  data<-Matrix::t(t(data)/(cells_sum/stats::median(cells_sum)))
  return(data)
}
################################################################################

# Z score calculation
Zscore<-function(data){ # row as genes and columns as cells
  # cell filtering cells with atleast 10% expressed genes
  data<-data[,Matrix::colSums(data>0)>200]
  # gene filtering genes with positive expression in atleast 3 cells
  data<-data[which(Matrix::rowSums(data > 0) > 1),]
  data[data<1]=1
  # log 2 transformation
  data<-log2(data)
  # zscore calculation
  data<-t(apply(data,1,function(x) (x-mean(x))/sd(x)))
  # removing genes with having null counts due to zscore calculation
  data<-data[rowSums(is.na(data))==0,]
  return(data)
}
################################################################################

# Stouffer score calculation
stoufferScore<-function(g,data){# row as genes and columns as cells
  # stouffer score calculation
  data=apply(data[g,],2,function(x) sum(x)/sqrt(length(x)))
  return(data)
}
################################################################################

# Regression plot wth correlation. having two sets of genes as signature to caculate stouffer score
corrPlot<-function(data,condition,condName,type,ct,DT,markers,clust1,clust2,XL,YL,LX,LY){ # rows as genes and columns as cells
  mat=data[,data@meta.data[[condition]]==condName & data@meta.data[[type]]==ct]
  mat=attr(mat@assays$RNA,DT)
  mat<-medianNorm(mat)
  mat<-Zscore(mat)
  g1=markers$gene[markers$cluster==clust1]
  g2=markers$gene[markers$cluster==clust2]
  g12=intersect(g1,g2)
  g1=setdiff(g1,g12)
  g2=setdiff(g2,g12)
  g1=intersect(g1,rownames(mat))
  g2=intersect(g2,rownames(mat))
  x=stoufferScore(g1,mat)
  y=stoufferScore(g2,mat)
  df=data.frame("x"=x,"y"=y)
  plot=ggscatter(df,x = "x",y = "y",add = "reg.line", add.params = list(color="blue",fill = "lightgray"), conf.int = TRUE)+
    stat_cor(method = "pearson", label.x = LX, label.y = LY)+
    labs(x = XL, y = YL)
  plot=rasterize(plot, layers='Point', dpi=50)
  return(plot)
}

################################################################################

# Exact-ish equivalent to Seurat::DotPlot, but uses facet_grid(rows=...) when split.by is set
DotPlot_rows <- function(
    object,
    features,
    cols = c("lightgrey", "blue"),
    col.min = NA,
    col.max = NA,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    assay = NULL,
    slot = "data",
    scale = TRUE,
    scale.by = c("radius","size")[1],   # kept for API parity
    scale.min = NA,
    scale.max = NA
) {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required")
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix is required")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr is required")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr is required")
  
  # subset by idents if requested
  if (!is.null(idents)) object <- Seurat::subset(object, idents = idents)
  
  # assay handling
  if (is.null(assay)) assay <- Seurat::DefaultAssay(object)
  Seurat::DefaultAssay(object) <- assay
  
  # group (x-axis)
  if (is.null(group.by)) {
    id <- Seurat::Idents(object)
  } else {
    if (!group.by %in% colnames(object@meta.data)) stop("group.by not in metadata: ", group.by)
    id <- object@meta.data[[group.by]]
  }
  id <- factor(as.character(id))
  
  # split (facets)
  split <- NULL
  if (!is.null(split.by)) {
    if (!split.by %in% colnames(object@meta.data)) stop("split.by not in metadata: ", split.by)
    split <- factor(as.character(object@meta.data[[split.by]]))
  }
  
  # features present in this assay
  features <- unique(features)
  present  <- intersect(features, rownames(object[[assay]]))
  if (length(present) == 0) stop("None of the requested features are present in assay: ", assay)
  missing  <- setdiff(features, present)
  if (length(missing) > 0) message("Dropping missing features: ", paste(missing, collapse = ", "))
  
  # pull data matrix (handle SeuratObject v5 layer=)
  so_ver <- tryCatch(utils::packageVersion("SeuratObject"), error = function(...) numeric_version("0"))
  mat <- if (!is.na(so_ver) && so_ver >= "5.0.0") {
    Seurat::GetAssayData(object, assay = assay, layer = slot)[present, , drop = FALSE]
  } else {
    Seurat::GetAssayData(object, assay = assay, slot  = slot)[present, , drop = FALSE]
  }
  
  # cells x features long frame with id(+split)
  df <- as.data.frame(t(as.matrix(mat)), stringsAsFactors = FALSE)
  df$..id.. <- id
  if (!is.null(split)) df$..split.. <- split
  
  df_long <- tidyr::pivot_longer(
    df,
    cols = all_of(present),
    names_to = "features.plot",
    values_to = "expr"
  )
  df_long$id <- df_long$..id..
  if (!is.null(split)) df_long$split <- df_long$..split..
  df_long$..id.. <- NULL; df_long$..split.. <- NULL
  
  # summarise like Seurat::DotPlot (avg.exp & pct.exp)
  if (is.null(split)) {
    summ <- df_long |>
      dplyr::group_by(id, features.plot) |>
      dplyr::summarise(
        avg.exp = mean(expr, na.rm = TRUE),
        pct.exp = 100 * mean(expr > 0, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    summ <- df_long |>
      dplyr::group_by(id, features.plot, split) |>
      dplyr::summarise(
        avg.exp = mean(expr, na.rm = TRUE),
        pct.exp = 100 * mean(expr > 0, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  # scaling across groups per feature (and per split, if present)
  if (scale) {
    if (is.null(split)) {
      summ <- summ |>
        dplyr::group_by(features.plot) |>
        dplyr::mutate(avg.exp.scaled = as.numeric(scale(avg.exp))) |>
        dplyr::ungroup()
    } else {
      summ <- summ |>
        dplyr::group_by(features.plot, split) |>
        dplyr::mutate(avg.exp.scaled = as.numeric(scale(avg.exp))) |>
        dplyr::ungroup()
    }
    if (!is.na(scale.min)) summ$avg.exp.scaled[summ$avg.exp.scaled < scale.min] <- scale.min
    if (!is.na(scale.max)) summ$avg.exp.scaled[summ$avg.exp.scaled > scale.max] <- scale.max
    color_var <- "avg.exp.scaled"
  } else {
    summ$avg.exp.scaled <- summ$avg.exp
    color_var <- "avg.exp.scaled"
  }
  
  # optional color clipping
  if (!is.na(col.min)) summ[[color_var]][summ[[color_var]] < col.min] <- col.min
  if (!is.na(col.max)) summ[[color_var]][summ[[color_var]] > col.max] <- col.max
  
  # dot threshold (accept fraction <=1 or percent >=1 like Seurat)
  thr <- if (dot.min <= 1) dot.min * 100 else dot.min
  summ$pct.exp[ summ$pct.exp < thr ] <- NA_real_
  
  # ordering
  summ$features.plot <- factor(summ$features.plot, levels = rev(present))
  summ$id <- factor(summ$id)
  
  # plot (no deprecated facets!)
  p <- ggplot2::ggplot(summ, ggplot2::aes(x = id, y = features.plot)) +
    ggplot2::geom_point(
      ggplot2::aes(size = pct.exp, color = .data[[color_var]]),
      stroke = 0, na.rm = TRUE
    ) +
    ggplot2::scale_size(range = c(0, dot.scale), limits = c(0, 100), name = "% cells") +
    ggplot2::scale_color_gradientn(colors = cols, name = if (scale) "Avg exp (scaled)" else "Avg exp") +
    ggplot2::theme_classic() +
    ggplot2::labs(x = NULL, y = NULL)
  
  if (!is.null(split)) p <- p + ggplot2::facet_grid(rows = ggplot2::vars(split))
  
  return(p)
}


DotPlot_1 <- function (object, features, assay = NULL, cols = c("lightgrey", 
                                                   "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
          idents = NULL, group.by = NULL, split.by = NULL, cluster.idents = FALSE, 
          scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) 
{
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% 
                                                   rownames(x = brewer.pal.info))
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(X = 1:length(features), 
                                        FUN = function(x) {
                                          return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                        }))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.", call. = FALSE, 
              immediate. = TRUE)
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, cells = colnames(object[[assay]]), 
                                        idents = idents))
  data.features <- FetchData(object = object, vars = features, 
                             cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  }else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- FetchData(object = object, vars = split.by)[cells, 
                                                          split.by]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop(paste0("Need to specify at least ", length(x = unique(x = splits)), 
                    " colors using the cols parameter"))
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(what = rbind, args = lapply(X = data.plot, 
                                               FUN = unlist))
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled", 
            call. = FALSE, immediate. = TRUE)
  }else if (ngroup < 5 & scale) {
    warning("Scaling data with a low number of groups may produce misleading results", 
            call. = FALSE, immediate. = TRUE)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             if (scale) {
                               data.use <- scale(x = log1p(data.use))
                               data.use <- MinMax(data = data.use, min = col.min, 
                                                  max = col.max)
                             }
                             else {
                               data.use <- log1p(x = data.use)
                             }
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- unlist(x = lapply(X = data.plot$id, FUN = function(x) sub(paste0(".*_(", 
                                                                                   paste(sort(unique(x = splits), decreasing = TRUE), 
                                                                                         collapse = "|"), ")$"), "\\1", x)))
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = split.colors, yes = "colors", no = "avg.exp.scaled")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], 
                                       levels = unique(x = feature.groups))
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", 
                                                        y = "id")) + geom_point(mapping = aes_string(size = "pct.exp", 
                                                                                                     color = color.by)) + scale.func(range = c(0, dot.scale), 
                                                                                                                                     limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(), 
                                                                                                                                                                               axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) + 
    labs(x = "Features", y = ifelse(test = is.null(x = split.by), 
                                    yes = "Identity", no = "Split Identity")) + theme_cowplot()
  if (!is.null(x = feature.groups)) {
    plot <- plot + facet_grid(cols   = vars(feature.groups), scales = "free_x", 
                              space = "free_x", switch = "y") + theme(panel.spacing = unit(x = 1, 
                                                                                           units = "lines"), strip.background = element_blank())
    # plot <- plot + facet_grid(facets = ~feature.groups, scales = "free_x", 
    #                           space = "free_x", switch = "y") + theme(panel.spacing = unit(x = 1, 
    #                                                                                        units = "lines"), strip.background = element_blank())
  }
  if (split.colors) {
    plot <- plot + scale_color_identity()
  }else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  }else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (!split.colors) {
    plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
  }
  return(plot)
}



############################################################
## Circular pathway–gene plot (works with df OR cp object) ##
############################################################

circular_enrich_plot <- function(
    res,
    score_col = c("p.adjust", "pvalue", "qvalue"),
    cutoff = 0.05,
    top_n = 10,
    pathway_size_factor = 1,
    pathway_colors = NULL,
    gene_color = "grey30",
    gene_label_expand = 1.04,
    term_label_expand = 1.08,
    gene_point_size = 2.2
) {
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr)
    library(tidygraph); library(ggraph)
    library(scales)
  })
  
  get_enrich_table <- function(res) {
    if (inherits(res, "data.frame")) {
      res
    } else if (!is.null(res@result)) {
      res@result
    } else {
      stop("Need clusterProfiler result or data.frame with ID, Description, geneID.")
    }
  }
  
  score_col <- match.arg(score_col)
  df <- tibble::as_tibble(get_enrich_table(res))
  
  need <- c("ID", "Description", "geneID")
  miss <- setdiff(need, colnames(df))
  if (length(miss) > 0) stop("Missing: ", paste(miss, collapse = ", "))
  
  avail_scores <- intersect(c("p.adjust", "pvalue", "qvalue"), colnames(df))
  if (length(avail_scores) == 0) {
    term_df <- df %>% slice_head(n = top_n)
  } else {
    if (!(score_col %in% colnames(df))) {
      score_col <- avail_scores[1]
      message("Chosen score not found, using: ", score_col)
    }
    term_df <- df %>%
      dplyr::filter(.data[[score_col]] <= cutoff) %>%
      arrange(.data[[score_col]]) %>%
      slice_head(n = top_n)
    if (nrow(term_df) == 0) {
      term_df <- df %>%
        arrange(.data[[score_col]]) %>%
        slice_head(n = top_n)
    }
  }
  if (nrow(term_df) == 0) stop("No enriched terms to plot.")
  
  # ---- FIX: normalize geneID to character (handles list-column geneID) ----
  term_df <- term_df %>%
    mutate(
      ID = as.character(ID),
      Description = as.character(Description),
      geneID = if (is.list(geneID)) {
        vapply(geneID, function(x) paste(x, collapse = "/"), character(1))
      } else {
        as.character(geneID)
      }
    )
  
  # edges
  edge_df <- term_df %>%
    select(ID, Description, geneID) %>%
    tidyr::separate_rows(geneID, sep = "/") %>%
    mutate(geneID = as.character(geneID)) %>%
    filter(!is.na(geneID), geneID != "") %>%
    rename(term_id = ID, term_name = Description, gene = geneID) %>%
    tibble::as_tibble()
  
  # term nodes
  term_nodes <- edge_df %>%
    dplyr::count(term_id, term_name, name = "n_genes") %>%
    mutate(
      name    = term_name,
      is_term = TRUE,
      grp     = term_name
    )
  
  # gene nodes
  gene_nodes <- edge_df %>%
    distinct(gene) %>%
    mutate(
      term_id   = NA_character_,
      term_name = NA_character_,
      n_genes   = NA_integer_,
      name      = gene,
      is_term   = FALSE,
      grp       = "GENE"
    )
  
  nodes <- bind_rows(term_nodes, gene_nodes) %>%
    mutate(node_id = dplyr::row_number())
  
  edges <- edge_df %>%
    mutate(
      from = match(term_name, nodes$name),
      to   = match(gene,      nodes$name),
      term_col = term_name
    ) %>%
    select(from, to, term_col)
  
  if (any(is.na(edges$from)) || any(is.na(edges$to))) {
    stop("Some edges couldn't match nodes.")
  }
  
  g <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
  
  pathway_names <- term_nodes$term_name
  n_terms <- length(pathway_names)
  
  if (is.null(pathway_colors)) {
    pal_terms <- hue_pal()(n_terms)
  } else {
    pal_terms <- rep(pathway_colors, length.out = n_terms)
  }
  names(pal_terms) <- pathway_names
  
  pal_nodes <- c(pal_terms, GENE = gene_color)
  
  gene_counts <- term_nodes$n_genes
  min_gc <- min(gene_counts, na.rm = TRUE)
  max_gc <- max(gene_counts, na.rm = TRUE)
  brks <- pretty(c(min_gc, max_gc), n = 3)
  brks <- brks[brks >= min_gc & brks <= max_gc]
  if (length(brks) < 3) brks <- unique(c(min_gc, round((min_gc + max_gc)/2), max_gc))
  
  g <- g %>%
    activate(nodes) %>%
    mutate(grp = factor(grp, levels = names(pal_nodes)))
  
  p <- ggraph(g, layout = "linear", circular = TRUE) +
    geom_edge_arc(
      aes(edge_colour = term_col),
      strength = 0.35,
      edge_alpha = 0.6,
      edge_width = 0.6,
      show.legend = FALSE
    ) +
    geom_node_point(aes(
      x = x, y = y,
      colour = grp,
      size = ifelse(is_term, n_genes * pathway_size_factor, gene_point_size)
    )) +
    geom_node_text(
      data = function(d) {
        d %>%
          filter(!is_term) %>%
          mutate(
            angle = atan2(y, x) * 180 / pi,
            x = x * gene_label_expand,
            y = y * gene_label_expand,
            hjust = ifelse(angle > 90 | angle < -90, 1, 0),
            angle = ifelse(angle > 90 | angle < -90, angle + 180, angle)
          )
      },
      aes(x = x, y = y, label = name, angle = angle, hjust = hjust),
      size = 3,
      colour = "black"
    ) +
    geom_node_text(
      data = function(d) d %>% filter(is_term),
      aes(x = x * term_label_expand,
          y = y * term_label_expand,
          label = name,
          colour = grp),
      fontface = "bold",
      size = 3,
      show.legend = FALSE
    ) +
    scale_colour_manual(
      values = pal_nodes,
      breaks = names(pal_nodes),
      labels = names(pal_nodes),
      name   = "Pathways name"
    ) +
    scale_edge_colour_manual(values = pal_terms, guide = "none") +
    scale_size_continuous(
      range  = c(2, 5),
      name   = "Genes count in pathway",
      breaks = brks,
      labels = as.character(brks),
      guide  = guide_legend(override.aes = list(size = scales::rescale(brks, to = c(2, 5))))
    ) +
    coord_equal() +
    theme_void() +
    theme(legend.position = "right", plot.margin = margin(10, 10, 10, 10))
  
  print(p)
  return(list("plot"=p, "pathwaysColor"=pal_terms))
}


# circular_enrich_plot(res=dataGO3, score_col="pvalue", cutoff = 0.05, top_n=10, pathway_size_factor = 1,
#                      pathway_colors=NULL,gene_color = "black", # pathway_colors=c("red","green","blue","orange","yellow","cyan")
#                      gene_label_expand = 1.08, term_label_expand = 1.08, gene_point_size = 2)

update_dataGO2_genes <- function(
    cell_type,
    dataGO2,
    pathway_gene_list,
    sep = "/",
    update_count = TRUE,
    update_gene_ratio = TRUE
) {
  # ---- checks ----
  req <- c("cell.type", "Description", "geneID")
  miss <- setdiff(req, colnames(dataGO2))
  if (length(miss) > 0) stop("dataGO2 missing columns: ", paste(miss, collapse = ", "))
  
  # ---- normalize to avoid factor/name mismatches ----
  cell_type <- trimws(as.character(cell_type))
  dataGO2$cell.type   <- trimws(as.character(dataGO2$cell.type))
  dataGO2$Description <- trimws(as.character(dataGO2$Description))
  dataGO2$geneID      <- as.character(dataGO2$geneID)
  
  L <- pathway_gene_list[[cell_type]]
  if (is.null(L) || length(L) == 0) stop("No pathways found for this cell_type in pathway_gene_list.")
  names(L) <- trimws(names(L))
  L <- lapply(L, function(x) unique(trimws(as.character(x))))
  
  # ---- rows to update ----
  idx <- dataGO2$cell.type == cell_type & dataGO2$Description %in% names(L)
  if (!any(idx)) return(dataGO2)
  
  desc_vec <- dataGO2$Description[idx]
  
  # ---- update geneID ----
  dataGO2$geneID[idx] <- vapply(
    desc_vec,
    function(d) paste(L[[d]], collapse = sep),
    character(1)
  )
  
  # ---- optional: update Count ----
  if (update_count && "Count" %in% colnames(dataGO2)) {
    dataGO2$Count[idx] <- vapply(desc_vec, function(d) length(L[[d]]), integer(1))
  }
  
  # ---- optional: update GeneRatio ----
  if (update_gene_ratio && all(c("GeneRatio", "Count") %in% colnames(dataGO2))) {
    denom <- sub("^\\d+/(\\d+)$", "\\1", as.character(dataGO2$GeneRatio[idx]))
    dataGO2$GeneRatio[idx] <- paste0(dataGO2$Count[idx], "/", denom)
  }
  
  dataGO2
}

# run:
# dataGO2 <- update_dataGO2_genes(cell_type="BC", dataGO2=dataGO2, pathway_gene_list=up_in_AC_pathway_gene_list)

# ########################################################################################################################
# ######################### for searching TFs of gene, not completed yet, and very slow in running
# ########################################################################################################################
# 
# jaspar_predicted_tfs <- function(gene_symbol, species = c("human", "mouse"), score_cutoff = 400, flank_bp = 2000) {
#   species <- match.arg(species)
# 
#   # ---- packages ----
#   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# 
#   pkgs <- c("rtracklayer", "GenomicRanges", "AnnotationDbi")
#   for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)
# 
#   if (species == "human") {
#     if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE))
#       BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", ask = FALSE)
#     if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
#       BiocManager::install("org.Hs.eg.db", ask = FALSE)
# 
#     txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#     orgdb <- org.Hs.eg.db::org.Hs.eg.db
#     bb_url <- "https://hgdownload.soe.ucsc.edu/gbdb/hg38/jaspar/JASPAR2026.bb"
#   } else {
#     if (!requireNamespace("TxDb.Mmusculus.UCSC.mm10.knownGene", quietly = TRUE))
#       BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene", ask = FALSE)
#     if (!requireNamespace("org.Mm.eg.db", quietly = TRUE))
#       BiocManager::install("org.Mm.eg.db", ask = FALSE)
# 
#     txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
#     orgdb <- org.Mm.eg.db::org.Mm.eg.db
#     bb_url <- "https://hgdownload.soe.ucsc.edu/gbdb/mm10/jaspar/JASPAR2024.bb"
#   }
# 
#   # ---- map symbol -> ENTREZ ----
#   eid <- AnnotationDbi::mapIds(orgdb,
#                                keys = gene_symbol,
#                                keytype = "SYMBOL",
#                                column = "ENTREZID",
#                                multiVals = "first")
#   if (is.na(eid) || length(eid) == 0) stop("Could not map SYMBOL to ENTREZID: ", gene_symbol)
# 
#   # ---- gene range -> TSS -> promoter window ----
#   g <- GenomicFeatures::genes(txdb)[eid]
#   if (length(g) == 0) stop("Gene not found in TxDb for: ", gene_symbol)
# 
#   tss <- if (as.character(GenomicRanges::strand(g)) == "+") GenomicRanges::start(g) else GenomicRanges::end(g)
#   region <- GenomicRanges::GRanges(GenomicRanges::seqnames(g),
#                                    IRanges::IRanges(start = tss - flank_bp, end = tss + flank_bp))
# 
#   # ---- import TFBS bigBed only for this region ----
#   tfbs <- rtracklayer::import(rtracklayer::BigBedFile(bb_url), which = region)
#   if (length(tfbs) == 0) return(character())
# 
#   m <- as.data.frame(GenomicRanges::mcols(tfbs))
# 
#   # ---- optional score filter (if available) ----
#   if (!is.null(score_cutoff) && "score" %in% colnames(m)) {
#     keep <- !is.na(m$score) & (m$score >= score_cutoff)
#     tfbs <- tfbs[keep]
#     if (length(tfbs) == 0) return(character())
#     m <- as.data.frame(GenomicRanges::mcols(tfbs))
#   }
# 
#   # ---- find TF gene-symbol column robustly ----
#   # UCSC bigBeds usually include an extra TF-name field; try common names first.
#   cand <- c("TFName", "tfName", "tf_name", "geneName", "name2", "symbol")
#   tf_col <- intersect(cand, colnames(m))
# 
#   tf_names <- NULL
#   if (length(tf_col) >= 1) {
#     tf_names <- m[[tf_col[1]]]
#   } else if ("name" %in% colnames(m)) {
#     # fallback: sometimes "name" may look like "MAxxxx.x:TF" -> take part after ":"
#     nm <- as.character(m$name)
#     tf_names <- ifelse(grepl(":", nm), sub("^.*:", "", nm), nm)
#   } else {
#     stop("Could not detect TF-name column. Inspect: colnames(as.data.frame(mcols(tfbs))).")
#   }
# 
#   tf_names <- sort(unique(na.omit(as.character(tf_names))))
#   tf_names
# }
# 
# # Example:
# tfs_ldlr <- jaspar_predicted_tfs(gene_symbol="LDLR", species="human", score_cutoff = 400, flank_bp = 1000)
# length(tfs_ldlr); head(tfs_ldlr, 30)
# tfs_ldlr <- jaspar_predicted_tfs(gene_symbol="Ldlr", species="mouse", score_cutoff = 400, flank_bp = 1000)
# length(tfs_ldlr); head(tfs_ldlr, 30)
# ########################################################################################################################


############################################### corrected functions
assign_ligands_to_celltype_fixed <- function(seuratObj, ligands, celltype_col, 
                                             func.agg = mean, 
                                             func.assign = function(x) mean(x) + sd(x), 
                                             condition_oi = NULL, condition_col = NULL, ...) {
  if (any(!is.na(condition_col), !is.na(condition_oi)) & 
      !all(!is.na(condition_col), !is.na(condition_oi))) {
    stop("Please input both condition_colname and condition_oi")
  }
  if (any(!ligands %in% rownames(seuratObj))) {
    stop("Not all ligands are in the Seurat object")
  }
  
  slot <- "data"
  if (length(list(...)) > 0) {
    if (any(grepl("slot|layer", names(list(...))))) {
      slot <- list(...)[[which(grepl("slot|layer", names(list(...))))]]
    } else {
      warning("No slot/layer provided even though extra argument was provided, using default slot = 'data'")
    }
  }
  
  seuratObj_subset <- subset(seuratObj, features = ligands)
  
  if (!is.null(condition_oi)) {
    cond_vec <- seuratObj_subset[[condition_col, drop = TRUE]]
    seuratObj_subset <- seuratObj_subset[, cond_vec == condition_oi]
  }
  
  celltype_vec <- seuratObj_subset[[celltype_col, drop = TRUE]]  # fixed: drop = TRUE
  
  avg_expression_ligands <- lapply(unique(celltype_vec), function(celltype) {
    cells_this_type <- celltype_vec == celltype  # now a plain logical vector
    if (slot == "data") {
      expm1(GetAssayData(seuratObj_subset[, cells_this_type], ...)) %>% apply(1, func.agg)
    } else {
      apply(GetAssayData(seuratObj_subset[, cells_this_type], ...), 1, func.agg)
    }
  }) %>% setNames(unique(celltype_vec)) %>% do.call(cbind, .)
  
  sender_ligand_assignment <- avg_expression_ligands %>% 
    apply(1, function(ligand_expression) ligand_expression > func.assign(ligand_expression)) %>% 
    t()
  
  sender_ligand_assignment <- sender_ligand_assignment %>% 
    apply(2, function(x) x[x == TRUE]) %>% 
    purrr::keep(function(x) length(x) > 0)
  
  all_assigned_ligands <- sender_ligand_assignment %>% lapply(names) %>% unlist()
  unique_ligands <- all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
  general_ligands <- ligands %>% setdiff(unique_ligands)
  
  ligand_type_indication_df <- lapply(names(sender_ligand_assignment), function(sender) {
    unique_ligands_sender <- names(sender_ligand_assignment[[sender]]) %>% setdiff(general_ligands)
    if (length(unique_ligands_sender) > 0) {
      return(data.frame(ligand_type = sender, ligand = unique_ligands_sender))
    }
  }) %>% bind_rows()
  
  ligand_type_indication_df <- bind_rows(
    ligand_type_indication_df, 
    data.frame(ligand = general_ligands) %>% mutate(ligand_type = "General")
  )
  
  return(ligand_type_indication_df)
}


