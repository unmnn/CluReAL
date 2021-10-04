#' Draw ideogram for clustering solution
#'
#' @inheritParams refine
#' @param similar_density_th Threshold of the (relative) difference between the
#' cluster with the highest density and the cluster with the lowest density for
#' the clustering to exhibit inter-cluster density differences.
#' @param radii_ratio_th Threshold of the extended-to-core ratio for the 
#' classification as "long-tailed" clusters.
#'
#' @export
#'
#' @examples
#' # TODO
draw_symbol <- function(cc, gv, rc,
                        similar_density_th = 3,
                        radii_ratio_th = 2) {
  # inputs
  #   cc: cluster context (ClusterContext)
  #   gv: goi validity indices (GValidty) 
  #   rc: cluster refinement context (RefinementContext)
  
  k <- cc$k
  outliers <- cc$outliers
  g_str <- gv$g_str
  g_rex <- gv$g_rex
  g_min <- gv$g_min
  vol_r <- mean(gv$vol_r)
  mm <- rc$mm
  kinship <- rc$kinship
  k_dens <- rc$k_dens
  
  child <- which(kinship == 1)
  dens_diff <- abs(max(k_dens) - min(k_dens)) / abs(min(max(k_dens), min(k_dens))) # !! ???
  # densdiff = np.absolute(np.nanmax(kdens)-np.nanmin(kdens))/np.absolute(np.minimum( np.nanmax(kdens),np.nanmin(kdens) ))  
  
  x_ec <- x_ec2 <- x_cc <- x_cc2 <- x_ccup <- x_ech <- y_ec <- y_cc <- 0
  v_ec <- v_ec2 <- v_cc <- v_cc2 <- v_ccup <- v_ech <- v_r1 <- v_l1 <- 
    v_eov <- v_ol <- v_om <- v_oh <- FALSE
  c_ec2 <- "black" 
  f_ec2 <- "transparent"
  
  if(any(mm)) {
    y_cc <- -0.08
    v_ccup <- TRUE
  } else {
    y_cc <- 0
  }
  
  if(length(child) > 0) {
    v_ech <- TRUE
  }
  
  if(g_min < 0 & k > 2 & g_str >= 0 & g_rex >= 1) {
    v_eov <- TRUE
  }
  
  if((k == 1) | (k == 2 & length(child) > 0)){
    x_ec <- y_ec <- 0
    v_ec <- v_cc <- TRUE
  } else {
    if (dens_diff >= similar_density_th) {
      c_ec2 <- f_ec2 <- "lightgray"
    }
    if(g_str > 1) {
      x_ec <- x_cc <- x_ccup <- -0.5
      x_ech <- -0.7
      v_ec <- v_cc <- TRUE
      x_ec2 <- x_cc2 <- 0.5
      v_ec2 <- v_cc2 <- TRUE
    } else if(g_str > 0) {
      if(g_rex > 1) {
        x_ec <- x_cc <- x_ccup <- -0.4
        x_ech <- -0.6
        v_ec <- v_cc <- TRUE
        x_ec2 <- x_cc2 <- 0.4
        v_ec2 <- v_cc2 <- TRUE
      } else {
        x_ec <- x_cc <- x_ccup <- -0.3
        x_ech <- -0.5
        v_ec <- v_cc <- TRUE
        x_ec2 <- x_cc2 <- 0.3
        v_ec2 <- v_cc2 <- TRUE
      }
    } else {
      if(g_rex > 1) {
        x_ec <- x_cc <- x_ccup <- -0.15
        x_ech <- -0.35
        v_ec <- v_cc <- TRUE
        x_ec2 <- x_cc2 <- 0.15
        v_ec2 <- v_cc2 <- TRUE
      } else if(g_rex > 0) {
        v_r1 <- TRUE
        x_cc <- x_ccup <- x_ech <- -0.2
        v_cc <- TRUE
        x_cc2 <- 0.2
        v_cc2 <- TRUE
      } else {
        v_l1 <- v_r1 <- TRUE
        v_ccup <- v_cc <- v_cc2 <- v_ech <- FALSE
      }
    }
  }
  
  if(vol_r < radii_ratio_th) {
    v_cc2 <- FALSE
    if(!any(mm)) {
      v_cc <- FALSE
    }
  }
  
  if(outliers > 0) {
    v_ol <- TRUE
    if(outliers > 0.05) {
      v_om <- TRUE
      if(outliers > 0.2) {
        v_oh <- TRUE
      }
    }
  }
  
  p <- ggplot() + 
    coord_equal(xlim = c(-1, 1), ylim = c(-0.6, 0.6), expand = FALSE) + 
    theme_void()
  if(v_r1) p <- p + create_rectangle(-0.6, -0.4, 1.2, 0.8, "lightgray", 
                                     "black", TRUE)
  if(v_ec2) p <- p + create_circle(x_ec2, 0, 0.4, f_ec2, c_ec2)
  if(v_ec) p <- p + create_circle(x_ec, y_ec, 0.4, NA, "black")
  if(v_eov) p <- p + create_circle(x_ec - 0.3, 0.25, 0.12, NA, "black")
  if(v_ech) p <- p + create_circle(x_ech, -0.2, 0.08, FALSE, "black")
  if(v_cc) p <- p + create_circle(x_cc, y_cc, 0.03, "black", "black")
  if(v_cc2) p <- p + create_circle(x_cc2, 0, 0.03, "black", "black")
  if(v_ccup) p <- p + create_circle(x_ccup, 0.08, 0.03, "black", "black")
  if(v_om) p <- p + create_circle(-0.15, -0.45, 0.02, "black", "black")
  if(v_oh) p <- p + create_circle(-0.3, -0.5, 0.02, "black", "black")
  if(v_om) p <- p + create_circle(0.15, -0.45, 0.02, "black", "black")
  if(v_oh) p <- p + create_circle(0.3, -0.5, 0.02, "black", "black")
  if(v_ol) p <- p + create_circle(0, -0.5, 0.02, "black", "black")

  p <- p + annotate(
    "text",
    x = 0, y = 0.45, label = k, size = 20 / .pt, vjust = "bottom"
  )
  
  p
  
}

create_rectangle <- function(x, y, width, height, fill, color, hatched) {
  
  out <- list(
    annotate(
      "rect", 
      xmin = x, ymin = y, xmax = x + width, ymax = y + height,
      fill = fill, color = color
    )
  )
  
  if(hatched) {
    l1_x <- x
    l1_xend <- l1_x + 0.1 * height 
    l1_y <- y + 0.1 * height
    l1_yend <- y
    l2_x <- x
    l2_xend <- x + height
    l2_y <- y + height 
    l2_yend <- y
    l3_x <- x + 0.9 * height
    l3_xend <- min(x + width, x + 1.9 * height)
    l3_y <- y + height
    l3_yend <- l3_y - (l3_xend - l3_x)
    out <- c(
      out,
      list(
        annotate("segment", x = l1_x, xend = l1_xend, y = l1_y, yend = l1_yend),
        annotate("segment", x = l2_x, xend = l2_xend, y = l2_y, yend = l2_yend),
        annotate("segment", x = l3_x, xend = l3_xend, y = l3_y, yend = l3_yend)
      )
    )
  }
  
  out
}

create_circle <- function(x, y, radius, fill, color) {
  annotate(
    "polygon",
    x = x + radius * cos(seq(0, 2 * pi, length.out = 100)),
    y = y + radius * sin(seq(0, 2 * pi, length.out = 100)),
    fill = fill, color = color
  )
}
