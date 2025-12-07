plot_NMDS <- function(MDS,
                      grp,
                      indic = NULL,
                      k = 2,
                      stat = T,
                      labs = NULL,
                      colors,
                      shading = "ellipse") {
  # Setup: Libraries -----------------------------------------------------------
  require(ggrepel)      # Adds geom_text_repel function
  require(ggordiplots)  # For generating hulls
  require(scales)       # provides default ggplot color palette
  require(colorspace)   # darken() color function
  require(RColorBrewer) # ggplot palettes

  # Extracts points
  MDS.pt <- data.frame(MDS$points) %>%
    rownames_to_column("Sample.ID") %>%
    left_join(data.frame(grp) %>% rownames_to_column("Sample.ID"),
              by = join_by(Sample.ID)) %>%
    rename(Groups = grp)

  # Extracts compound (species) scores
  if (!is.null(indic)) {
    MDS.Cmpd <- data.frame(MDS$species) %>%
      rownames_to_column("Cmpd") %>%
      inner_join(indic, by = join_by(Cmpd))
  }

  # Use custom labels for groups, if specified.
  # Otherwise, plot names exactly as they appear in data (default behavior)
  if (is.null(labs)) {
    labs <- levels(grp)
  }

  # Define shapes to use for points
  shapes <- c(1, 2, 3, 4, 5)

  p <- list()

  # Generate the plot ------------------------------------------------------------
  # Iterate through each set of axes
  for (i in 1:(k - 1)) {
    # Generate & extract ellipses from gg_ordiplot()
    MDS.ell <- gg_ordiplot(
      ord     = MDS,
      groups  = grp,
      scaling = 1,
      choices = c(i, i + 1),
      kind    = "sd",
      conf    = 0.95,
      show.groups = "all",
      ellipse = T,
      label   = F,
      hull    = F,
      spiders = F,
      plot    = F
    )$df_ellipse %>%
      rename(Groups = Group) %>%
      # Groups col should be a factor just like in MDS.pt so ggplot creates 1 legend
      mutate(Groups = factor(
        Groups,
        levels = levels(MDS.pt$Groups),
        ordered = T
      ))

    # Get the correct point coordinates from MDS.pt based on the axis we're plotting
    X <- paste0("MDS", 1)
    Y <- paste0("MDS", i + 1)

    # See https://forum.posit.co/t/string-as-column-name-in-ggplot/155588/2
    p[[i]] <- ggplot(data = MDS.pt, aes(x = .data[[X]], y = .data[[Y]])) +

      # Conditionally draw hulls
      list(if (shading == "hull") {
        geom_polygon(
          data = MDS.pt %>% group_by(Groups) %>% slice(chull(.data[[X]], .data[[Y]])),
          aes(fill = Groups),
          alpha = 1
        )
      }) +

      # Conditionally draws ellipses
      list(if (shading == "ellipse") {
        geom_polygon(data = MDS.ell,
                     aes(x = x, y = y, fill = Groups),
                     alpha = 1)
      }) +

      # Draw points
      geom_point(aes(color = Groups, shape = Groups),
                 size = 2,
                 stroke = 1) +

      # Label weighted av. location of samples in which indicator compounds occur
      list(if (!is.null(indic)) {
        geom_text_repel(
          data = MDS.Cmpd,
          aes(
            label = ifelse(is.na(Groups), paste0(Cmpd, "*"), Cmpd),
            color = Groups
          ),
          point.size = NA,
          bg.color = "white",
          bg.r = 0.1
        )
      }) +

      # Use custom shapes for points
      scale_shape_manual(values = head(shapes, n_distinct(pull(MDS.pt, Groups))),
                         labels = labs) +

      # Fill for ellipse & hull
      scale_fill_manual(values = colors, labels = labs) +

      # Color of points
      scale_color_manual(
        values = darken(colors, 0.5),
        labels = labs,
        guide  = "none",
        na.value = "grey30"
      ) +

      # Label axes
      xlab(paste("NMDS Axis", 1)) +
      ylab(paste("NMDS Axis", i + 1)) +

      # Clean theme w/ no axis lines.
      # Remove numbers and ticks since NMDS axes are biologically meaningless.
      theme_classic() +
      theme(
        axis.line  = element_line(),
        axis.ticks = element_blank(),
        axis.text  = element_blank()
      ) +

      # Conditionally label with K and stress
      list(if (stat == T) {
        labs(caption = paste0("k = ", MDS$ndim, ", stress = ", signif(MDS$stress, 2)))
      })
  }

  return(invisible(p))

}