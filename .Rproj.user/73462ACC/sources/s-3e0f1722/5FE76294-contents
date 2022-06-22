plotAnnotatedStatic <- function(annotated, max_legend = 100){
  annotated <- mutate(annotated,
                      homologue_id = if_else(is.na(homologue_id), 0L, homologue_id),
                      homologue_id = as.factor(homologue_id))
  if (length(unique(annotated$homologue_id)) > max_legend) {
    legend_setting = "none"
  } else {
    legend_setting = "right"
  }
  colourCount = length(unique(annotated$homologue_id))
  getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
  ncolor = length(unique(annotated$homologue_id))
  annotated <- arrange(annotated, desc(homologue_id))
  g <- ggplot(annotated, aes(group = homologue_id)) +
    geom_point(aes(x = rt, y = mz, color = homologue_id, size = homologue_id,
                   shape = homologue_id, alpha = homologue_id)) +
    geom_line(data = filter(annotated, homologue_id != 0),
              aes(x = rt, y = mz, group = homologue_id), alpha = 0.4, size = 0.1) +
    ggtitle("Annotated Peak Table") +
    scale_colour_manual(values = c("black", getPalette(colourCount-1))) +
    scale_alpha_manual(values =  c(0.2, rep(0.85, times = ncolor))) +
    scale_size_manual(values = c(0.1, rep(1.5, times = ncolor))) +
    scale_shape_manual(values = c(1, rep(19, times = ncolor))) +
    xlab("Retention Time (s)") +
    ylab("Mass to Charge Ratio") +
    #coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) +
    theme(legend.position=legend_setting, text = element_text(family="mono"))
  return(g) # use with # ggsave("test.pdf", device = "pdf", width = 20, height = 12, units = "cm", dpi = 300)
}
# Usage Example:
# annotated <- readRDS("example-runs/flasch-msdial-results.rds")
# g <- plotAnnotatedStatic(annotated)
# g
# ggsave(plot = g, "test.pdf", device = "pdf", width = 20, height = 12, units = "cm", dpi = 300)


plotNontarget <- function(annotated_peak_table) {
  #require(Cairo)   # For nicer ggplot2 output when deployed on Linux
  annotated_peak_table <- mutate(annotated_peak_table,
                                 homologue_id = if_else(is.na(homologue_id), 0L, homologue_id),
                                 homologue_id = as.factor(homologue_id),
                                 n_assignments = as.factor(n_assignments))
  if (length(unique(annotated_peak_table$homologue_id)) > 100) {
    legend_setting = "none"
  } else {
    legend_setting = "right"
  }
  colourCount = length(unique(annotated_peak_table$homologue_id))
  getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
  ncolor = length(unique(annotated_peak_table$homologue_id))
  nassignments <- length(unique(annotated_peak_table$n_assignments))
  annotated_peak_table <- arrange(annotated_peak_table, desc(homologue_id))
  ui <- fluidPage(
    plotOutput("plot", width = '100%', height = 600,
               dblclick = dblclickOpts("plot_dblclick"),
               brush = brushOpts(
                 id = "plot_brush",
                 resetOnNew = TRUE)),
    tableOutput("annotated_peak_table")
  )

  server <- function(input, output) {
    ranges <- reactiveValues(x = NULL, y = NULL) # Null defaults
    output$plot <- renderPlot({
      ggplot(annotated_peak_table, aes(group = homologue_id)) +
        geom_point(aes(x = rt, y = mz, color = n_assignments, size = homologue_id,
                       shape = n_assignments, alpha = n_assignments)) +
        geom_line(data = filter(annotated_peak_table, homologue_id != 0),
                  aes(x = rt, y = mz, group = homologue_id), alpha = 0.5) +
        ggtitle("Peak Table") +
        scale_colour_manual(values = c("grey20", rep("red", times = nassignments))) +
        scale_size_manual(values = c(0.1, rep(4, times = ncolor))) +
        scale_shape_manual(values = c(19, rep(17, times = nassignments) )) +
        scale_alpha_manual(values = c(0.1, rep(1, times = nassignments))) +
        coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) +
        theme(legend.position=legend_setting)
    })
    # When a double-click happens, check if there's a brush on the plot.
    # If so, zoom to the brush bounds; if not, reset the zoom.
    observeEvent(input$plot_dblclick, {
      brush <- input$plot_brush
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)
      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    })
    output$annotated_peak_table <- renderTable({
      brushedPoints(annotated_peak_table, input$plot_brush)
    })
  }
  shinyApp(ui, server)
}

createPseudoPeakTableFromNonTarget <- function(nontarget_out){
  peak_table <- as_tibble(nontarget_out$`Peaks in homologue series`)
  names(peak_table) <- stringr::str_replace_all(names(peak_table), " |/", "_")
  all(peak_table$HS_IDs == peak_table$HS_cluster)
  # unnesting creates duplicate rows immediately! So for each peak with n series assignments, there will now be n rows.
  peak_table <- peak_table %>%
    mutate(., series_ids = str_split(HS_cluster, pattern = "/")) %>%
    unnest(., series_ids)



  peak_table <- peak_table %>%
    rename(., mz = mz, rt = RT, homologue_id = series_ids, nontarget_peak_id = peak_ID) %>%
    mutate(., homologue_id = as.integer(homologue_id), nontarget_peak_id = as.integer(nontarget_peak_id)) %>%
    mutate(., homologue_id = if_else(homologue_id == 0, as.integer(NA), homologue_id)) %>%
    select(., mz, rt, homologue_id, nontarget_peak_id, intensity)

  peak_table <- peak_table %>%
    arrange(., mz, rt) %>%
    group_by(., homologue_id) %>%
    mutate(., within_series_id = row_number()) %>%
    mutate(., within_series_id = if_else(homologue_id == 0L, as.integer(NA), within_series_id)) %>%
    ungroup(.)

  peak_table <- peak_table %>%
    group_by(., nontarget_peak_id) %>%
    add_tally(name = "n_assignments") %>%
    ungroup(.)
  return (peak_table)
}


display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ..., disable.logging = TRUE)
  grid.draw(venn_object)
}

createSeriesStrings <- function(out){
  summary <- out %>%
    na.omit(.) %>%
    group_by(homologue_id) %>%
    arrange(., peak_id) %>%
    summarize(., series_string = paste(peak_id, collapse = "-->"))
  strings <- pull(summary, series_string)
  return(strings)
}

createAnnotedPeakIdList<- function(out){
  ids <- out %>%
    na.omit(.) %>%
    pull(peak_id) %>%
    unique(.)
  return(ids)
}
