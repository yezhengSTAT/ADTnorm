#' Prompt Shiny browser to manually customize peaks and valleys locations.
#'
#' This function will launch a shiny app allowing the user to set the location of peaks and valleys manually. The function will output the landmark positions that the user has set.
#' @param cell_x_adt_sample Matrix of ADT counts of the selected marker, with columns of sample and batch information for each row of cells.
#' @param landmark_pos Matrix of landmark location including peaks and valleys.
#' @param bw Bandwidth for the density plot.
#' @param adt_marker_select_name The ADT marker needed to be manually processed to set the landmarks.
#' @param brewer_palettes Set the color scheme of the color brewer. The default is "Set1".
#' @examples
#' \dontrun{
#' get_customize_landmark(
#' cell_x_adt_sample, 
#' landmark_pos, 
#' bw, 
#' adt_marker_select_name, 
#' brewer_palettes = "Set1"
#' )}
#' @export
#' @import ggplot2 ggridges shiny
get_customize_landmark = function(cell_x_adt_sample, landmark_pos, bw, adt_marker_select_name, brewer_palettes = "Set1"){

  max_value = ceiling(max(cell_x_adt_sample[, 1], na.rm = TRUE)) + 2
  # cell_x_adt_sample_filter = cell_x_adt_sample %>% dplyr::filter(!is.na(adt))
  sample_num = levels(cell_x_adt_sample$sample) %>% length
  plot_height = paste0(sample_num * 50, "px")
  cell_x_adt_sample$sample_rev = factor(cell_x_adt_sample$sample, levels = rev(levels(cell_x_adt_sample$sample)))
  ## Create an example UI
  ui <- create_ui(landmark_pos, max_value, bw, plot_height)
  server <- create_server(landmark_pos, cell_x_adt_sample, bw, adt_marker_select_name, brewer_palettes, max_value)

  ## Return user input, vals can change
  res <- shiny::runApp(shinyApp(ui, server))
  return(res)

}

flat_table <- function(tab) {
  tmp <- lapply(1:nrow(tab), function(i) {
    data.frame(
      y = rownames(tab)[i],
      x = sapply(1:ncol(tab), function(j) {
        tab[i,j]
      }),
      name = colnames(tab),
      type = sub("\\d+$", "", colnames(tab))
    )
  })
  do.call(rbind, tmp)
}

create_ui <- function(landmark_pos, max_value, bw, plot_height) {
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        selectInput(inputId = "landmark_select",
                    label = "Select the sample to tune:",
                    choices = rownames(landmark_pos)),
        uiOutput(outputId = "landmark_slider_ui")
      ),
      mainPanel(
        sliderInput(inputId = "bandwidth_slider",
                    label = "Bandwidth for Density Visualization Below:",
                    min = 0,
                    max = 3,
                    value = bw,
                    step = 0.05, width = '1000px'),
        plotOutput("plot", height = plot_height, width = '1000px'),
        DT::dataTableOutput("table"),
        shiny::actionButton("done", "Step 2: Record User Input and Exit")
      )
    ))
  return(ui)
}


create_server <- function(landmark_pos, cell_x_adt_sample, bw = 0.1, adt_marker_select_name, brewer_palettes = "Set1", max_value) {
  
  fillColor <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, brewer_palettes))(length(unique(cell_x_adt_sample$batch)))

  
  server <- function(input, output, session) {
    vals <- shiny::reactiveValues(
      landmark_pos = landmark_pos,
      slider_values = flat_table(landmark_pos),
      fig = NULL
    )
    
    output$table <- DT::renderDataTable(vals$landmark_pos %>% DT::datatable() %>% DT::formatRound(columns = colnames(landmark_pos), digits = 3))

      output$plot <- renderPlot({
        bw_value <- input$bandwidth_slider
        # no_mis_idx = !is.na(cell_x_adt_sample$adt)
        
        vals$fig <- ggplot(cell_x_adt_sample, aes_string(x = "adt", y = "sample_rev")) + 
          ggridges::geom_density_ridges2(aes(fill = factor(batch)), bandwidth = bw_value) +
          theme_bw(base_size = 20) +
          ggpubr::rotate_x_text(angle = 90) +
          ggpubr::rremove("legend") +
          xlab(paste0(adt_marker_select_name, " ADT Counts")) +
          ylab("Sample") +
          scale_fill_manual(values = fillColor) +
          scale_x_continuous(breaks = seq(0, max_value, 0.5)) +
          geom_point(data = vals$slider_values, aes_string(x = "x", y = "y", shape = "type"), size = 5)
          

        vals$fig
      })

    
    output$landmark_slider_ui <- renderUI({
      if (!is.null(input$landmark_select)) {
        fluidRow(
          column(12,
                 h2(input$landmark_select),
                 lapply(1:ncol(landmark_pos), function(j) {
                   i <- which(rownames(landmark_pos) == input$landmark_select)
                   val <- vals$landmark_pos[i, j]
                   na_val <- is.na(val)
                   val <- ifelse(na_val, 0, val)
                   fluidRow(
                     checkboxInput(
                       inputId = paste0("na_checkbox", j),
                       label = "Set as NA",
                       value = na_val
                     ),
                     sliderInput(
                       inputId = paste0(colnames(landmark_pos)[j]),
                       label = colnames(landmark_pos)[j],
                       min = 0,
                       max = max_value,
                       value = val,
                       step = 0.0001
                     )
                   )
                 })),
          shiny::actionButton("applied", "Step 1: Set Customized Location")
        )
      }
    })
    
    observeEvent(input$applied, {
      i <- which(rownames(landmark_pos) == input$landmark_select)
      for (j in 1:ncol(landmark_pos)) {
        if (input[[paste0("na_checkbox", j)]]) {
          vals$landmark_pos[i, j] <- NA
        } else {
          vals$landmark_pos[i, j] <- input[[paste0(colnames(landmark_pos)[j])]]
        }
      }
      output$table <- DT::renderDataTable(vals$landmark_pos %>% DT::datatable() %>% DT::formatRound(columns = colnames(landmark_pos), digits = 3))
      
      slider_values <- flat_table(vals$landmark_pos)
      
        output$plot <- renderPlot({
          bw_value <- input$bandwidth_slider
          
          vals$fig <- ggplot(cell_x_adt_sample, aes_string(x = "adt", y = "sample_rev")) + 
            ggridges::geom_density_ridges2(aes(fill = factor(batch)), bandwidth = bw_value) +
            theme_bw(base_size = 20) +
            ggpubr::rotate_x_text(angle = 90) +
            ggpubr::rremove("legend") +
            xlab(paste0(adt_marker_select_name, " ADT Counts")) +
            ylab("Sample") +
            scale_fill_manual(values = fillColor) +
            scale_x_continuous(breaks = seq(0, max_value, 0.5)) + 
            geom_point(data = slider_values, aes_string(x = "x", y = "y", shape = "type"), size = 5) 
          vals$fig
          })
    })
    
    shiny::observeEvent(input$done, {
      shiny::stopApp(vals$landmark_pos)
    })
  }
  
  return(server)

}





