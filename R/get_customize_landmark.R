#' Get the user-customized peaks and valleys location for each sample
#'
#' This function will launch a shiny app to allow user to manually set the peaks and valleys location. The function will output the landmark positions that user has set.
#' @param cell_x_adt_sample Matrix of adt counts of seleted marker, with columns of sample and batch information for each row of cells.
#' @param landmark_pos Matrix of landmark location including peaks and valleys.
#' @param bw Bandwidth for the density plot.
#' @export
#' @examples
#' \dontrun{
#' get_peak_mode(
#'   cell_x_adt = cell_x_adt,
#'   cell_x_feature = cell_x_feature,
#'   adt_marker_select = "CD3",
#'   neg_candidate_thres = asinh(6/5+1)
#' )
#' }
# require(ggplot2)
# require(ggridges)
# require(shiny)
# require(DT)

flat_table <- function(tab) {
  tmp <- lapply(1:nrow(tab), function(i) {
    data.frame(
      y = rownames(tab)[i],
      x = sapply(1:ncol(tab), function(j) {
        tab[i,j]
      }),
      name = colnames(tab)
    )
  })
  do.call(rbind, tmp)
}

create_ui <- function(landmark_pos, max_value) {
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        selectInput(inputId = "landmark_select",
                    label = "Select the sample to tune:",
                    choices = rownames(landmark_pos)),
        uiOutput(outputId = "landmark_slider_ui")
      ),
      mainPanel(
        plotOutput("plot", height='600px', width='1000px'),
        DT::dataTableOutput("table"),
        shiny::actionButton("done", "Step 2: Record User Input and Exit")
        
        # tabsetPanel(id='plot_table',
        #             tabPanel("tab_plot", 
        #                      plotOutput("plot"),
        #                      shiny::actionButton("done", "Record User Input!"),
        #                      ),
        #             tabPanel("tab_table",
        #                      DT::dataTableOutput("table")
        #                      )
        #             )
      )
    ))
  return(ui)
}


create_server <- function(landmark_pos, cell_x_adt_sample, bw = 1, brewer_palettes = "Set1", max_value) {
  fillColor = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, brewer_palettes))(length(unique(cell_x_adt_sample$batch)))

  server <- function(input, output, session) {
    vals <- shiny::reactiveValues(
      landmark_pos=landmark_pos,
      slider_values=flat_table(landmark_pos),
      fig = ggplot(cell_x_adt_sample, aes(x = adt, y = sample)) + 
        ggridges::geom_density_ridges2(aes(fill = factor(batch)), bandwidth = bw) +
        theme_bw(base_size = 20) +
        ggpubr::rotate_x_text(angle = 90) +
        ggpubr::rremove("legend") +
        xlab("ADT Counts") +
        ylab("Sample") +
        scale_fill_manual(values = fillColor)  +
        scale_x_continuous(breaks = seq(0, max_value, 0.5))

    )
    output$table <- DT::renderDataTable(vals$landmark_pos %>% DT::datatable() %>% DT::formatRound(columns = colnames(landmark_pos), digits=3))
    
    output$plot <- renderPlot(
        vals$fig +
          geom_point(data=vals$slider_values, aes(x = x, y = y), size=5)
      )
    
    output$landmark_slider_ui <- renderUI({
      if (!is.null(input$landmark_select)) {
        fluidRow(
          column(12,
                 h2(input$landmark_select),
                 lapply(1:ncol(landmark_pos), function(j) {
                   i = which(rownames(landmark_pos) == input$landmark_select)
                   val = vals$landmark_pos[i,j]
                   na_val = is.na(val)
                   val = ifelse(na_val, 0, val)
                   fluidRow(
                     checkboxInput(inputId = paste0("na_checkbox",j),
                                   label = "Set as NA",
                                   value = na_val),
                     sliderInput(
                       inputId = paste0(colnames(landmark_pos)[j]),
                       label = colnames(landmark_pos)[j],
                       min = 0,
                       max = max_value,
                       value = val,
                       step=0.0001
                     )
                   )
                 }),
                 shiny::actionButton("applied", "Step 1: Set Customized Location")
          ))
      }
    })
    
    observeEvent(input$applied, {
      i = which(rownames(landmark_pos) == input$landmark_select)
      for (j in 1:ncol(landmark_pos)) {
        if (input[[paste0("na_checkbox",j)]]) {
          vals$landmark_pos[i, j] <- NA
        } else {
          vals$landmark_pos[i, j] <- input[[paste0(colnames(landmark_pos)[j])]]
        }
      }
      output$table <- DT::renderDataTable(vals$landmark_pos %>% DT::datatable() %>% DT::formatRound(columns = colnames(landmark_pos), digits=3))
      
      slider_values = flat_table(vals$landmark_pos)
      output$plot <- renderPlot(
        vals$fig +
          geom_point(data=slider_values, aes(x = x, y = y), size=5)
      )
    })
    
    shiny::observeEvent(input$done, {
      shiny::stopApp(vals$landmark_pos)
    })
  }
  
  return(server)
}


get_customize_landmark = function(cell_x_adt_sample, landmark_pos, bw, brewer_palettes = "Set1"){

  max_value = ceiling(max(cell_x_adt_sample[, 1], na.rm = TRUE)) + 2

  ## Create an example UI
  ui <- create_ui(landmark_pos, max_value)
  server <- create_server(landmark_pos, cell_x_adt_sample, bw, brewer_palettes, max_value)

  ## Return user input, vals can change
  res <- shiny::runApp(shinyApp(ui, server))

  return(res)

}


