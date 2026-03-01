options(future.globals.maxSize = 100e9)
library(shiny)
library(Seurat)
library(ggplot2)
library(pins)
library(dplyr)
library(RColorBrewer)
library(bslib)
options(shiny.maxRequestSize = 1000 * 1024^2)

# /Users/sameet/Data/other-projects/Samir/analysis-jan092026/data_board/Bone.16.Met_subtypev2_stromal_combined_annotation_2026-02-10

source("utils.R")

# Define the UI for your app
ui <- fluidPage(
  titlePanel("Xenium Data Explorer"),

  sidebarLayout(
    sidebarPanel(
      textInput(
        "select_board",
        "Choose a Pin Board",
        placeholder = "Give a path to folder"
      ),
      actionButton("get_pins", "Get Pins"),
      selectInput(
        "select_pin",
        "Choose a Pin",
        choices = list(placeholder = "Type a pin name"),
        multiple = FALSE,
      ),
      actionButton("load_data", "Load Data"),
      br(),
      selectInput(
        "use_column",
        "Use Parameter",
        choices = c(
          "coarse_bone",
          "fine_bone",
          "coarse_prostate",
          "fine_prostate"
        ),
        multiple = FALSE,
        selected = "coarse_bone"
      ),
      actionButton("update_plot", "Update Plot"),

      selectInput(
        "select_gene",
        "Select Gene",
        choices = NULL,
        multiple = FALSE,
        selected = NULL
      ),
      actionButton("update_gene", "Update Gene Plot"),

      br(),
      br(),

      wellPanel(
        tags$h4("Zoom Bounding Box Coordinates"),
        fluidRow(
          column(4, numericInput("x_min", "X Min", 0, max = 6000)),
          column(4, numericInput("x_max", "X Max", 0, max = 6000)),
        ),
        fluidRow(
          column(4, numericInput("y_min", "Y Min", 0, max = 6000)),
          column(4, numericInput("y_max", "Y Max", 0, max = 6000)),
        )
      ),
      actionButton("zoom_plot", "Zoom to Bounding Box"),
      #   conditionalPanel(
      #     condition = !(is.null(sce)),
      #     selectInput(
      #       "select_column",
      #       choices = (sce@meta.data |> colnames())[grep(
      #         "bone|prostate",
      #         (sce@meta.data |> colnames())
      #       )]
      #     )
      #   )
    ),
    mainPanel(
      layout_column_wrap(
        uiOutput("sample_info"),
        uiOutput("brush_info")
      ),
      card(
        full_screen = TRUE,
        card_header(
          class = "d-flex justify-content-between align-items-center",
          "Main Image Plot",
          div(
            style = "display: flex; gap: 5px;",
            downloadButton("download_main_png", "PNG", class = "btn-sm"),
            downloadButton("download_main_pdf", "PDF", class = "btn-sm")
          )
        ),
        plotOutput(
          outputId = "plot_image"
        )
      ),
      # conditionalPanel(
      #   condition = "output.plot_image",
      #   card(
      #     card_header("Zoom Coordinates"),
      #     layout_column_wrap(
      #       width = 1 / 4,
      #       numericInput("x1", "x1", value = 0),
      #       numericInput("x2", "x2", value = 10000),
      #       numericInput("y1", "y1", value = 0),
      #       numericInput("y2", "y2", value = 10000)
      #     )
      #   )
      # ),
      card(
        full_screen = TRUE,
        card_header(
          class = "d-flex justify-content-between align-items-center",
          "Zoomed Image Plot",
          div(
            style = "display: flex; gap: 5px;",
            downloadButton("download_zoom_png", "PNG", class = "btn-sm"),
            downloadButton("download_zoom_pdf", "PDF", class = "btn-sm")
          )
        ),
        plotOutput(outputId = "plot_zoom")
      ),
      card(
        card_header(
          class = "d-flex justify-content-between align-items-center",
          "Gene Expression Plot",
          plotOutput(
            outputId = "plot_expression"
          ),
          plotOutput(
            outputId = "plot_expression_zoomed"
          )
        ),
        # conditionalPanel(
        #   #   condition = "output.plot_image",
        #   #   card(
        #   #     full_screen = TRUE,
        #   #     card_header(
        #   #       class = "d-flex justify-content-between align-items-center",
        #   #       "Gene Expression Plot",
        #   #       div(
        #   #         style = "display: flex; gap: 5px;",
        #   #         selectInput("select_gene", "Select Gene", choices = NULL),
        #   #         downloadButton("download_gene_png", "PNG", class = "btn-sm"),
        #   #         downloadButton("download_gene_pdf", "PDF", class = "btn-sm")
        #   #       )
        # ),
        # plotOutput(outputId = "plot_gene")
      )
    )
    # )
    # ),
  )
)

server <- function(input, output, session) {
  # Use reactiveValues to store data that changes
  data_store <- reactiveValues(
    sce = NULL,
    use_board = NULL
  )

  observeEvent(
    input$get_pins,
    {
      data_store$use_board <- make_board(get_dir_name(input$select_board))
      pin_lists <- data_store$use_board |> pin_list()
      pin_lists <- pin_lists[grepl("2026", pin_lists)]
      updateSelectInput(session, "select_pin", choices = pin_lists)
    }
  )

  observeEvent(
    input$load_data,
    {
      req(data_store$use_board)
      data_store$sce <- data_store$use_board |> pin_read(input$select_pin)
      updateSelectInput(
        session,
        "select_gene",
        choices = rownames(data_store$sce),
        selected = sample(rownames(data_store$sce), 1L)
      )
    }
  )

  output$sample_info <- renderUI({
    req(data_store$sce)
    value_box(
      title = "Sample Description",
      value = paste0(dim(data_store$sce), collapse = ", "),
      theme = "bg-gradient-blue-purple"
    )
  })

  output$brush_info <- renderUI({
    req(input$x_min, input$x_max, input$y_min, input$y_max)
    value_box(
      title = "Zoom Coordinates",
      value = sprintf(
        "X: (%.1f, %.1f) | Y: (%.1f, %.1f)",
        input$x_min,
        input$x_max,
        input$y_min,
        input$y_max
      ),
      theme = "bg-gradient-orange-red"
    )
  })

  # Reactive for the main plot (ggplot object)
  main_plot_obj <- eventReactive(input$update_plot, {
    req(data_store$sce)
    use_cols <- make_colors(data_store$sce, input$use_column)

    ImageDimPlot(
      data_store$sce,
      boundaries = "segmentations",
      group.by = input$use_column,
      dark.background = TRUE,
      border.color = NA,
      border.size = 0.1,
      axes = TRUE
    ) +
      scale_fill_manual(
        values = use_cols,
        breaks = names(use_cols) |> sort()
      )
  })
  output$plot_image <- renderPlot({
    req(main_plot_obj())
    main_plot_obj()
  })

  main_plot_expression <- eventReactive(input$update_gene, {
    req(data_store$sce)

    ImageFeaturePlot(
      data_store$sce,
      features = input$select_gene,
      boundaries = "segmentations",
      border.color = NA,
      border.size = 0.1,
      axes = TRUE,
      scale = "feature"
    )
  })
  output$plot_expression <- renderPlot({
    req(main_plot_expression())
    main_plot_expression()
  })

  # Reactive for the zoomed plot (ggplot object)
  zoom_plot_obj <- eventReactive(input$zoom_plot, {
    req(main_plot_obj(), input$x_min, input$y_min, input$x_max, input$y_max)
    # req(input$x_min, input$x_max, input$y_min, input$y_max)
    output$plot_image <- renderPlot({
      req(main_plot_obj())
      main_plot_obj() +
        annotate(
          "rect",
          xmin = input$x_min,
          xmax = input$x_max,
          ymin = input$y_min,
          ymax = input$y_max,
          alpha = 0.2,
          fill = "red"
        )
    })

    p <- main_plot_obj()[[1]]
    p +
      # coord_cartesian(
      xlim(input$x_min, input$x_max) +
      ylim(input$y_min, input$y_max) +
      theme_void() +
      theme(plot.background = element_rect(fill = "white"))
    # expand = FALSE
    # )
  })
  output$plot_zoom <- renderPlot({
    req(zoom_plot_obj())
    zoom_plot_obj()
  })

  zoom_gene_obj <- eventReactive(input$zoom_plot, {
    req(
      main_plot_expression(),
      input$x_min,
      input$x_max,
      input$y_min,
      input$y_max
    )
    p <- main_plot_expression()[[1]]
    p +
      xlim(input$x_min, input$x_max) +
      ylim(input$y_min, input$y_max) +
      theme_void() +
      theme(plot.background = element_rect(fill = "white"))
  })
  output$plot_expression_zoomed <- renderPlot({
    req(zoom_gene_obj())
    zoom_gene_obj()
  })

  # Reactive for gene expression plot
  gene_plot_obj <- reactive({
    req(data_store$sce, input$select_gene)
    ImageDimPlot(
      data_store$sce,
      boundaries = "segmentations",
      features = input$select_gene,
      dark.background = TRUE,
      border.color = NA,
      border.size = 0.1
    )
  })

  output$plot_gene <- renderPlot({
    req(gene_plot_obj())
    gene_plot_obj()
  })

  # Download Handlers for Main Plot
  output$download_main_png <- downloadHandler(
    filename = function() {
      paste0("main_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      ggsave(
        file,
        plot = main_plot_obj(),
        device = "png",
        width = 12,
        height = 10,
        units = "in"
      )
    }
  )

  output$download_main_pdf <- downloadHandler(
    filename = function() {
      paste0("main_plot_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      ggsave(
        file,
        plot = main_plot_obj(),
        device = "pdf",
        width = 12,
        height = 10,
        units = "in"
      )
    }
  )

  # Download Handlers for Zoomed Plot
  output$download_zoom_png <- downloadHandler(
    filename = function() {
      paste0("zoom_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      ggsave(
        file,
        plot = zoom_plot_obj(),
        device = "png",
        width = 12,
        height = 10,
        units = "in"
      )
    }
  )

  output$download_zoom_pdf <- downloadHandler(
    filename = function() {
      paste0("zoom_plot_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      ggsave(
        file,
        plot = zoom_plot_obj(),
        device = "pdf",
        width = 12,
        height = 10,
        units = "in"
      )
    }
  )

  # Download Handlers for Gene Plot
  output$download_gene_png <- downloadHandler(
    filename = function() {
      paste0("gene_plot_", input$select_gene, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      ggsave(
        file,
        plot = gene_plot_obj(),
        device = "png",
        width = 12,
        height = 10,
        units = "in"
      )
    }
  )

  output$download_gene_pdf <- downloadHandler(
    filename = function() {
      paste0("gene_plot_", input$select_gene, "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      ggsave(
        file,
        plot = gene_plot_obj(),
        device = "pdf",
        width = 12,
        height = 10,
        units = "in"
      )
    }
  )
}

shinyApp(ui, server)
