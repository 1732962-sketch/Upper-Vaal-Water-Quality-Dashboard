# UPPER VAAL WATER QUALITY DASHBOARD
# Run from RStudio: shiny::runApp("R/App_New.R")

required <- c("shiny", "bslib", "sf", "dplyr", "leaflet", "ggplot2", "DT", "terra", "gstat")
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install packages: ", paste(missing_packages, collapse = ", "))

library(shiny)
library(bslib)
library(sf)
library(dplyr)
library(leaflet)
library(ggplot2)
library(DT)
library(terra)
library(gstat)


# Project directory for local and online use
project_dir <- getwd() 
project_file <- function(...) file.path(project_dir, ...)
load_csv <- function(...) {
  path <- project_file(...)
  if (!file.exists(path)) stop("Missing file: ", path)
  read.csv(path, stringsAsFactors = FALSE)
}
load_rds <- function(...) {
  path <- project_file(...)
  if (!file.exists(path)) stop("Missing file: ", path)
  readRDS(path)
}
load_shape <- function(filename) {
  path <- project_file("data", "raw", "Clip Features", filename)
  if (!file.exists(path)) {
    warning("Optional shapefile not found: ", path)
    return(NULL)
  }
  tryCatch(st_transform(st_read(path, quiet = TRUE), 4326),
           error = function(e) {
             warning("Could not load ", filename, ": ", conditionMessage(e))
             NULL
           })
}

water_data <- load_csv("data", "raw", "upper_vaal_annual_medians.csv")
global_results <- load_csv("data", "processed", "global_moran_results.csv")
lisa_results <- st_transform(load_rds("data", "processed", "lisa_all_results.rds"), 4326)
kriging_results <- load_rds("data", "processed", "kriging_all_results.rds")
kriging_summary <- load_csv("data", "processed", "kriging_summary.csv")

shapes <- list(
  "Upper Vaal Boundary" = load_shape("WMA_Clip.shp"),
  "Primary Catchments" = load_shape("Primary_catch_clip.shp"),
  "Secondary Catchments" = load_shape("Secondary_catch_clip.shp"),
  "Quaternary Catchments" = load_shape("Quaternary_catch.shp"),
  "Rivers" = load_shape("Rivers_clip.shp"),
  "Dams" = load_shape("Dams_clip.shp"),
  "Land Use" = load_shape("Landuse_NGI clip.shp")
)

parameter_choices <- c("pH" = "pH_median", "Electrical Conductivity (EC)" = "EC_median", "Sulphate (SO4)" = "SO4_median")
parameter_label <- function(p) names(parameter_choices)[match(p, parameter_choices)]
methods <- c("Global Moran's I", "Local Moran's I (LISA)", "Ordinary Kriging")
# ==========================================
# NEW LAND USE / LAND COVER RASTER
# ==========================================

# Load the optimised land-use raster
lulc_raster <- terra::rast(
  "data/processed/lulc_web.tif"
)

# Convert the categorical raster to its numeric class codes
lulc_raster <- as.numeric(lulc_raster)

# Land-use class codes
lulc_values <- c(
  20, 21, 32, 38, 39, 40, 41, 66, 69, 71, 72
)

# Land-use colours
lulc_colors <- c(
  "#2166AC",  # 20 - Artificial sewage ponds
  "#D73027",  # 21 - Artificial flooded mine pits
  "#F46D43",  # 32 - Permanent orchards
  "#FFD92F",  # 38 - Pivot irrigated crops
  "#1A9850",  # 39 - Non-pivot irrigated crops
  "#A6D96A",  # 40 - Rain-fed crops
  "#66C2A5",  # 41 - Subsistence crops
  "#3288BD",  # 66 - Industrial
  "#984EA3",  # 69 - Mining extraction pits
  "#E7298A",  # 71 - Mine tailings
  "#E41A1C"   # 72 - Landfills
)

# Create a colour palette for categorical values
lulc_palette <- leaflet::colorFactor(
  palette = lulc_colors,
  domain = lulc_values,
  levels = lulc_values,
  na.color = "transparent"
)

# Map overlay helpers: avoid errors if an optional layer is missing
add_base_layers <- function(map, layer_names) {
  if ("Land Use" %in% layer_names) {
    
    # Display the new Land Use raster
    map <- leaflet::addRasterImage(
      map,
      lulc_raster,
      colors = lulc_palette,
      opacity = 0.75,
      project = TRUE,
      maxBytes = 20 * 1024 * 1024,
      group = "Land Use"
    )
    
    # Add the Land Use classification legend
    map <- leaflet::addLegend(
      map,
      position = "bottomleft",
      colors = lulc_colors,
      labels = c(
        "Artificial sewage ponds",
        "Artificial flooded mine pits",
        "Cultivated commercial permanent orchards",
        "Commercial annual crops - pivot irrigated",
        "Commercial annual crops - non-pivot irrigated",
        "Commercial annual crops - rain-fed/dryland",
        "Subsistence/small-scale annual crops",
        "Industrial",
        "Mines: extraction pits, quarries",
        "Mine: tailings and resource dumps",
        "Landfills"
      ),
      title = "Land Use Classification",
      opacity = 1,
      group = "Land Use"
    )
    
  }
  for (nm in c("Primary Catchments", "Secondary Catchments", "Quaternary Catchments")) {
    if (nm %in% layer_names && !is.null(shapes[[nm]])) {
      col <- switch(nm, "Primary Catchments" = "#467D60",
                    "Secondary Catchments" = "#8D6C3C", "Quaternary Catchments" = "#999999")
      map <- addPolygons(map, data = shapes[[nm]], color = col, weight = 1,
                         fill = FALSE, group = nm)
    }
  }
  if ("Dams" %in% layer_names && !is.null(shapes[["Dams"]])) {
    dams <- shapes[["Dams"]]
    geom <- unique(as.character(st_geometry_type(dams)))
    if (all(grepl("POLYGON", geom))) {
      map <- addPolygons(map, data = dams, color = "#287BA8", fillColor = "#67B7D7",
                         fillOpacity = 0.75, weight = 1, group = "Dams")
    } else if (all(grepl("POINT", geom))) {
      map <- addCircleMarkers(map, data = dams, radius = 4, color = "#287BA8",
                              group = "Dams")
    } else if (all(grepl("LINESTRING", geom))) {
      map <- addPolylines(map, data = dams, color = "#287BA8", weight = 2,
                          group = "Dams")
    }
  }
  if ("Rivers" %in% layer_names && !is.null(shapes[["Rivers"]])) {
    map <- addPolylines(map, data = shapes[["Rivers"]], color = "#2A8DC4",
                        weight = 1.2, opacity = 0.8, group = "Rivers")
  }
  if ("Upper Vaal Boundary" %in% layer_names && !is.null(shapes[["Upper Vaal Boundary"]])) {
    map <- addPolygons(map, data = shapes[["Upper Vaal Boundary"]],
                       color = "#205844", weight = 2.5, fill = FALSE,
                       group = "Upper Vaal Boundary")
  }
  map
}
# Shared layer selection; Shiny filters determine which layers are initially shown
available_layers <- names(shapes)[!vapply(shapes, is.null, logical(1))]
default_layers <- intersect(c("Upper Vaal Boundary", "Rivers", "Dams"), available_layers)

ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  tags$head(tags$style(HTML("\n    body {background: #f4f7f8;}\n    .dashboard-title {background:#416c70; color:white; padding:18px 24px; margin:-15px -15px 16px;}\n    .side-panel {background:#edf2f2; padding:16px; border-radius:8px; margin-bottom:15px;}\n    .metric {background:white; border:1px solid #dce5e7; border-radius:8px; padding:16px; margin-bottom:14px;}\n    .metric strong {font-size:26px; display:block; color:#243e42;}\n    .tab-content {background:white; padding:14px; border:1px solid #e4e9ec; border-top:0;}
    .intro-copy p {font-size:14px; margin-bottom:9px;}
    .intro-copy li {font-size:13px; margin-bottom:5px;}
    .intro-copy h4 {margin-top:9px;}\n  "))),
  div(class = "dashboard-title",
      h2("Upper Vaal Water Quality Dashboard"),
      p(style = "margin:4px 0 0; font-size:15px; color:#f1f8f8;",
        strong("Project team: "),
        "Unarine Tshilande | Cokisa Gabavana | Kanukani Khorombi | Neo Hlalele")),
  sidebarLayout(
    sidebarPanel(width = 3, class = "side-panel",
                 h4("Analysis Controls"),
                 selectInput("year", "Select Year", choices = 2016:2020, selected = 2019),
                 selectInput("parameter", "Water Quality Parameter", choices = parameter_choices),
                 selectInput("method", "Statistical Method", choices = methods, selected = "Ordinary Kriging"),
                 helpText("Choose a year, parameter and one of the three spatial statistical methods.")
    ),
    mainPanel(width = 9,
              fluidRow(
                column(4, div(class = "metric", "Study Period", tags$strong("2016–2020"))),
                column(4, div(class = "metric", "Water Quality Parameters", tags$strong("3"))),
                column(4, div(class = "metric", "Monitoring Stations", tags$strong(length(unique(water_data$dws_dataset_id)))))
              ),
              tabsetPanel(id = "main_tab",
                          tabPanel("Introduction", br(),
                                   fluidRow(
                                     column(7, div(class = "intro-copy",
                                                   h3("Upper Vaal Catchment"),
                                                   p("The Upper Vaal Catchment is part of the Vaal River Basin, spanning parts of Gauteng, Mpumalanga and the Free State. It includes the Vaal Dam and supports domestic, industrial, mining, agricultural and energy-related water needs."),
                                                   p("Mining, wastewater, industry and agricultural runoff can affect pH, electrical conductivity and sulphate. Understanding where and when these parameters vary can help prioritise monitoring and further investigation."),
                                                   h4("Why Spatial Statistical Analysis?"),
                                                   
                                                   p("We use spatial statistical analysis to understand how water quality varies across the Upper Vaal catchment, identify areas with similar water-quality conditions, and estimate values where monitoring stations are absent."),
                                                   
                                                   tags$ul(
                                                     tags$li(
                                                       tags$b("Global Moran's I: "),
                                                       "Determines whether water-quality values are spatially clustered, dispersed or randomly distributed across the study area."
                                                     ),
                                                     
                                                     tags$li(
                                                       tags$b("Local Moran's I (LISA): "),
                                                       "Identifies specific locations with high or low water-quality clusters and spatial outliers."
                                                     ),
                                                     
                                                     tags$li(
                                                       tags$b("Ordinary Kriging: "),
                                                       "Predicts water-quality values at locations without monitoring stations, producing a continuous prediction map."
                                                     )
                                                   ),
                                                   p(tags$b("Study data: "), "45 stations, 2016–2020, annual station medians for pH, EC and sulphate."),
                                                   tags$small("Spatial associations and predictions do not establish pollution causes or regulatory exceedances.")
                                     )),
                                     column(5, h4("Study area and stations"), leafletOutput("intro_map", height = "310px"))
                                   )
                          ),
                          tabPanel("Overview", br(), h4(textOutput("selection_title")),
                                   tableOutput("overview_summary"),
                                   h4("Monitoring stations"), DTOutput("overview_table")),
                          tabPanel("Spatial Map", br(), h4(textOutput("map_title")),
                                   fluidRow(
                                     column(9, leafletOutput("analysis_map", height = 590)),
                                     column(3, div(class = "side-panel",
                                                   h4("Map layers"),
                                                   checkboxGroupInput("layers", NULL, choices = available_layers,
                                                                      selected = default_layers),
                                                   checkboxInput("show_stations", "Monitoring stations", value = TRUE),
                                                   tags$small("Turn GIS layers on or off. These are map overlays, not statistical methods.")
                                     ))
                                   )),
                          tabPanel("Statistical Results", br(), h4(textOutput("stat_title")),
                                   uiOutput("stat_explanation"),
                                   h4("Selected year — detailed findings"),
                                   tableOutput("stat_table"),
                                   h4("Five-year comparison (2016–2020)"),
                                   plotOutput("global_year_plot", height = 280),
                                   tableOutput("five_year_table"),
                                   uiOutput("five_year_interpretation"),
                                   conditionalPanel(
                                     condition = "input.method == 'Ordinary Kriging'",
                                     hr(), h4("Experimental variogram and fitted model"),
                                     p("Variogram diagnostic calculated from the selected year and parameter. The curve is fitted for visual explanation and may differ from the model used to produce the saved Kriging predictions."),
                                     plotOutput("variogram_plot", height = 330),
                                     textOutput("variogram_note")
                                   )),
                          tabPanel("Trend Analysis", br(),
                                   h3("Annual Trend (2016–2020)"),
                                   p("Each annual value is the median of station-level annual medians, not the median of all individual samples."),
                                   plotOutput("annual_trend", height = 300),
                                   hr(), h3("Individual Monitoring Station Trend"),
                                   selectInput("station", "Select Monitoring Station", choices = NULL),
                                   plotOutput("station_trend", height = 300),
                                   tableOutput("station_values")),
                          tabPanel("Data Explorer", br(), DTOutput("explorer_table")),
                          tabPanel("Conclusion", br(),
                                   h3("Conclusion and key findings"),
                                   p("This dashboard explores spatial and temporal variation in pH, electrical conductivity (EC) and sulphate across 45 Upper Vaal monitoring stations between 2016 and 2020."),
                                   h4("Global spatial patterns"),
                                   p("Global Moran's I evaluates whether similar water-quality values tend to occur near one another across the study area. Results vary by parameter and year; consult the Statistical Results tab for the actual values and significance tests."),
                                   h4("Local patterns"),
                                   p("Local Moran's I identifies potential local clusters and spatial outliers. These are exploratory patterns rather than proof of a pollution source."),
                                   h4("Predicted water quality"),
                                   p("Ordinary Kriging estimates water-quality values between monitoring stations. The prediction surface is a modelled estimate, not a set of direct measurements. Variograms describe spatial dependence, while cross-validation RMSE provides a measure of prediction performance."),
                                   h4("Limitations"),
                                   tags$ul(
                                     tags$li("Only 45 monitoring stations and the 2016–2020 period are included."),
                                     tags$li("Annual station medians summarise multiple samples and may hide short-term changes."),
                                     tags$li("Kriging uncertainty increases where monitoring coverage is sparse; the smoother display does not improve prediction accuracy."),
                                     tags$li("Spatial association does not establish causation or regulatory non-compliance.")
                                   ),
                                   h4("Recommendation"),
                                   p("Use the spatial results to guide further monitoring, targeted investigation and data-quality checks, alongside hydrological and catchment information."))
              )
    )
  )
)

server <- function(input, output, session) {
  observe({
    stations <- water_data |>
      distinct(dws_dataset_id, map_label) |>
      arrange(map_label)
    choices <- stats::setNames(as.character(stations$dws_dataset_id),
                               paste0(stations$map_label, " (", stations$dws_dataset_id, ")"))
    updateSelectInput(session, "station", choices = choices,
                      selected = if (length(choices)) unname(choices[1]) else character(0))
  })
  
  selected_data <- reactive({
    water_data |> filter(year == as.integer(input$year))
  })
  selected_lisa <- reactive({
    lisa_results |> filter(year == as.integer(input$year), Parameter == input$parameter)
  })
  selected_kriging <- reactive({
    key <- paste(input$year, input$parameter, sep = "_")
    req(kriging_results[[key]])
    st_transform(kriging_results[[key]], 4326)
  })
  selected_global <- reactive({
    global_results |> filter(Year == as.integer(input$year), Parameter == input$parameter)
  })
  selected_model <- reactive({
    kriging_summary |> filter(Year == as.integer(input$year), Parameter == input$parameter)
  })
  
  # Variogram diagnostic: derived from observed annual station medians.
  # This fitted curve is for explanation and is NOT substituted into saved predictions.
  output$variogram_plot <- renderPlot({
    req(input$method == "Ordinary Kriging")
    stations <- selected_data()
    stations <- stations[is.finite(stations[[input$parameter]]) &
                           is.finite(stations$longitude) & is.finite(stations$latitude), ]
    validate(need(nrow(stations) >= 8, "Insufficient valid stations for variogram."))
    pts <- sf::st_as_sf(stations, coords = c("longitude", "latitude"), crs = 4326)
    pts <- sf::st_transform(pts, 32735)
    v <- gstat::variogram(stats::as.formula(paste(input$parameter, "~ 1")), pts)
    validate(need(nrow(v) >= 3, "Insufficient distance bins for variogram."))
    fitted <- tryCatch({
      initial <- gstat::vgm(psill = max(v$gamma, na.rm = TRUE),
                            model = "Sph", range = max(v$dist, na.rm = TRUE) / 3,
                            nugget = min(v$gamma, na.rm = TRUE))
      gstat::fit.variogram(v, initial)
    }, error = function(e) NULL)
    chart <- ggplot(v, aes(x = dist / 1000, y = gamma)) +
      geom_point(color = "#247C8A", size = 2.7) +
      labs(x = "Distance (km)", y = "Semivariance",
           title = paste("Experimental variogram:", parameter_label(input$parameter), input$year)) +
      theme_minimal()
    if (!is.null(fitted) && all(is.finite(fitted$psill)) && all(fitted$psill >= 0)) {
      line <- gstat::variogramLine(fitted, maxdist = max(v$dist), n = 200)
      chart <- chart + geom_line(data = line,
                                 aes(x = dist / 1000, y = gamma),
                                 inherit.aes = FALSE, color = "#D1783B", linewidth = 1)
    }
    chart
  })
  output$variogram_note <- renderText({
    "Dots: experimental semivariance; orange line (if available): illustrative fitted spherical model. This diagnostic is not a replacement for the original saved Kriging model."
  })
  
  output$selection_title <- renderText(paste("Water quality overview:", parameter_label(input$parameter), input$year))
  output$map_title <- renderText(paste(input$method, "—", parameter_label(input$parameter), input$year))
  output$stat_title <- renderText(paste(input$method, "—", parameter_label(input$parameter), input$year))
  
  output$overview_summary <- renderTable({
    x <- selected_data()[[input$parameter]]
    data.frame(Stations = length(unique(selected_data()$dws_dataset_id)),
               Minimum = round(min(x, na.rm = TRUE), 3),
               Median = round(median(x, na.rm = TRUE), 3),
               Maximum = round(max(x, na.rm = TRUE), 3))
  })
  output$overview_table <- renderDT({
    datatable(st_drop_geometry(selected_lisa()) |>
                select(dws_dataset_id, map_label, year, all_of(input$parameter)),
              options = list(pageLength = 8))
  })
  
  # The introduction map is independent of analysis controls
  output$intro_map <- renderLeaflet({
    x <- water_data |> filter(year == 2016)
    points <- st_as_sf(x, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
    map <- leaflet() |> addTiles()
    map <- add_base_layers(map, default_layers)
    map <- addCircleMarkers(
      map, data = points, radius = 4, color = "#B43434",
      fillOpacity = 0.9, group = "Monitoring Stations",
      popup = ~paste0("<b>", map_label, "</b>")
    )
    addLegend(
      map, position = "bottomright",
      colors = "#B43434", labels = "Monitoring stations",
      title = "Study area", opacity = 0.9
    )
  })
  
  output$analysis_map <- renderLeaflet({
    req(input$method, input$year, input$parameter)
    map <- leaflet() |> addTiles()
    if (input$method == "Ordinary Kriging") {
      pred <- selected_kriging()
      vals <- pred$var1.pred
      pal <- colorNumeric("YlGnBu", domain = range(vals, na.rm = TRUE))
      # Display the SAME saved predictions as a continuous colour surface.
      # Bilinear resampling smooths only the display, not the Kriging model.
      map <- tryCatch({
        projected <- kriging_results[[paste(input$year, input$parameter, sep = "_")]]
        xy <- sf::st_coordinates(projected)
        xyz <- data.frame(x = xy[, 1], y = xy[, 2], prediction = projected$var1.pred)
        
        # Create the Kriging prediction raster
        surface <- terra::rast(
          xyz,
          type = "xyz",
          crs = sf::st_crs(projected)$wkt
        )
        
        # Smooth the display
        surface <- terra::disagg(
          surface,
          fact = 3,
          method = "bilinear"
        )
        
        # Get the Upper Vaal catchment boundary
        boundary <- shapes[["Upper Vaal Boundary"]]
        
        if (!is.null(boundary)) {
          
          # Convert the boundary to a raster-compatible format
          boundary_vect <- terra::vect(boundary)
          
          # Match the boundary to the raster coordinate system
          boundary_vect <- terra::project(
            boundary_vect,
            terra::crs(surface)
          )
          
          # Make cells outside the catchment transparent
          surface <- terra::mask(
            surface,
            boundary_vect
          )
        }
        
        # Display the smooth Kriging surface
        map <- leaflet::addRasterImage(
          map,
          surface,
          colors = pal,
          opacity = 0.78,
          project = TRUE,
          group = "Kriging Predictions"
        )
      }, error = function(e) {
        warning("Raster display unavailable; using original prediction markers: ", conditionMessage(e))
        leaflet::addCircleMarkers(map, data = pred, radius = 3, stroke = FALSE,
                                  fillColor = ~pal(var1.pred), fillOpacity = 0.75,
                                  group = "Kriging Predictions")
      })
      map <- addLegend(map, position = "bottomright", pal = pal,
                       values = vals,
                       title = paste("Kriging prediction:", parameter_label(input$parameter)),
                       opacity = 0.9)
    }
    map <- add_base_layers(map, input$layers)
    pts <- selected_lisa()
    if (input$method == "Local Moran's I (LISA)") {
      colors <- c("High-High" = "#D7191C", "Low-Low" = "#2C7BB6",
                  "High-Low" = "#FDAE61", "Low-High" = "#ABD9E9",
                  "Not significant" = "#A9A9A9")
      pal <- colorFactor(unname(colors), domain = names(colors))
      map <- addCircleMarkers(map, data = pts, radius = 7,
                              color = ~pal(Cluster), fillColor = ~pal(Cluster),
                              fillOpacity = 0.9, group = "LISA Clusters",
                              popup = ~paste0("Station: ", map_label,
                                              "<br>Cluster: ", Cluster,
                                              "<br>Local I: ", round(Local_I, 3),
                                              "<br>p-value: ", signif(P_value, 3)))
      map <- addLegend(map, position = "bottomleft", colors = unname(colors),
                       labels = names(colors), title = "LISA Clusters", opacity = 0.9)
    } else if (isTRUE(input$show_stations)) {
      pts$display_value <- pts[[input$parameter]]
      map <- addCircleMarkers(map, data = pts, radius = 5, color = "#123E52",
                              fillColor = "#42B8C4", fillOpacity = 0.95,
                              group = "Monitoring Stations",
                              popup = ~paste0("Station: ", map_label,
                                              "<br>Value: ", round(display_value, 3)))
    }
    groups <- unique(c(input$layers,
                       if (input$method == "Ordinary Kriging") "Kriging Predictions" else
                         if (input$method == "Local Moran's I (LISA)") "LISA Clusters" else NULL,
                       if (isTRUE(input$show_stations) && input$method != "Local Moran's I (LISA)") "Monitoring Stations" else NULL))
    # A second, compact legend explains only the GIS layers currently enabled.
    layer_symbols <- c(
      "Upper Vaal Boundary" = "#205844",
      "Primary Catchments" = "#467D60",
      "Secondary Catchments" = "#8D6C3C",
      "Quaternary Catchments" = "#999999",
      "Rivers" = "#2A8DC4",
      "Dams" = "#67B7D7",
      "Land Use" = "#C8D9A4"
    )
    active <- intersect(names(layer_symbols), input$layers)
    if (isTRUE(input$show_stations) && input$method != "Local Moran's I (LISA)") {
      active <- c(active, "Monitoring Stations")
    }
    if (length(active)) {
      symbol_html <- vapply(active, function(nm) {
        col <- if (nm == "Monitoring Stations") "#42B8C4" else layer_symbols[[nm]]
        kind <- if (nm == "Monitoring Stations") "border-radius:50%;" else ""
        paste0('<div style="display:flex;align-items:center;gap:7px;margin:3px 0;">',
               '<span style="display:inline-block;width:14px;height:10px;flex-shrink:0;',
               'background:', col, ';border:1px solid #65727a;', kind, '"></span>',
               htmltools::htmlEscape(nm), '</div>')
      }, character(1))
      map <- addControl(map,
                        html = paste0('<div style="background:rgba(255,255,255,.96);padding:8px 10px;',
                                      'border-radius:5px;box-shadow:0 1px 6px #999;font-size:11px;',
                                      'max-height:190px;overflow-y:auto;">',
                                      '<b>Visible GIS layers</b>', paste(symbol_html, collapse = ''), '</div>'),
                        position = "topright")
    }
    if (input$method == "Global Moran's I") {
      g <- selected_global()
      if (nrow(g)) {
        msg <- if (g$P_value[1] < 0.05) "Significant (p &lt; 0.05)" else "Not significant"
        map <- addControl(map,
                          html = paste0('<div style="background:white;padding:9px;border-radius:5px;',
                                        'box-shadow:0 1px 6px #999;font-size:12px;">',
                                        '<b>Global Moran&#39;s I</b><div>I = ',
                                        round(g$Morans_I[1], 3), '</div><div>p = ',
                                        signif(g$P_value[1], 3), '</div><div>', msg, '</div></div>'),
                          position = "bottomright")
      }
    }
    if (length(groups) > 0) {
      map <- addLayersControl(
        map, overlayGroups = groups,
        options = layersControlOptions(collapsed = TRUE)
      )
    }
    map
  })
  
  # Five-year summaries are recalculated from the saved analysis outputs.
  lisa_comparison <- reactive({
    x <- st_drop_geometry(lisa_results) |>
      filter(Parameter == input$parameter) |>
      count(year, Cluster, name = "Stations")
    x
  })
  kriging_comparison <- reactive({
    x <- kriging_summary |> filter(Parameter == input$parameter)
    x$Pred_Min <- vapply(x$Year, function(y) {
      z <- kriging_results[[paste(y, input$parameter, sep = "_")]]
      if (is.null(z)) return(NA_real_)
      min(z$var1.pred, na.rm = TRUE)
    }, numeric(1))
    x$Pred_Max <- vapply(x$Year, function(y) {
      z <- kriging_results[[paste(y, input$parameter, sep = "_")]]
      if (is.null(z)) return(NA_real_)
      max(z$var1.pred, na.rm = TRUE)
    }, numeric(1))
    x
  })
  output$stat_explanation <- renderUI({
    if (input$method == "Global Moran's I") {
      x <- selected_global()
      req(nrow(x) == 1)
      significant <- x$P_value[1] < 0.05
      tagList(
        p(sprintf("In %s, Global Moran's I for %s was %.3f (two-sided p = %.3g).",
                  input$year, parameter_label(input$parameter), x$Morans_I[1], x$P_value[1])),
        p(if (significant && x$Morans_I[1] > 0) "This indicates statistically significant positive spatial autocorrelation: similar values tend to occur near one another."
          else if (significant && x$Morans_I[1] < 0) "This indicates statistically significant negative spatial autocorrelation: neighbouring values tend to differ."
          else "This result does not provide sufficient evidence of catchment-wide spatial autocorrelation at the 5% level."),
        tags$small("The significance test uses unadjusted two-sided p-values; association does not establish a cause.")
      )
    } else if (input$method == "Local Moran's I (LISA)") {
      x <- selected_lisa()
      counts <- table(factor(x$Cluster, levels = c("High-High", "Low-Low", "High-Low", "Low-High", "Not significant")))
      tagList(p(sprintf("For %s in %s, %s of 45 stations have locally significant clusters or outlier patterns (unadjusted p < 0.05).",
                        parameter_label(input$parameter), input$year, sum(counts[1:4]))),
              p(sprintf("High–High: %s; Low–Low: %s; High–Low: %s; Low–High: %s; not significant: %s.",
                        counts[1], counts[2], counts[3], counts[4], counts[5])),
              tags$small("Local significance is exploratory; multiple testing can produce false positives."))
    } else {
      x <- selected_model()
      pred <- selected_kriging()
      req(nrow(x) == 1)
      tagList(p(sprintf("The selected %s %s Kriging model is %s, with leave-one-out RMSE %.3f.",
                        input$year, parameter_label(input$parameter), x$Model[1], x$RMSE[1])),
              p(sprintf("Predictions across the grid range from %.3f to %.3f. These are model estimates, not observed measurements.",
                        min(pred$var1.pred, na.rm = TRUE), max(pred$var1.pred, na.rm = TRUE))),
              tags$small("RMSE units follow the selected parameter. Extreme predictions require validation; uncertainty is represented by Kriging variance."))
    }
  })
  
  output$stat_table <- renderTable({
    
    if (input$method == "Global Moran's I") {
      
      x <- selected_global()
      
      x$Significant_5pct <- ifelse(
        x$P_value < 0.05, "Yes", "No"
      )
      
      return(x)
    }
    
    if (input$method == "Local Moran's I (LISA)") {
      
      x <- selected_lisa()
      
      clusters <- c(
        "High-High",
        "Low-Low",
        "High-Low",
        "Low-High",
        "Not significant"
      )
      
      counts <- table(
        factor(x$Cluster, levels = clusters)
      )
      
      return(
        data.frame(
          Cluster = clusters,
          Stations = as.integer(counts)
        )
      )
    }
    
    selected_model()
  })
  
  output$five_year_table <- renderTable({
    if (input$method == "Global Moran's I") {
      x <- global_results |> filter(Parameter == input$parameter) |> arrange(Year)
      return(data.frame(Year = x$Year, Morans_I = round(x$Morans_I, 4),
                        P_value = signif(x$P_value, 3),
                        Significant = ifelse(x$P_value < 0.05, "Yes", "No")))
    }
    if (input$method == "Local Moran's I (LISA)") {
      x <- lisa_comparison()
      tab <- xtabs(Stations ~ factor(year, levels = 2016:2020) +
                     factor(Cluster, levels = c("High-High", "Low-Low", "High-Low", "Low-High", "Not significant")), data = x)
      colnames(tab) <- c("High-High", "Low-Low", "High-Low", "Low-High", "Not significant")
      return(data.frame(Year = as.integer(rownames(tab)), tab[, c("High-High", "Low-Low", "High-Low", "Low-High", "Not significant"), drop = FALSE], check.names = FALSE))
    }
    x <- kriging_comparison()
    data.frame(Year = x$Year, Model = x$Model, RMSE = round(x$RMSE, 3),
               Predicted_Min = round(x$Pred_Min, 3), Predicted_Max = round(x$Pred_Max, 3))
  }, striped = TRUE)
  output$global_year_plot <- renderPlot({
    if (input$method == "Global Moran's I") {
      x <- global_results |> filter(Parameter == input$parameter)
      ggplot(x, aes(Year, Morans_I)) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        geom_line(color = "#317F96", linewidth = 1) +
        geom_point(aes(color = P_value < 0.05), size = 3) +
        scale_color_manual(values = c("TRUE" = "#16734C", "FALSE" = "#A76C32"),
                           labels = c("FALSE" = "Not significant", "TRUE" = "Significant"), name = "p < 0.05") +
        scale_x_continuous(breaks = 2016:2020) +
        labs(x = "Year", y = "Global Moran's I", title = "Global spatial autocorrelation over time") + theme_minimal()
    } else if (input$method == "Local Moran's I (LISA)") {
      x <- lisa_comparison()
      ggplot(x, aes(x = factor(year), y = Stations, fill = Cluster)) +
        geom_col() +
        scale_fill_manual(values = c("High-High" = "#D7191C", "Low-Low" = "#2C7BB6", "High-Low" = "#FDAE61", "Low-High" = "#ABD9E9", "Not significant" = "#A9A9A9")) +
        labs(x = "Year", y = "Monitoring stations", title = "Local cluster classifications by year") + theme_minimal()
    } else {
      x <- kriging_comparison()
      ggplot(x, aes(x = Year, y = RMSE)) + geom_line(color = "#317F96", linewidth = 1) +
        geom_point(size = 3, color = "#185A6D") +
        scale_x_continuous(breaks = 2016:2020) +
        labs(x = "Year", y = "Leave-one-out RMSE", title = "Kriging cross-validation error by year") + theme_minimal()
    }
  })
  output$five_year_interpretation <- renderUI({
    if (input$method == "Global Moran's I") {
      x <- global_results |> filter(Parameter == input$parameter)
      n <- sum(x$P_value < 0.05)
      peak <- x[which.max(x$Morans_I), ]
      p(sprintf("Across 2016–2020, %s showed significant global spatial autocorrelation in %s of 5 years. The highest Moran's I was %.3f in %s. A higher Moran's I suggests stronger positive spatial association, but year-to-year differences have not been formally tested.",
                parameter_label(input$parameter), n, peak$Morans_I, peak$Year))
    } else if (input$method == "Local Moran's I (LISA)") {
      x <- lisa_comparison() |> filter(Cluster != "Not significant") |>
        group_by(year) |> summarise(Significant_stations = sum(Stations), .groups = "drop")
      p("The chart compares local cluster and outlier counts for each year. Changes in counts are descriptive; local p-values are unadjusted and comparisons do not establish changes in pollution causes.")
    } else {
      x <- kriging_comparison()
      best <- x[which.min(x$RMSE), ]
      p(sprintf("For %s, the lowest cross-validation RMSE across the five years was %.3f in %s (%s model). RMSE reflects prediction error at sampled stations, not guaranteed accuracy at unsampled locations. Compare RMSE only within the same parameter and units.",
                parameter_label(input$parameter), best$RMSE, best$Year, best$Model))
    }
  })
  
  output$annual_trend <- renderPlot({
    x <- water_data |>
      group_by(year) |>
      summarise(value = median(.data[[input$parameter]], na.rm = TRUE), .groups = "drop")
    ggplot(x, aes(x = year, y = value)) +
      geom_line(color = "#247C8A", linewidth = 1) +
      geom_point(color = "#124D59", size = 3) +
      scale_x_continuous(breaks = 2016:2020) +
      labs(x = "Year", y = paste("Median", parameter_label(input$parameter)),
           title = paste("Annual trend:", parameter_label(input$parameter))) + theme_minimal()
  })
  output$station_trend <- renderPlot({
    req(input$station)
    x <- water_data |> filter(as.character(dws_dataset_id) == input$station) |> arrange(year)
    validate(need(nrow(x) > 0, "No observations for this station."))
    ggplot(x, aes(x = year, y = .data[[input$parameter]])) +
      geom_line(color = "#247C8A", linewidth = 1) +
      geom_point(color = "#124D59", size = 3) +
      scale_x_continuous(breaks = 2016:2020) +
      labs(x = "Year", y = parameter_label(input$parameter),
           title = paste("Station:", x$map_label[1])) + theme_minimal()
  })
  output$station_values <- renderTable({
    req(input$station)
    water_data |>
      filter(as.character(dws_dataset_id) == input$station) |>
      select(year, all_of(input$parameter)) |>
      arrange(year)
  })
  output$explorer_table <- renderDT({
    datatable(water_data, filter = "top", options = list(pageLength = 10, scrollX = TRUE))
  })
}

shinyApp(ui, server)
