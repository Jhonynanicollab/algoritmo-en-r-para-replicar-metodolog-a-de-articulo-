##################################################################################
proj_sf_path <- tryCatch(system.file("proj", package = "sf"), error = function(e) "")
gdal_sf_path <- tryCatch(system.file("gdal", package = "sf"), error = function(e) "")

if (nzchar(proj_sf_path)) {
  Sys.setenv(PROJ_LIB = proj_sf_path)
  Sys.setenv(PROJ_DATA = proj_sf_path)
}
if (nzchar(gdal_sf_path)) {
  Sys.setenv(GDAL_DATA = gdal_sf_path)
}


ruta_actual    <- Sys.getenv("PATH")
partes         <- strsplit(ruta_actual, ";")[[1]]
partes_limpias <- partes[!grepl("PostgreSQL", partes, ignore.case = TRUE)]
Sys.setenv(PATH = paste(partes_limpias, collapse = ";"))
##################################################################################

#EJECUTAR DESDE AQUI

library(shiny)
library(leaflet)
library(dplyr)
library(readr)
library(stringr)
library(sf)
library(spatstat)
library(raster)
library(classInt)
library(RColorBrewer)


RUTA_POR_DEFECTO <- "RENIPRESS_30-04-2026.csv"   # busca este archivo junto al script
CRS_UTM  <- "EPSG:32719"   # UTM WGS84 zona 19 Sur (igual que el paper)
CRS_WGS  <- "EPSG:4326"

tipos <- c("Hospital / Clínica", "Centro de Salud", "Puesto de Salud",
           "Policlínico", "Otro (odontológico, consultorios, dx, etc.)")

etiquetas_riesgo <- c("Bajo", "Bajo Moderado", "Moderado", "Alto Moderado", "Alto")
paleta_colores    <- c("#1a9850", "#91cf60", "#fee08b", "#fc8d59", "#d73027")


clasificar_tipo <- function(x) {
  case_when(
    str_detect(x, "HOSPITAL")            ~ "Hospital / Clínica",
    str_detect(x, "CENTRO(S)? DE SALUD") ~ "Centro de Salud",
    str_detect(x, "PUESTO")              ~ "Puesto de Salud",
    str_detect(x, "POLICLINICO")         ~ "Policlínico",
    TRUE                                  ~ "Otro (odontológico, consultorios, dx, etc.)"
  )
}


limpiar_datos <- function(datos_raw) {
  
  validate(need(all(c("ESTADO","NORTE","ESTE","CLASIFICACION","NOMBRE") %in%
                      names(datos_raw)),
                "El CSV no tiene las columnas esperadas (ESTADO, NORTE, ESTE, CLASIFICACION, NOMBRE)."))
  
  datos <- datos_raw %>%
    filter(ESTADO == "ACTIVO") %>%
    mutate(
      NORTE = suppressWarnings(as.numeric(NORTE)),
      ESTE  = suppressWarnings(as.numeric(ESTE))
    ) %>%
    filter(!is.na(NORTE), !is.na(ESTE)) %>%
    distinct(NOMBRE, NORTE, ESTE, .keep_all = TRUE)
  
 
  mad_lat <- mad(datos$NORTE, constant = 1.4826)
  mad_lon <- mad(datos$ESTE, constant = 1.4826)
  z_lat <- if (mad_lat > 0) abs(datos$NORTE - median(datos$NORTE)) / mad_lat else 0
  z_lon <- if (mad_lon > 0) abs(datos$ESTE - median(datos$ESTE)) / mad_lon else 0
  es_outlier <- z_lat > 3.5 | z_lon > 3.5
  
  if (any(es_outlier)) {
    showNotification(
      paste0("Se removieron ", sum(es_outlier),
             " registro(s) con coordenadas atípicas: ",
             paste(datos$NOMBRE[es_outlier], collapse = "; ")),
      type = "warning", duration = 10
    )
  }
  datos <- datos[!es_outlier, ]
  
  validate(need(nrow(datos) >= 5,
                "Muy pocos registros válidos tras la limpieza (revisa el CSV)."))
  
  datos %>% mutate(TIPO_EQUIP = clasificar_tipo(CLASIFICACION))
}

# ---- Ancho de banda (fórmula del paper, regla de Silverman/ESRI) ----------
#   SR = 0.9 * min(SD, sqrt(1/ln(2)) * Dm) * n^(-0.2)
calcular_SR <- function(coords) {
  n    <- nrow(coords)
  xbar <- mean(coords[, 1]); ybar <- mean(coords[, 2])
  SD   <- sqrt(sum((coords[, 1] - xbar)^2) / n + sum((coords[, 2] - ybar)^2) / n)
  d    <- sqrt((coords[, 1] - xbar)^2 + (coords[, 2] - ybar)^2)
  Dm   <- median(d)
  0.9 * min(SD, sqrt(1 / log(2)) * Dm) * n^(-0.2)
}

# ---- Densidad de Kernel por categoría -> mapa parcial (raster) ------------
generar_raster_densidad <- function(puntos_utm, tipo, ventana, resolucion) {
  sub <- puntos_utm %>% filter(TIPO_EQUIP == tipo)
  if (nrow(sub) < 2) {
    showNotification(paste0("'", tipo, "' tiene menos de 2 puntos: se omite del análisis."),
                     type = "warning", duration = 8)
    return(NULL)
  }
  coords <- st_coordinates(sub)
  SR   <- calcular_SR(coords)
  pp   <- ppp(coords[, 1], coords[, 2], window = ventana)
  dens <- density.ppp(pp, sigma = SR, eps = resolucion)
  r <- raster(dens)
  crs(r) <- CRS_UTM
  r
}


reclasificar_raster <- function(r, n_clases = 5) {
  if (is.null(r)) return(NULL)
  vals <- values(r)
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return(NULL)
  
  obtener_brks <- function(style) {
    brks <- tryCatch(
      classIntervals(vals, n = n_clases, style = style)$brks,
      error = function(e) NULL
    )
    if (!is.null(brks)) brks <- unique(brks)
    brks
  }
  
  brks <- obtener_brks("quantile")
  if (is.null(brks) || length(brks) < 3) brks <- obtener_brks("equal")
  
  if (is.null(brks) || length(brks) < 3) {
    # No hay variabilidad suficiente: todo el raster va a una sola clase (3 = moderado)
    r_out <- r; values(r_out) <- ifelse(is.na(values(r)), NA, 3)
    return(r_out)
  }
  
  n_reales <- length(brks) - 1
  brks[1] <- -Inf; brks[length(brks)] <- Inf
  m <- matrix(c(brks[-length(brks)], brks[-1], seq_len(n_reales)), ncol = 3)
  reclassify(r, m, include.lowest = TRUE)
}

# ---- Pipeline completo: de datos limpios a mapa de riesgo final -----------
calcular_riesgo <- function(datos_limpios, buffer_m, resolucion) {
  
  puntos_sf  <- st_as_sf(datos_limpios, coords = c("ESTE", "NORTE"), crs = 4326)
  puntos_utm <- st_transform(puntos_sf, crs = CRS_UTM)
  
  hull            <- st_convex_hull(st_union(puntos_utm))
  area_estudio    <- st_buffer(hull, dist = buffer_m)
  area_estudio_sf <- st_sf(id = 1, geometry = area_estudio)
  ventana         <- as.owin(area_estudio_sf)
  
  rasters_parciales <- lapply(tipos, generar_raster_densidad,
                              puntos_utm = puntos_utm, ventana = ventana,
                              resolucion = resolucion)
  names(rasters_parciales) <- tipos
  
  rasters_reclasif <- lapply(rasters_parciales, reclasificar_raster)
  
  rasters_validos <- Filter(Negate(is.null), rasters_reclasif)
  validate(need(length(rasters_validos) >= 1,
                "Ninguna categoría tuvo suficientes puntos para estimar densidad."))
  
  base_r <- rasters_validos[[1]]
  rasters_alineados <- lapply(rasters_validos, function(r) resample(r, base_r, method = "ngb"))
  raster_suma   <- Reduce(`+`, rasters_alineados)
  raster_riesgo <- reclasificar_raster(raster_suma, 5)
  
  # Áreas por clase (m² y %)
  f <- as.data.frame(freq(raster_riesgo, useNA = "no"))
  res_m2 <- prod(res(raster_riesgo))
  f$area_m2 <- round(f$count * res_m2, 2)
  f$pct     <- round(100 * f$area_m2 / sum(f$area_m2), 1)
  f$clase   <- etiquetas_riesgo[f$value]
  tabla_areas <- f[, c("clase", "area_m2", "pct")]
  
  # Reproyección a WGS84 para leaflet
  raster_riesgo_wgs <- projectRaster(raster_riesgo, crs = CRS_WGS, method = "ngb")
  rasters_reclasif_wgs <- lapply(rasters_reclasif, function(r) {
    if (is.null(r)) return(NULL)
    projectRaster(r, crs = CRS_WGS, method = "ngb")
  })
  
 
  area_estudio_wgs <- st_transform(area_estudio_sf, crs = CRS_WGS)
  bbox_estudio <- st_bbox(area_estudio_wgs)
  
  list(
    puntos_sf = puntos_sf,
    rasters_reclasif_wgs = rasters_reclasif_wgs,
    raster_riesgo = raster_riesgo,
    raster_riesgo_wgs = raster_riesgo_wgs,
    tabla_areas = tabla_areas,
    area_estudio_wgs = area_estudio_wgs,
    bbox_estudio = bbox_estudio
  )
}


ui <- fluidPage(
  titlePanel("Zonas de riesgo de contagio - Establecimientos de salud, Puno"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      fileInput("archivo", "Sube el CSV de RENIPRESS/SUSALUD (;-separado):",
                accept = ".csv"),
      helpText(if (file.exists(RUTA_POR_DEFECTO))
        paste0("Se cargó por defecto: ", RUTA_POR_DEFECTO)
        else "No se encontró un archivo por defecto: sube tu CSV arriba."),
      sliderInput("buffer", "Buffer del área de estudio (m):",
                  min = 100, max = 1500, value = 300, step = 50),
      sliderInput("resolucion", "Resolución del raster (m/píxel):",
                  min = 5, max = 100, value = 20, step = 5),
      actionButton("recalcular", "Recalcular mapas", class = "btn-primary"),
      hr(),
      downloadButton("descargar_raster", "Descargar mapa de riesgo (.tif)"),
      br(), br(),
      downloadButton("descargar_tabla", "Descargar tabla de áreas (.csv)")
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("1. Equipamiento de salud (puntos)",
                 br(), leafletOutput("mapa_puntos", height = 600)),
        tabPanel("2. Mapas parciales (Densidad de Kernel)",
                 br(),
                 selectInput("categoria", "Tipo de equipamiento:", choices = tipos),
                 leafletOutput("mapa_parcial", height = 600)),
        tabPanel("3. Mapa de riesgo (Álgebra de mapas)",
                 br(), leafletOutput("mapa_riesgo", height = 600)),
        tabPanel("4. Resultados (áreas)",
                 br(),
                 tableOutput("tabla_resultados"),
                 uiOutput("area_total"))
      )
    )
  )
)


server <- function(input, output, session) {
  
  # ---- Fuente de datos: archivo subido o archivo por defecto ----
  datos_raw <- reactive({
    if (!is.null(input$archivo)) {
      read_delim(input$archivo$datapath, delim = ";",
                 locale = locale(encoding = "UTF-8"), show_col_types = FALSE)
    } else {
      validate(need(file.exists(RUTA_POR_DEFECTO),
                    "Sube un archivo CSV para iniciar el análisis."))
      read_delim(RUTA_POR_DEFECTO, delim = ";",
                 locale = locale(encoding = "UTF-8"), show_col_types = FALSE)
    }
  })
  
  datos_limpios <- reactive({
    limpiar_datos(datos_raw())
  })
  
  # ---- Recalcula solo al presionar el botón (o en la primera carga) -------
  resultado <- eventReactive(list(input$recalcular, datos_limpios()), {
    withProgress(message = "Calculando densidad de Kernel y zonas de riesgo...", {
      calcular_riesgo(datos_limpios(), input$buffer, input$resolucion)
    })
  }, ignoreNULL = FALSE)
  
  paleta_tipos <- colorFactor(brewer.pal(length(tipos), "Set1"), domain = tipos)
  pal_clases   <- colorFactor(paleta_colores, domain = 1:5)
  

  mapa_base <- function(res) {
    bbox <- res$bbox_estudio
    leaflet(options = leafletOptions(minZoom = 3)) %>%
      addProviderTiles("Esri.WorldImagery", group = "Satelital") %>%
      addProviderTiles("OpenStreetMap", group = "Calles (OSM)") %>%
      addProviderTiles("Esri.WorldTopoMap", group = "Topográfico") %>%
      addProviderTiles("CartoDB.PositronOnlyLabels", group = "Satelital") %>%
      addPolygons(
        data = res$area_estudio_wgs,
        fill = FALSE, color = "#FF0000", weight = 3, dashArray = "6 4",
        group = "Área de estudio",
        popup = "Área de estudio - Puno, Perú"
      ) %>%
      addScaleBar(position = "bottomleft") %>%
      fitBounds(bbox[["xmin"]], bbox[["ymin"]], bbox[["xmax"]], bbox[["ymax"]]) %>%
      addLayersControl(
        baseGroups = c("Satelital", "Calles (OSM)", "Topográfico"),
        overlayGroups = "Área de estudio",
        options = layersControlOptions(collapsed = FALSE)
      )
  }
  
  output$mapa_puntos <- renderLeaflet({
    req(resultado())
    mapa_base(resultado()) %>%
      addCircleMarkers(
        data = resultado()$puntos_sf,
        radius = 6, stroke = TRUE, weight = 1, opacity = 1,
        color = "#000000", fillColor = ~paleta_tipos(TIPO_EQUIP),
        fillOpacity = 0.9,
        popup = ~paste0("<b>", NOMBRE, "</b><br>", TIPO_EQUIP)
      ) %>%
      addLegend(pal = paleta_tipos, values = tipos, title = "Tipo de equipamiento")
  })
  
  output$mapa_parcial <- renderLeaflet({
    req(resultado())
    r <- resultado()$rasters_reclasif_wgs[[input$categoria]]
    validate(need(!is.null(r), "Esta categoría no tuvo suficientes puntos para estimar densidad."))
    mapa_base(resultado()) %>%
      addRasterImage(r, colors = pal_clases, opacity = 0.65, group = "Riesgo") %>%
      addLegend(colors = paleta_colores, labels = etiquetas_riesgo,
                title = paste("Riesgo -", input$categoria))
  })
  
  output$mapa_riesgo <- renderLeaflet({
    req(resultado())
    mapa_base(resultado()) %>%
      addRasterImage(resultado()$raster_riesgo_wgs, colors = pal_clases,
                     opacity = 0.65, group = "Riesgo") %>%
      addLegend(colors = paleta_colores, labels = etiquetas_riesgo,
                title = "Riesgo de contagio COVID-19")
  })
  
  output$tabla_resultados <- renderTable({
    req(resultado())
    resultado()$tabla_areas
  })
  
  output$area_total <- renderUI({
    req(resultado())
    p("Área total analizada (m²): ", strong(round(sum(resultado()$tabla_areas$area_m2), 2)))
  })
  
  output$descargar_raster <- downloadHandler(
    filename = function() "mapa_riesgo_salud.tif",
    content = function(file) {
      writeRaster(resultado()$raster_riesgo, file, format = "GTiff", overwrite = TRUE)
    }
  )
  
  output$descargar_tabla <- downloadHandler(
    filename = function() "tabla_areas_riesgo.csv",
    content = function(file) {
      write_csv(resultado()$tabla_areas, file)
    }
  )
}

shinyApp(ui, server)