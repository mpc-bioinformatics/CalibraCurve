


install_error <- NULL

if (identical(R.version$os, "emscripten")) {
    install_error <- tryCatch({
        webr::install(c("CalibraCurve", "colourpicker"),
                      repos = c("https://mpc-bioinformatics.r-universe.dev",
                                "https://bioc.r-universe.dev",
                                "https://repo.r-wasm.org"))
        NULL
    }, error = function(e) conditionMessage(e))
}

library(shiny)

if (!is.null(install_error)) {
    shinyApp(
        ui = fluidPage(h3("Package installation failed"),
                       verbatimTextOutput("err")),
        server = function(input, output, session) {
            output$err <- renderText(install_error)
        }
    )
} else {
    library("CalibraCurve", character.only = TRUE)
    CalibraCurveApp()   # must be the last expression
}



# library(shiny)
#
# msg <- tryCatch({
#     webr::install("CalibraCurve",
#                   repos = c("https://mpc-bioinformatics.r-universe.dev",
#                             "https://bioc.r-universe.dev",
#                             "https://repo.r-wasm.org"))
#     pkg <- system.file(package = "CalibraCurve")
#     paste(c(
#         paste("Version:", as.character(packageVersion("CalibraCurve"))),
#         paste("system.file extdata:", system.file("extdata", package = "CalibraCurve")),
#         "Files under extdata:",
#         head(list.files(file.path(pkg, "extdata"), recursive = TRUE,
#                         include.dirs = TRUE), 100)
#     ), collapse = "\n")
# }, error = function(e) paste("Failed:", conditionMessage(e)))
#
# ui <- fluidPage(verbatimTextOutput("out"))
# server <- function(input, output, session) output$out <- renderText(msg)
# shinyApp(ui, server)
