#
# library(shiny)
#
# msg <- tryCatch({
#     webr::install("CalibraCurve",
#                   repos = c("https://mpc-bioinformatics.r-universe.dev",
#                             "https://bioc.r-universe.dev",
#                             "https://repo.r-wasm.org"))
#     paste0("Version: ", as.character(packageVersion("CalibraCurve")), "\n",
#            "Exports: ",
#            paste(getNamespaceExports("CalibraCurve"), collapse = ", "))
# }, error = function(e) paste("Install failed:", conditionMessage(e)),
# warning = function(w) paste("Install warning:", conditionMessage(w)))
#
# ui <- fluidPage(h4(R.version.string), verbatimTextOutput("out"))
# server <- function(input, output, session) output$out <- renderText(msg)
# shinyApp(ui, server)



install_error <- NULL

if (identical(R.version$os, "emscripten")) {
    install_error <- tryCatch({
        webr::install("CalibraCurve",
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
