


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


