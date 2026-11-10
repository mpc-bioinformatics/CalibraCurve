if (identical(R.version$os, "emscripten")) {
    if (!requireNamespace("CalibraCurve", quietly = TRUE)) {
        webr::install("CalibraCurve",
                      repos = c("https://bioc.r-universe.dev",
                                "https://repo.r-wasm.org"))
    }
}

library(shiny)
library(CalibraCurve)

calibraCurveApp()   # must be the last expression, returning the app object
