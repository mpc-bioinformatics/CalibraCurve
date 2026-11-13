#' Upload tab module UI
#'
#' @keywords internal
#' @importFrom tools file_ext
#' @importFrom utils read.table
#' @importFrom openxlsx read.xlsx
.read_uploaded_columns <- function(path, text_separator = "auto", header = TRUE, file_ext = NULL) {
    ext <- if (is.null(file_ext)) tolower(tools::file_ext(path)) else tolower(file_ext)

    if (ext %in% c("csv", "txt", "tsv")) {
        sep_value <- .resolve_effective_separator(path, ext, text_separator, header = header)
        header_cols <- .extract_header_columns(path, sep_value)
        if (length(header_cols) > 1L) {
            return(header_cols)
        }
        tmp <- utils::read.table(
            path,
            header = isTRUE(header),
            sep = sep_value,
            check.names = FALSE,
            stringsAsFactors = FALSE,
            nrows = 10
        )
        return(colnames(tmp))
    }

    if (ext == "xlsx") {
        tmp <- openxlsx::read.xlsx(path, colNames = TRUE, sheet = 1, rowNames = FALSE)
        return(colnames(tmp))
    }

    if (ext == "rds") {
        return(c("amount_fmol", "Substance"))
    }

    character(0)
}

.detect_text_separator <- function(path, fallback = ",") {
    first_lines <- tryCatch(
        readLines(path, n = 5L, warn = FALSE),
        error = function(e) character(0)
    )
    if (length(first_lines) == 0L) {
        return(fallback)
    }

    candidates <- c(",", ";", "\t", "|")
    split_scores <- vapply(candidates, function(sep) {
        sum(vapply(strsplit(first_lines, sep, fixed = TRUE), length, integer(1)))
    }, numeric(1))

    best_idx <- which.max(split_scores)
    best <- candidates[[best_idx]]
    if (length(best) == 0L || is.na(best) || max(split_scores) <= length(first_lines)) {
        return(fallback)
    }
    best
}

.default_text_separator <- function(ext) {
    if (identical(ext, "tsv")) {
        return("\t")
    }
    ","
}

.extract_header_columns <- function(path, separator) {
    header_line <- tryCatch(readLines(path, n = 1L, warn = FALSE), error = function(e) character(0))
    if (length(header_line) == 0L) {
        return(character(0))
    }

    # Remove UTF-8 BOM if present.
    header_line <- sub("^\\ufeff", "", header_line[1])
    parts <- strsplit(header_line, separator, fixed = TRUE)[[1]]
    parts <- trimws(parts)
    parts <- gsub('^"|"$', "", parts)
    parts <- gsub("^'|'$", "", parts)
    parts
}

.select_best_text_separator <- function(path, ext, preferred_separator = NULL, header = TRUE) {
    candidates <- unique(c(preferred_separator, .default_text_separator(ext), ",", ";", "\t", "|"))
    candidates <- candidates[!is.na(candidates) & nzchar(candidates)]

    if (length(candidates) == 0L) {
        return(.default_text_separator(ext))
    }

    scores <- vapply(candidates, function(sep) {
        ncols <- tryCatch({
            tmp <- utils::read.table(
                path,
                header = isTRUE(header),
                sep = sep,
                check.names = FALSE,
                stringsAsFactors = FALSE,
                fill = TRUE,
                quote = "\"'",
                nrows = 10
            )
            ncol(tmp)
        }, error = function(e) 0L)
        as.integer(ncols)
    }, integer(1))

    if (length(scores) == 0L || all(scores <= 0L)) {
        return(.default_text_separator(ext))
    }

    best_idx <- which(scores == max(scores))[1]
    candidates[[best_idx]]
}

.resolve_text_separator <- function(path, ext, text_separator = "auto") {
    if (!identical(text_separator, "auto")) {
        return(text_separator)
    }

    fallback <- .default_text_separator(ext)
    .detect_text_separator(path, fallback = fallback)
}

.resolve_effective_separator <- function(path, ext, text_separator = "auto", header = TRUE) {
    preferred_separator <- .resolve_text_separator(path, ext, text_separator)
    .select_best_text_separator(
        path = path,
        ext = ext,
        preferred_separator = preferred_separator,
        header = header
    )
}

.find_matching_column <- function(choices, target) {
    if (length(choices) == 0L || is.null(target) || !nzchar(target)) {
        return(NULL)
    }
    idx_exact <- which(choices == target)
    if (length(idx_exact) > 0L) {
        return(choices[idx_exact[1L]])
    }
    idx_ci <- which(tolower(choices) == tolower(target))
    if (length(idx_ci) > 0L) {
        return(choices[idx_ci[1L]])
    }
    NULL
}

.locate_example_data_dir <- function(example_dataset) {
    folder_name <- switch(example_dataset,
        "lipidomics" = "SCFA",
        "proteomics" = "MSQC1_xlsx",
        NULL
    )
    if (is.null(folder_name)) {
        stop("Please choose an example dataset.")
    }

    pkg_path <- system.file("extdata", folder_name, package = "CalibraCurve")
    if (nzchar(pkg_path) && dir.exists(pkg_path)) {
        return(pkg_path)
    }

    dev_path <- file.path("inst", "extdata", folder_name)
    if (dir.exists(dev_path)) {
        return(dev_path)
    }

    stop("Could not locate the selected example data folder.")
}

.default_example_columns <- function(example_dataset, choices) {
    conc_target <- if (identical(example_dataset, "proteomics")) "amount" else "Concentration"
    meas_target <- if (identical(example_dataset, "proteomics")) "Area" else "Total_area"

    conc_sel <- .find_matching_column(choices, conc_target)
    meas_sel <- .find_matching_column(choices, meas_target)

    if (is.null(conc_sel) && length(choices) > 0L) {
        conc_sel <- choices[1L]
    }
    if (is.null(meas_sel)) {
        if (length(choices) > 1L) {
            meas_sel <- choices[2L]
        } else if (length(choices) > 0L) {
            meas_sel <- choices[1L]
        }
    }

    list(conc = conc_sel, meas = meas_sel)
}

.import_example_dataset <- function(example_dataset, conc_col, meas_col) {
    example_dir <- .locate_example_data_dir(example_dataset)
    file_type <- if (identical(example_dataset, "proteomics")) "xlsx" else "csv"

    all_files <- list.files(example_dir, full.names = TRUE)
    if (length(all_files) == 0L) {
        stop("The selected example folder does not contain any files.")
    }

    file_paths <- all_files[tolower(tools::file_ext(all_files)) == file_type]
    if (length(file_paths) == 0L) {
        stop("No matching example files were found for the selected dataset.")
    }

    file_names <- basename(file_paths)
    .import_multiple_files(
        paths = file_paths,
        file_names = file_names,
        conc_col = conc_col,
        meas_col = meas_col,
        text_separator = if (identical(file_type, "csv")) "," else "auto",
        decimal_separator = ".",
        header = TRUE
    )
}

#' @importFrom tools file_ext
#' @importFrom stats setNames
.import_single_file <- function(path, conc_col, meas_col, subst_col = NULL,
                               substance_name = NULL,
                               file_ext = NULL,
                               text_separator = "auto",
                               decimal_separator = ".",
                               header = TRUE) {
    ext <- if (is.null(file_ext)) tolower(tools::file_ext(path)) else tolower(file_ext)
    if (ext == "xlsx") {
        result <- readDataTable(
            dataPath = path,
            fileType = "xlsx",
            concCol = conc_col,
            measCol = meas_col,
            substCol = subst_col
        )
    } else if (ext %in% c("csv", "txt", "tsv")) {
        sep_value <- .resolve_effective_separator(path, ext, text_separator, header = header)
        dec_value <- if (identical(decimal_separator, "auto")) "." else decimal_separator
        file_type <- if (ext == "csv") "csv" else "txt"
        result <- readDataTable(
            dataPath = path,
            fileType = file_type,
            concCol = conc_col,
            measCol = meas_col,
            substCol = subst_col,
            sep = sep_value,
            dec = dec_value,
            header = header
        )
    } else if (ext == "rds") {
        result <- readDataSE(
            dataPath = path,
            concColName = "amount_fmol",
            substColName = "Substance",
            assayNumber = 1
        )
    } else {
        stop("Unsupported file type. Use CSV, TSV, TXT, XLSX, or RDS.")
    }

    if (is.data.frame(result)) {
        if (is.null(substance_name)) {
            substance_name <- "substance"
        }
        return(stats::setNames(list(result), substance_name))
    }
    if (is.list(result) && length(result) > 0L) {
        return(result)
    }
    stop("The uploaded data did not contain any usable calibration rows.")
}

#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom stats setNames
.import_multiple_files <- function(paths, file_names, conc_col, meas_col,
                                 text_separator = "auto", decimal_separator = ".",
                                 header = TRUE) {
    if (length(paths) == 0L) {
        stop("Please select at least one file.")
    }

    exts <- unique(tolower(tools::file_ext(file_names)))
    if (length(exts) != 1L) {
        stop("Please upload files with the same format only.")
    }

    file_type <- if (exts[1] == "xlsx") "xlsx" else if (exts[1] == "csv") "csv" else "txt"

    result_list <- lapply(seq_along(paths), function(i) {
        path <- paths[[i]]
        file_name <- tools::file_path_sans_ext(basename(file_names[[i]]))
        if (file_type == "xlsx") {
            dat <- readDataTable(
                dataPath = path,
                fileType = file_type,
                concCol = conc_col,
                measCol = meas_col
            )
        } else {
            current_ext <- tolower(tools::file_ext(file_names[[i]]))
            sep_value <- .resolve_effective_separator(path, current_ext, text_separator, header = header)
            dec_value <- if (identical(decimal_separator, "auto")) "." else decimal_separator
            dat <- readDataTable(
                dataPath = path,
                fileType = file_type,
                concCol = conc_col,
                measCol = meas_col,
                sep = sep_value,
                dec = dec_value,
                header = header
            )
        }
        if (!is.data.frame(dat)) {
            stop("The uploaded file '", file_name, "' did not yield a valid calibration table.")
        }
        stats::setNames(list(dat), file_name)
    })

    unlist(result_list, recursive = FALSE)
}

.format_import_error_message <- function(error, file_type = NULL) {
    paste0(
        "Import failed: Please check if the correct concentration and measurement columns are selected, adjust if needed and press 'Import data' again. ",
        "For .tsv, .txt and .csv files only, a wrong decimal separator or text separator can also cause this problem. ",
        "Please adjust the advanced upload settings if necessary."
    )
}

.sanitize_export_name <- function(x) {
    x <- gsub("[[:space:]]+", "_", x)
    x <- gsub("[^A-Za-z0-9_.-]", "_", x)
    x <- gsub("_+", "_", x)
    x <- sub("^_", "", x)
    x <- sub("_$", "", x)
    if (!nzchar(x)) {
        x <- "substance"
    }
    x
}

.extract_first_ggplot <- function(x) {
    if (inherits(x, "ggplot")) {
        return(x)
    }
    if (is.list(x)) {
        for (i in seq_along(x)) {
            found <- .extract_first_ggplot(x[[i]])
            if (!is.null(found)) {
                return(found)
            }
        }
    }
    NULL
}

.select_calibration_dataset <- function(calibration_result, selected_name = NULL) {
    if (is.null(calibration_result) || is.null(calibration_result$RES)) {
        return(calibration_result)
    }

    res_names <- names(calibration_result$RES)
    if (length(res_names) == 0L) {
        return(calibration_result)
    }
    if (is.null(selected_name) || !selected_name %in% res_names) {
        selected_name <- res_names[1L]
    }

    selected_idx <- match(selected_name, res_names)
    summary_tab <- calibration_result$summary_tab
    if (!is.null(summary_tab) && !is.null(summary_tab$substance)) {
        summary_tab <- summary_tab[summary_tab$substance %in% selected_name, , drop = FALSE]
    }

    plot_CC_list <- calibration_result$plot_CC_list
    if (!is.null(plot_CC_list) && length(plot_CC_list) > 1L) {
        plot_CC_list <- plot_CC_list[selected_idx]
    }

    plot_RF_list <- calibration_result$plot_RF_list
    if (!is.null(plot_RF_list) && length(plot_RF_list) > 1L) {
        plot_RF_list <- plot_RF_list[selected_idx]
    }

    list(
        RES = calibration_result$RES[selected_idx],
        summary_tab = summary_tab,
        plot_CC_list = plot_CC_list,
        plot_RF_list = plot_RF_list
    )
}

.default_export_plot_dimensions <- function(plot_type = "single_plots") {
    if (!is.null(plot_type) && plot_type %in% c("multiplot", "all_in_one")) {
        return(list(width = 25, height = 15))
    }
    list(width = 20, height = 10)
}

.plot_for_export <- function(plot_list, idx) {
    if (is.null(plot_list) || length(plot_list) == 0L) {
        return(NULL)
    }
    if (inherits(plot_list, "ggplot")) {
        return(plot_list)
    }
    if (!is.list(plot_list)) {
        return(plot_list)
    }
    if (length(plot_list) == 1L) {
        return(plot_list[[1L]])
    }
    if (idx > length(plot_list)) {
        return(NULL)
    }
    plot_list[[idx]]
}

.calibration_dataset_names <- function(calibration_result) {
    if (is.null(calibration_result) || is.null(calibration_result$RES)) {
        return(character(0))
    }
    names(calibration_result$RES)
}

.calibration_selection_tables <- function(calibration_result) {
    if (is.null(calibration_result) || is.null(calibration_result$RES) ||
        length(calibration_result$RES) == 0L) {
        return(list())
    }

    current_result <- calibration_result$RES[[1]]
    tables <- list()
    if (!is.null(current_result$result_table_conc_levels)) {
        tables[["Concentration levels"]] <- current_result$result_table_conc_levels
    }
    if (!is.null(current_result$result_table_obs)) {
        tables[["Observation-level results"]] <- current_result$result_table_obs
    }
    tables
}


#' @importFrom openxlsx write.xlsx
#' @importFrom utils write.csv zip
#' @importFrom grDevices pdf png dev.off
.export_calibration_results_zip <- function(calibration_result,
    zip_path = tempfile(fileext = ".zip"),
    table_format = c("xlsx", "csv"),
    plot_format = c("png", "pdf"),
    plot_type = NULL,
    plot_width = 14,
    plot_height = 10,
    plot_dpi = 300,
    quiet_zip = FALSE) {
    table_format <- match.arg(table_format)
    plot_format <- match.arg(plot_format)

    if (is.null(calibration_result) || is.null(calibration_result$RES)) {
        stop("No calibration result available to export.")
    }

    tmp_dir <- tempfile(pattern = "calibra_curve_export_")
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

    summary_tab <- calibration_result$summary_tab
    if (!is.null(summary_tab)) {
        if (table_format == "xlsx") {
            openxlsx::write.xlsx(summary_tab, file = file.path(tmp_dir, "summary_tab.xlsx"))
        } else {
            utils::write.csv(summary_tab, file = file.path(tmp_dir, "summary_tab.csv"), row.names = FALSE)
        }
    }

    substance_names <- names(calibration_result$RES)
    if (is.null(substance_names)) {
        substance_names <- paste0("substance_", seq_along(calibration_result$RES))
    }
    substance_names[!nzchar(substance_names)] <- paste0("substance_", which(!nzchar(substance_names)))
    substance_export_names <- make.unique(vapply(substance_names, .sanitize_export_name, character(1)))
    combined_plot_mode <- !is.null(plot_type) && plot_type %in% c("multiplot", "all_in_one")

    for (i in seq_along(calibration_result$RES)) {
        substance_suffix <- substance_export_names[i]
        res_obj <- calibration_result$RES[[i]]

        if (!is.null(res_obj$result_table_conc_levels)) {
            if (table_format == "xlsx") {
                openxlsx::write.xlsx(
                    res_obj$result_table_conc_levels,
                    file = file.path(tmp_dir, paste0("result_table_conc_levels_", substance_suffix, ".xlsx"))
                )
            } else {
                write.csv(
                    res_obj$result_table_conc_levels,
                    file = file.path(tmp_dir, paste0("result_table_conc_levels_", substance_suffix, ".csv")),
                    row.names = FALSE
                )
            }
        }
        if (!is.null(res_obj$result_table_obs)) {
            if (table_format == "xlsx") {
                openxlsx::write.xlsx(
                    res_obj$result_table_obs,
                    file = file.path(tmp_dir, paste0("result_table_obs_", substance_suffix, ".xlsx"))
                )
            } else {
                write.csv(
                    res_obj$result_table_obs,
                    file = file.path(tmp_dir, paste0("result_table_obs_", substance_suffix, ".csv")),
                    row.names = FALSE
                )
            }
        }

        if (!combined_plot_mode && !is.null(calibration_result$plot_CC_list)) {
            CC_plot <- .plot_for_export(calibration_result$plot_CC_list, i)
            if (!is.null(CC_plot)) {
                plot_file <- file.path(tmp_dir, paste0("calibration_curve_", substance_suffix, ".", plot_format))
                if (plot_format == "pdf") {
                    grDevices::pdf(file = plot_file, width = plot_width, height = plot_height)
                    print(CC_plot)
                    grDevices::dev.off()
                } else {
                    grDevices::png(file = plot_file, width = plot_width * plot_dpi / 2.54,
                        height = plot_height * plot_dpi / 2.54, units = "px",
                        res = plot_dpi)
                    print(CC_plot)
                    grDevices::dev.off()
                }
            }
        }

        if (!combined_plot_mode && !is.null(calibration_result$plot_RF_list)) {
            RF_plot <- .plot_for_export(calibration_result$plot_RF_list, i)
            if (!is.null(RF_plot)) {
                plot_file <- file.path(tmp_dir, paste0("response_factor_plot_", substance_suffix, ".", plot_format))
                if (plot_format == "pdf") {
                    grDevices::pdf(file = plot_file, width = plot_width, height = plot_height)
                    print(RF_plot)
                    grDevices::dev.off()
                } else {
                    grDevices::png(file = plot_file, width = plot_width * plot_dpi / 2.54,
                        height = plot_height * plot_dpi / 2.54, units = "px",
                        res = plot_dpi)
                    print(RF_plot)
                    grDevices::dev.off()
                }
            }
        }
    }

    if (combined_plot_mode) {
        CC_plot <- .plot_for_export(calibration_result$plot_CC_list, 1L)
        if (!is.null(CC_plot)) {
            plot_file <- file.path(tmp_dir, paste0("calibration_curve.", plot_format))
            if (plot_format == "pdf") {
                grDevices::pdf(file = plot_file, width = plot_width, height = plot_height)
                print(CC_plot)
                grDevices::dev.off()
            } else {
                grDevices::png(file = plot_file, width = plot_width * plot_dpi / 2.54,
                    height = plot_height * plot_dpi / 2.54, units = "px",
                    res = plot_dpi)
                print(CC_plot)
                grDevices::dev.off()
            }
        }

        RF_plot <- .plot_for_export(calibration_result$plot_RF_list, 1L)
        if (!is.null(RF_plot)) {
            plot_file <- file.path(tmp_dir, paste0("response_factor_plot.", plot_format))
            if (plot_format == "pdf") {
                grDevices::pdf(file = plot_file, width = plot_width, height = plot_height)
                print(RF_plot)
                grDevices::dev.off()
            } else {
                grDevices::png(file = plot_file, width = plot_width * plot_dpi / 2.54,
                    height = plot_height * plot_dpi / 2.54, units = "px",
                    res = plot_dpi)
                print(RF_plot)
                grDevices::dev.off()
            }
        }
    }

    files_to_zip <- list.files(tmp_dir, recursive = TRUE, full.names = FALSE)
    if (length(files_to_zip) == 0L) {
        stop("No files were produced for export.")
    }

    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(tmp_dir)
    zip_flags <- if (isTRUE(quiet_zip)) "-rq9X" else "-r9X"
    utils::zip(zipfile = zip_path, files = files_to_zip, flags = zip_flags)
    invisible(zip_path)
}

normalize_uploaded_dataset <- function(x) {
    if (is.data.frame(x)) {
        return(x)
    }
    if (is.list(x) && length(x) > 0L) {
        frames <- Filter(function(z) is.data.frame(z), x)
        if (length(frames) == 0L) {
            stop("Uploaded data did not match the expected CalibraCurve table format.")
        }
        if (length(frames) == 1L) {
            return(frames[[1]])
        }

        combined <- lapply(seq_along(frames), function(i) {
            df <- frames[[i]]
            if (nrow(df) == 0L) {
                return(NULL)
            }
            rownames(df) <- NULL
            df$Dataset <- rep(if (is.null(names(frames)[i])) paste0("Dataset_", i) else names(frames)[i], nrow(df))
            df
        })
        combined <- Filter(Negate(is.null), combined)
        if (length(combined) == 0L) {
            stop("Uploaded data did not match the expected CalibraCurve table format.")
        }
        return(do.call(rbind, combined))
    }
    if (is.matrix(x)) {
        return(as.data.frame(x, stringsAsFactors = FALSE))
    }
    stop("Uploaded data did not match the expected CalibraCurve table format.")
}


#' @importFrom shiny NS tagList fluidRow column radioButtons conditionalPanel
#' @importFrom shiny fileInput helpText checkboxInput uiOutput actionButton
#' @importFrom shiny wellPanel div textOutput tableOutput strong
mod_data_upload_ui <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::fluidRow(
            shiny::column(
                width = 4,
                shiny::radioButtons(
                    ns("data_source"),
                    "Data source",
                    choices = c(
                        "Upload your own file(s)" = "upload",
                        "Use example data" = "example"
                    ),
                    selected = "upload"
                ),
                shiny::conditionalPanel(
                    condition = "input.data_source == 'example'",
                    ns = ns,
                    shiny::selectInput(
                        ns("example_dataset"),
                        "Example dataset",
                        choices = c(
                            "Lipidomics (SCFA)" = "lipidomics",
                            "Proteomics (MSQC1_xlsx)" = "proteomics"
                        ),
                        selected = "lipidomics"
                    )
                ),
                shiny::conditionalPanel(
                    condition = "input.data_source == 'upload'",
                    ns = ns,
                    shiny::radioButtons(
                    ns("upload_mode"),
                    "Choose upload type",
                    choices = c(
                        "Single file, single substance" = "single_single",
                        "Single file, multiple substances" = "single_multiple",
                        "Multiple files" = "multiple_files"
                    ),
                    selected = "single_single"
                )
                ),
                shiny::conditionalPanel(
                    condition = "input.data_source == 'upload' && (input.upload_mode == 'single_single' || input.upload_mode == 'single_multiple')",
                    ns = ns,
                    shiny::fileInput(
                        ns("data_file"),
                        "Upload calibration file",
                        accept = c(".csv", ".tsv", ".txt", ".xlsx", ".rds"),
                        multiple = FALSE
                    )
                ),
                shiny::conditionalPanel(
                    condition = "input.data_source == 'upload' && input.upload_mode == 'multiple_files'",
                    ns = ns,
                    shiny::fileInput(
                        ns("data_files"),
                        "Upload calibration files",
                        accept = c(".csv", ".tsv", ".txt", ".xlsx"),
                        multiple = TRUE
                    )
                ),
                shiny::helpText("Supported upload file types: .csv, .tsv, .txt, .xlsx, .rds."),
                shiny::checkboxInput(ns("show_advanced_settings"), "Advanced upload settings", value = FALSE),
                shiny::uiOutput(ns("advanced_upload_controls")),
                shiny::uiOutput(ns("upload_controls")),
                shiny::actionButton(ns("import_data"), "Import data")
            ),
            shiny::column(
                width = 8,
                shiny::wellPanel(
                    shiny::helpText("The data are automatically re-ordered by increasing concentration before preview and calibration."),
                    shiny::uiOutput(ns("status")),
                    shiny::div(
                        style = "display: flex; gap: 8px; align-items: center;",
                        shiny::actionButton(ns("prev_page"), "Previous 10 rows"),
                        shiny::actionButton(ns("next_page"), "Next 10 rows")
                    ),
                    shiny::textOutput(ns("preview_info")),
                    shiny::tableOutput(ns("preview")),
                    shiny::strong("Please check the data preview. If everything looks correct, click 'Next' to continue to the calibration tab."),
                    shiny::uiOutput(ns("next_button"))
                )
            )
        )
    )
}

#' Upload tab module server
#'
#' @keywords internal
mod_data_upload_server <- function(id) {
    shiny::moduleServer(id, function(input, output, session) {
        preview_page <- shiny::reactiveVal(1L)
        import_status <- shiny::reactiveVal("No data imported yet.")

        column_choices <- shiny::reactive({
            if (identical(input$data_source, "example")) {
                example_dataset <- if (is.null(input$example_dataset)) "lipidomics" else input$example_dataset
                example_dir <- .locate_example_data_dir(example_dataset)
                file_type <- if (identical(example_dataset, "proteomics")) "xlsx" else "csv"
                paths <- list.files(example_dir, full.names = TRUE)
                paths <- paths[tolower(tools::file_ext(paths)) == file_type]
                if (length(paths) == 0L) {
                    return(character(0))
                }
                return(.read_uploaded_columns(
                    path = paths[1],
                    text_separator = if (identical(file_type, "csv")) "," else "auto",
                    header = TRUE,
                    file_ext = file_type
                ))
            }

            mode <- input$upload_mode
            if (is.null(mode)) {
                return(character(0))
            }

            file_type <- file_type_for_selected_upload()
            text_sep <- if (is.null(input$text_separator)) .default_text_separator(file_type) else input$text_separator

            if (mode %in% c("single_single", "single_multiple")) {
                if (is.null(input$data_file)) {
                    return(character(0))
                }
                return(.read_uploaded_columns(
                    path = input$data_file$datapath,
                    text_separator = text_sep,
                    header = TRUE,
                    file_ext = tools::file_ext(input$data_file$name)
                ))
            }

            if (is.null(input$data_files)) {
                return(character(0))
            }
            if (length(input$data_files$datapath) == 0L) {
                return(character(0))
            }
            .read_uploaded_columns(
                path = input$data_files$datapath[1],
                text_separator = text_sep,
                header = TRUE,
                file_ext = tools::file_ext(input$data_files$name[1])
            )
        })

        file_type_for_selected_upload <- shiny::reactive({
            if (identical(input$data_source, "example")) {
                if (identical(input$example_dataset, "proteomics")) {
                    return("xlsx")
                }
                return("csv")
            }

            mode <- input$upload_mode
            if (is.null(mode)) {
                return("unknown")
            }

            if (mode %in% c("single_single", "single_multiple")) {
                if (is.null(input$data_file)) {
                    return("unknown")
                }
                return(tolower(tools::file_ext(input$data_file$name)))
            }

            if (is.null(input$data_files) || length(input$data_files$name) == 0L) {
                return("unknown")
            }
            tolower(tools::file_ext(input$data_files$name[1]))
        })

        output$advanced_upload_controls <- shiny::renderUI({
            if (!isTRUE(input$show_advanced_settings)) {
                return(NULL)
            }

            file_type <- file_type_for_selected_upload()
            if (file_type %in% c("csv", "txt", "tsv")) {
                default_text_separator <- .default_text_separator(file_type)
                return(shiny::tagList(
                    shiny::selectInput(
                        session$ns("decimal_separator"),
                        "Decimal separator",
                        choices = c("Point (.)" = ".", "Comma (, )" = ","),
                        selected = "."
                    ),
                    shiny::selectInput(
                        session$ns("text_separator"),
                        "Text separator",
                        choices = c("Comma" = ",", "Semicolon" = ";", "Tab" = "\t", "Pipe" = "|"),
                        selected = default_text_separator
                    )
                ))
            }

            shiny::helpText("No additional upload parsing options are needed for this file type.")
        })

        output$upload_controls <- shiny::renderUI({
            if (identical(input$data_source, "example")) {
                choices <- column_choices()
                example_dataset <- if (is.null(input$example_dataset)) "lipidomics" else input$example_dataset
                defaults <- .default_example_columns(example_dataset, choices)

                return(shiny::tagList(
                    shiny::selectInput(
                        session$ns("conc_col_example"),
                        "Concentration column",
                        choices = choices,
                        selected = defaults$conc
                    ),
                    shiny::selectInput(
                        session$ns("meas_col_example"),
                        "Measurement column",
                        choices = choices,
                        selected = defaults$meas
                    )
                ))
            }

            mode <- input$upload_mode
            if (is.null(mode)) {
                mode <- "single_single"
            }
            choices <- column_choices()

            if (mode == "single_single") {
                return(shiny::tagList(
                    shiny::selectInput(
                        session$ns("conc_col"),
                        "Concentration column",
                        choices = choices,
                        selected = if (length(choices) > 0L) choices[1] else NULL
                    ),
                    shiny::selectInput(
                        session$ns("meas_col"),
                        "Measurement column",
                        choices = choices,
                        selected = if (length(choices) > 1L) choices[2] else if (length(choices) > 0L) choices[1] else NULL
                    )
                ))
            }

            if (mode == "single_multiple") {
                return(shiny::tagList(
                    shiny::selectInput(
                        session$ns("conc_col"),
                        "Concentration column",
                        choices = choices,
                        selected = if (length(choices) > 0L) choices[1] else NULL
                    ),
                    shiny::selectInput(
                        session$ns("meas_col"),
                        "Measurement column",
                        choices = choices,
                        selected = if (length(choices) > 1L) choices[2] else if (length(choices) > 0L) choices[1] else NULL
                    ),
                    shiny::selectInput(
                        session$ns("subst_col"),
                        "Substance column",
                        choices = choices,
                        selected = if (length(choices) > 2L) choices[3] else if (length(choices) > 0L) choices[1] else NULL
                    )
                ))
            }

            shiny::tagList(
                shiny::selectInput(
                    session$ns("conc_col_multi"),
                    "Concentration column",
                    choices = choices,
                    selected = if (length(choices) > 0L) choices[1] else NULL
                ),
                shiny::selectInput(
                    session$ns("meas_col_multi"),
                    "Measurement column",
                    choices = choices,
                    selected = if (length(choices) > 1L) choices[2] else if (length(choices) > 0L) choices[1] else NULL
                )
            )
        })

        uploaded_data <- shiny::eventReactive(input$import_data, {
            if (identical(input$data_source, "example")) {
                return(tryCatch({
                    shiny::req(input$example_dataset)
                    shiny::req(input$conc_col_example)
                    shiny::req(input$meas_col_example)
                    result <- .import_example_dataset(
                        example_dataset = input$example_dataset,
                        conc_col = input$conc_col_example,
                        meas_col = input$meas_col_example
                    )
                    preview_page(1L)
                    import_status("Example data import successful.")
                    result
                }, error = function(e) {
                    preview_page(1L)
                    import_status(.format_import_error_message(e))
                    NULL
                }))
            }

            mode <- input$upload_mode
            if (is.null(mode)) {
                stop("Please choose an upload mode first.")
            }

            tryCatch({
                if (mode == "single_single") {
                    shiny::req(input$data_file)
                    file_type <- tolower(tools::file_ext(input$data_file$name))
                    text_sep <- if (is.null(input$text_separator)) .default_text_separator(file_type) else input$text_separator
                    dec_sep <- if (is.null(input$decimal_separator)) "." else input$decimal_separator
                    result <- .import_single_file(
                        path = input$data_file$datapath,
                        conc_col = input$conc_col,
                        meas_col = input$meas_col,
                        subst_col = NULL,
                        substance_name = tools::file_path_sans_ext(basename(input$data_file$name)),
                        file_ext = tools::file_ext(input$data_file$name),
                        text_separator = if (file_type %in% c("csv", "txt", "tsv")) text_sep else ",",
                        decimal_separator = if (file_type %in% c("csv", "txt", "tsv")) dec_sep else ".",
                        header = TRUE
                    )
                    preview_page(1L)
                    import_status("Import successful.")
                    return(result)
                }

                if (mode == "single_multiple") {
                    shiny::req(input$data_file)
                    file_type <- tolower(tools::file_ext(input$data_file$name))
                    text_sep <- if (is.null(input$text_separator)) .default_text_separator(file_type) else input$text_separator
                    dec_sep <- if (is.null(input$decimal_separator)) "." else input$decimal_separator
                    result <- .import_single_file(
                        path = input$data_file$datapath,
                        conc_col = input$conc_col,
                        meas_col = input$meas_col,
                        subst_col = input$subst_col,
                        substance_name = tools::file_path_sans_ext(basename(input$data_file$name)),
                        file_ext = tools::file_ext(input$data_file$name),
                        text_separator = if (file_type %in% c("csv", "txt", "tsv")) text_sep else ",",
                        decimal_separator = if (file_type %in% c("csv", "txt", "tsv")) dec_sep else ".",
                        header = TRUE
                    )
                    preview_page(1L)
                    import_status("Import successful.")
                    return(result)
                }

                shiny::req(input$data_files)
                file_type <- tolower(tools::file_ext(input$data_files$name[1]))
                text_sep <- if (is.null(input$text_separator)) .default_text_separator(file_type) else input$text_separator
                dec_sep <- if (is.null(input$decimal_separator)) "." else input$decimal_separator
                result <- .import_multiple_files(
                    paths = input$data_files$datapath,
                    file_names = input$data_files$name,
                    conc_col = input$conc_col_multi,
                    meas_col = input$meas_col_multi,
                    text_separator = if (file_type %in% c("csv", "txt", "tsv")) text_sep else ",",
                    decimal_separator = if (file_type %in% c("csv", "txt", "tsv")) dec_sep else ".",
                    header = TRUE
                )
                preview_page(1L)
                import_status("Import successful.")
                return(result)
            }, error = function(e) {
                preview_page(1L)
                import_status(.format_import_error_message(e))
                NULL
            })
        })

        shiny::observeEvent(input$prev_page, {
            if (preview_page() > 1L) {
                preview_page(preview_page() - 1L)
            }
        })

        shiny::observeEvent(input$next_page, {
            dat <- uploaded_data()
            if (is.null(dat)) {
                return(NULL)
            }
            dat <- normalize_uploaded_dataset(dat)
            max_page <- max(1L, ceiling(nrow(dat) / 10))
            if (preview_page() < max_page) {
                preview_page(preview_page() + 1L)
            }
        })

        output$status <- shiny::renderUI({
            status_text <- import_status()
            if (grepl("^Import failed", status_text, ignore.case = TRUE)) {
                shiny::div(style = "color: #b22222; font-weight: bold;", status_text)
            } else {
                shiny::div(status_text)
            }
        })

        output$preview_info <- shiny::renderText({
            dat <- uploaded_data()
            if (is.null(dat)) {
                return("No data available for preview.")
            }
            dat <- normalize_uploaded_dataset(dat)
            start_row <- (preview_page() - 1L) * 10L + 1L
            end_row <- min(nrow(dat), start_row + 9L)
            paste0("Showing rows ", start_row, " to ", end_row, " of ", nrow(dat), ".")
        })

        output$next_button <- shiny::renderUI({
            data_is_valid <- !is.null(uploaded_data())
            outer_style <- if (data_is_valid) {
                "margin-top: 12px;"
            } else {
                "margin-top: 12px; pointer-events: none; opacity: 0.6;"
            }
            button_style <- if (data_is_valid) {
                "background-color: #2e7d32; color: white;"
            } else {
                "background-color: #d9d9d9; color: #666;"
            }
            shiny::div(
                style = outer_style,
                shiny::actionButton(
                    inputId = session$ns("go_to_calibration"),
                    label = "Next",
                    style = button_style
                )
            )
        })

        shiny::observeEvent(input$go_to_calibration, {
            if (is.null(uploaded_data())) {
                return(NULL)
            }
            shiny::showNotification("Please check the data preview and continue to the calibration tab.", type = "message")
        })

        output$preview <- shiny::renderTable({
            dat <- uploaded_data()
            if (is.null(dat)) {
                return(NULL)
            }
            dat <- normalize_uploaded_dataset(dat)
            start_row <- (preview_page() - 1L) * 10L + 1L
            end_row <- min(nrow(dat), start_row + 9L)
            dat[start_row:end_row, , drop = FALSE]
        }, rownames = TRUE)

        list(
            data = uploaded_data,
            next_clicked = shiny::reactive(input$go_to_calibration)
        )
    })
}

#' Calibration tab module UI
#'
#' @keywords internal
mod_calibration_ui <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::fluidRow(
            shiny::column(
                width = 4,
                shiny::wellPanel(
                    shiny::uiOutput(ns("plot_type_ui")),
                    shiny::checkboxInput(ns("show_advanced_plotting"), "Advanced plotting settings", value = FALSE),
                    shiny::conditionalPanel(
                        condition = "input.show_advanced_plotting == true",
                        ns = ns,
                        shiny::tagList(
                            shiny::checkboxInput(ns("show_regression_info"), "Show regression info", FALSE),
                            shiny::checkboxInput(ns("show_linear_range"), "Show linear range", TRUE),
                            shiny::checkboxInput(ns("show_data_points"), "Show data points", TRUE),
                            shiny::checkboxInput(ns("plot_zeros"), "Show blank panel", FALSE),
                            shiny::conditionalPanel(
                                condition = "input.plot_type == 'multiplot'",
                                ns = ns,
                                shiny::selectInput(
                                    ns("multiplot_scales"),
                                    "Multiplot axis scaling",
                                    choices = c(
                                        "Fix x and y axis" = "fixed",
                                        "Fix x axis only" = "free_y",
                                        "Fix y axis only" = "free_x",
                                        "Free x and y axis" = "free"
                                    ),
                                    selected = "free"
                                )
                            ),
                            shiny::selectInput(
                                ns("curve_colour"),
                                "Curve colour",
                                choices = c("red" = "red", "blue" = "blue", "black" = "black", "green" = "green"),
                                selected = "red"
                            ),
                            shiny::selectInput(
                                ns("point_colour"),
                                "Point colour",
                                choices = c("black" = "black", "red" = "red", "blue" = "blue", "green" = "green"),
                                selected = "black"
                            ),
                            shiny::selectInput(
                                ns("linear_range_colour"),
                                "Linear range colour",
                                choices = c("black" = "black", "grey" = "grey", "blue" = "blue", "red" = "red"),
                                selected = "black"
                            ),
                            shiny::checkboxInput(ns("log10_x_axis"), "Use log10 x-axis", TRUE),
                            shiny::selectInput(
                                ns("log_tick_label_style"),
                                "Log10 tick label style",
                                choices = c("Auto" = "auto", "Math" = "math", "Plain" = "plain"),
                                selected = "auto"
                            ),
                            shiny::numericInput(
                                ns("base_size"),
                                "Base font size",
                                value = 11,
                                min = 1,
                                step = 1
                            ),
                            shiny::textInput(
                                ns("RF_colour_threshold"),
                                "RF threshold line colour",
                                value = "orange"
                            ),
                            shiny::textInput(
                                ns("RF_colour_within"),
                                "RF colour within thresholds",
                                value = "#00BFC4"
                            ),
                            shiny::textInput(
                                ns("RF_colour_outside"),
                                "RF colour outside thresholds",
                                value = "#F8766D"
                            )
                        )
                    ),
                    shiny::checkboxInput(ns("show_advanced_algorithm"), "Advanced algorithm settings", value = FALSE),
                    shiny::conditionalPanel(
                        condition = "input.show_advanced_algorithm == true",
                        ns = ns,
                        shiny::tagList(
                            shiny::numericInput(
                                ns("minReplicates"),
                                "Minimum replicates",
                                value = 3,
                                min = 1,
                                step = 1
                            ),
                            shiny::numericInput(
                                ns("cvThres"),
                                "CV threshold (%)",
                                value = 20,
                                min = 0,
                                step = 1
                            ),
                            shiny::checkboxInput(
                                ns("calcContinuousPrelimRanges"),
                                "Require continuous preliminary range",
                                FALSE
                            ),
                            shiny::selectInput(
                                ns("weightingMethod"),
                                "Weighting method",
                                choices = c("1/x^2" = "1/x^2", "1/x" = "1/x", "None" = "None"),
                                selected = "1/x^2"
                            ),
                            shiny::selectInput(
                                ns("centralTendencyMeasure"),
                                "Bias central tendency",
                                choices = c("mean" = "mean", "median" = "median"),
                                selected = "mean"
                            ),
                            shiny::numericInput(
                                ns("perBiasThres"),
                                "Percent bias threshold (%)",
                                value = 20,
                                min = 0,
                                step = 1
                            ),
                            shiny::checkboxInput(ns("considerPerBiasCV"), "Consider per-bias CV", TRUE),
                            shiny::numericInput(
                                ns("perBiasDistThres"),
                                "Bias distance threshold (%)",
                                value = 10,
                                min = 0,
                                step = 1
                            ),
                            shiny::numericInput(
                                ns("RfThresL"),
                                "Lower RF threshold (%)",
                                value = 80,
                                min = 0,
                                step = 1
                            ),
                            shiny::numericInput(
                                ns("RfThresU"),
                                "Upper RF threshold (%)",
                                value = 120,
                                min = 0,
                                step = 1
                            )
                        )
                    ),
                    shiny::actionButton(ns("run_calibration"), "Run calibration")
                )
            ),
            shiny::column(
                width = 8,
                shiny::wellPanel(
                    style = "background-color: #f5f5f5; min-height: 900px;",
                    shiny::uiOutput(ns("dataset_selector")),
                    shiny::div(
                        style = "display: flex; align-items: center; gap: 12px; margin-bottom: 12px;",
                        shiny::downloadButton(ns("download_results"), "Download all results as ZIP"),
                        shiny::selectInput(
                            ns("table_export_format"),
                            "Table output format",
                            choices = c("XLSX" = "xlsx", "CSV" = "csv"),
                            selected = "xlsx",
                            width = "180px"
                        ),
                        shiny::selectInput(
                            ns("plot_export_format"),
                            "Plot output format",
                            choices = c("PNG" = "png", "PDF" = "pdf"),
                            selected = "png",
                            width = "180px"
                        )
                    ),
                    shiny::div(
                        style = "font-size: 12px; color: #666; margin-top: -6px; margin-bottom: 8px;",
                        "This may take up to a minute."
                    ),
                    shiny::checkboxInput(ns("show_advanced_export_settings"), "Advanced export settings", value = FALSE),
                    shiny::conditionalPanel(
                        condition = "input.show_advanced_export_settings == true",
                        ns = ns,
                        shiny::tagList(
                            shiny::numericInput(
                                ns("export_plot_width"),
                                "Plot width (cm)",
                                value = 14,
                                min = 1,
                                step = 1
                            ),
                            shiny::numericInput(
                                ns("export_plot_height"),
                                "Plot height (cm)",
                                value = 10,
                                min = 1,
                                step = 1
                            ),
                            shiny::numericInput(
                                ns("export_plot_dpi"),
                                "Plot resolution (dpi)",
                                value = 300,
                                min = 72,
                                step = 25
                            )
                        )
                    ),
                    shiny::h4("Calibration Curve(s)"),
                    shiny::plotOutput(ns("curve_plot"), height = "420px"),
                    shiny::hr(),
                    shiny::h4("Response factor plot(s)"),
                    shiny::plotOutput(ns("rf_plot"), height = "420px"),
                    shiny::hr(),
                    shiny::h4("Summary table"),
                    shiny::tableOutput(ns("summary_tab")),
                    shiny::hr(),
                    shiny::h4("Detailed result tables"),
                    shiny::uiOutput(ns("detail_tables"))
                )
            )
        )
    )
}

#' Calibration tab module server
#'
#' @keywords internal
mod_calibration_server <- function(id, data_input) {
    shiny::moduleServer(id, function(input, output, session) {
        calibration_input_count <- shiny::reactive({
            dat <- data_input()
            if (is.null(dat)) {
                return(0L)
            }
            if (is.data.frame(dat)) {
                return(1L)
            }
            length(dat)
        })

        calibration_result <- shiny::eventReactive(input$run_calibration, {
            shiny::req(data_input())

            calibration_input <- data_input()
            if (is.data.frame(calibration_input)) {
                substance_name <- "substance"
            } else {
                substance_name <- names(calibration_input)[1]
            }

            plot_type_value <- if (!is.null(input$plot_type) && nzchar(input$plot_type)) {
                input$plot_type
            } else if (calibration_input_count() > 1L) {
                "single_plots"
            } else {
                "single_plots"
            }

            CalibraCurve(
                D_list = calibration_input,
                substance = substance_name,
                minReplicates = input$minReplicates,
                cvThres = input$cvThres,
                calcContinuousPrelimRanges = isTRUE(input$calcContinuousPrelimRanges),
                weightingMethod = input$weightingMethod,
                centralTendencyMeasure = input$centralTendencyMeasure,
                perBiasThres = input$perBiasThres,
                considerPerBiasCV = input$considerPerBiasCV,
                perBiasDistThres = input$perBiasDistThres,
                RfThresL = input$RfThresL,
                RfThresU = input$RfThresU,
                RF_colour_threshold = if (!is.null(input$RF_colour_threshold) && nzchar(input$RF_colour_threshold)) input$RF_colour_threshold else "orange",
                RF_colour_within = if (!is.null(input$RF_colour_within) && nzchar(input$RF_colour_within)) input$RF_colour_within else "#00BFC4",
                RF_colour_outside = if (!is.null(input$RF_colour_outside) && nzchar(input$RF_colour_outside)) input$RF_colour_outside else "#F8766D",
                plot_type = plot_type_value,
                point_colour = if (!is.null(input$point_colour)) input$point_colour else "black",
                curve_colour = if (!is.null(input$curve_colour)) input$curve_colour else "red",
                linear_range_colour = if (!is.null(input$linear_range_colour)) input$linear_range_colour else "black",
                multiplot_scales = if (!is.null(input$multiplot_scales) && nzchar(input$multiplot_scales)) input$multiplot_scales else "free",
                base_size = if (!is.null(input$base_size)) input$base_size else 11,
                log10_x_axis = isTRUE(input$log10_x_axis),
                log_tick_label_style = if (!is.null(input$log_tick_label_style) && nzchar(input$log_tick_label_style)) input$log_tick_label_style else "auto",
                show_regression_info = isTRUE(input$show_regression_info),
                show_linear_range = isTRUE(input$show_linear_range),
                show_data_points = isTRUE(input$show_data_points),
                plot_zeros = isTRUE(input$plot_zeros),
                verbose = FALSE,
                output_path = NULL
            )
        })

        output$plot_type_ui <- shiny::renderUI({
            if (calibration_input_count() <= 1L) {
                return(NULL)
            }
            shiny::selectInput(
                session$ns("plot_type"),
                "Calibration plot type",
                choices = c(
                    "single_plots" = "single_plots",
                    "multiplot" = "multiplot",
                    "all_in_one" = "all_in_one"
                ),
                selected = if (is.null(input$plot_type)) "single_plots" else input$plot_type
            )
        })

        shiny::observeEvent(calibration_result(), {
            dataset_names <- .calibration_dataset_names(calibration_result())
            if (length(dataset_names) == 0L) {
                return(NULL)
            }
            selected <- if (is.null(input$selected_dataset) || !(input$selected_dataset %in% dataset_names)) dataset_names[1L] else input$selected_dataset
            shiny::updateSelectInput(session, "selected_dataset", choices = dataset_names, selected = selected)
        }, ignoreNULL = FALSE)

        output$dataset_selector <- shiny::renderUI({
            dataset_names <- .calibration_dataset_names(calibration_result())
            if (length(dataset_names) <= 1L) {
                return(NULL)
            }
            shiny::selectInput(
                session$ns("selected_dataset"),
                "Select dataset / substance",
                choices = dataset_names,
                selected = dataset_names[1L]
            )
        })

        selected_calibration_result <- shiny::reactive({
            shiny::req(calibration_result())
            .select_calibration_dataset(calibration_result(), input$selected_dataset)
        })

        output$summary_tab <- shiny::renderTable({
            shiny::req(selected_calibration_result())
            selected_calibration_result()$summary_tab
        }, rownames = FALSE, bordered = TRUE, hover = TRUE, striped = TRUE)

        output$curve_plot <- shiny::renderPlot({
            shiny::req(selected_calibration_result())
            plot_obj <- .extract_first_ggplot(selected_calibration_result()$plot_CC_list)
            if (is.null(plot_obj)) {
                return(NULL)
            }
            print(plot_obj)
        })

        output$rf_plot <- shiny::renderPlot({
            shiny::req(selected_calibration_result())
            plot_obj <- .extract_first_ggplot(selected_calibration_result()$plot_RF_list)
            if (is.null(plot_obj)) {
                return(NULL)
            }
            print(plot_obj)
        })

        output$detail_tables <- shiny::renderUI({
            shiny::req(selected_calibration_result())
            tables <- .calibration_selection_tables(selected_calibration_result())
            if (length(tables) == 0L) {
                return(shiny::p("No additional tables are available for the selected dataset."))
            }

            shiny::div(
                style = "max-height: 700px; overflow-y: auto;",
                lapply(seq_along(tables), function(i) {
                    nm <- names(tables)[i]
                    shiny::tagList(
                        shiny::h5(nm),
                        shiny::renderTable({
                            tables[[i]]
                        }, rownames = FALSE, bordered = TRUE, hover = TRUE, striped = TRUE)
                    )
                })
            )
        })

        output$download_results <- shiny::downloadHandler(
            filename = function() {
                "calibration_results.zip"
            },
            content = function(file) {
                shiny::req(calibration_result())

                plot_type_for_export <- if (!is.null(input$plot_type) && nzchar(input$plot_type)) {
                    input$plot_type
                } else {
                    "single_plots"
                }

                default_dims <- .default_export_plot_dimensions(plot_type_for_export)
                use_advanced_size <- isTRUE(input$show_advanced_export_settings)
                export_plot_width <- if (use_advanced_size) input$export_plot_width else default_dims$width
                export_plot_height <- if (use_advanced_size) input$export_plot_height else default_dims$height

                .export_calibration_results_zip(
                    calibration_result = calibration_result(),
                    zip_path = file,
                    table_format = input$table_export_format,
                    plot_format = input$plot_export_format,
                    plot_type = plot_type_for_export,
                    plot_width = export_plot_width,
                    plot_height = export_plot_height,
                    plot_dpi = input$export_plot_dpi
                )
            },
            contentType = "application/octet-stream"
        )

        calibration_result
    })
}

#' Prediction tab module UI
#'
#' @keywords internal
mod_prediction_ui <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::textInput(
            ns("new_values"),
            "Intensity values (comma-separated)",
            value = "1000000, 2000000, 3000000"
        ),
        shiny::actionButton(ns("predict_btn"), "Predict concentrations"),
        shiny::tableOutput(ns("prediction_table"))
    )
}

#' Prediction tab module server
#'
#' @keywords internal
mod_prediction_server <- function(id, calibration_result) {
    shiny::moduleServer(id, function(input, output, session) {
        prediction_result <- shiny::eventReactive(input$predict_btn, {
            shiny::req(calibration_result())
            vals <- scan(text = input$new_values, sep = ",", quiet = TRUE)
            if (length(vals) == 0L) {
                return(NULL)
            }

            first_curve <- calibration_result()$RES[[1]]
            if (is.null(first_curve)) {
                return(NULL)
            }

            predictConcentration(
                CC_res = first_curve,
                newdata = vals,
                verbose = TRUE
            )
        })

        output$prediction_table <- shiny::renderTable({
            shiny::req(prediction_result())
            prediction_result()
        }, rownames = FALSE, bordered = TRUE, hover = TRUE, striped = TRUE)

        prediction_result
    })
}

#' Generate the CalibraCurve Shiny app UI
#'
#' @keywords internal
calibra_curve_app_ui <- function() {
    shiny::navbarPage(
        id = "main_tabs",
        title = "CalibraCurve",
        shiny::tabPanel(
            "1. Data upload",
            shiny::fluidPage(
                mod_data_upload_ui("upload")
            )
        ),
        shiny::tabPanel(
            "2. Calibration curve",
            shiny::fluidPage(
                mod_calibration_ui("calibration")
            )
        ),
        shiny::tabPanel(
            "3. Prediction",
            shiny::fluidPage(
                mod_prediction_ui("prediction")
            )
        )
    )
}

#' Generate the server for the CalibraCurve Shiny app
#'
#' @keywords internal
calibra_curve_app_server <- function(input, output, session) {
    upload_state <- mod_data_upload_server("upload")
    upload_data <- upload_state$data

    shiny::observeEvent(upload_state$next_clicked(), {
        if (is.null(upload_data())) {
            return(NULL)
        }
        shiny::updateNavbarPage(session = session, inputId = "main_tabs", selected = "2. Calibration curve")
    })

    calibration_results <- mod_calibration_server(
        id = "calibration",
        data_input = upload_data
    )
    mod_prediction_server("prediction", calibration_result = calibration_results)
}

#' Launch the CalibraCurve Shiny app
#'
#' This starts a small app that wraps the package's public
#' \code{\link{CalibraCurve}} workflow with a 3-tab UI.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' CalibraCurveApp()
#' }
CalibraCurveApp <- function() {
    shiny::shinyApp(
        ui = calibra_curve_app_ui(),
        server = calibra_curve_app_server
    )
}

#' Run the CalibraCurve Shiny app locally
#' @param ... further parameters
#' @export
runCalibraCurveApp <- function(...) {
    shiny::runApp(CalibraCurveApp(), ...)
}

#' #' Alias for the app launcher
#' #'
#' #' @export
#' run_app <- runCalibraCurveApp

