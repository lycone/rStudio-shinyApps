#' Starts the UMOD App
#'
#' Starts the UMOD App in shiny.
#' @export


launch_umod_app <- function (){

  appDir <- system.file("shinyApp", package = "umod")

  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `umod`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
