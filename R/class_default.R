#' Default method seeker_gen_pathway
#'
#' This method return a error message when the function can use the class in the first parameter
#'
#' @param x
#'
#' @return
#' A message of error
#' @export
#'
#' @examples
#' This method return an error message
seeker_gen_pathway.default <- function(x) {
  stop(
    "Don't know how to make seeker_gen_pathway <",
    class(x)[[1]], ">",
    call. = FALSE
  )
}
