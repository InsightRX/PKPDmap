#' Replace list elements by name
#' 
#' @param lst original list
#' @param replacement replacement list
#' 
#' @details
#' Finds and replaces list elements by name and throws an error if an 
#'    element is not available in the original list. This is a local duplicate
#'    of the PKPDmisc copy for the VPC package to reduce dependency on PKPDmisc
#'    at this time.
#'    
#' @examples 
#' \dontrun{
#' lst <- list(ipred = "ipred", dv = "dv", idv = "idv", "pred" = "pred")
#' replacement <- list(dv = "conc", idv = "time")
#' lst <- replace_list_elements(lst, replacement)
#' }
#' 
#' @export
#' 
replace_list_elements <- function(lst, replacement) {
  missing <- which(!names(replacement) %in% names(lst))
  if(length(missing) != 0) {
    warning(paste("Nothing named: ", paste(names(replacement)[missing], collapse= ", ", "found to replace") ))
    replacement <- replacement[-missing]
  }
  lst[names(replacement)] <- lapply(names(replacement), function(x) lst[[x]] <- replacement[[x]])
  lst
}