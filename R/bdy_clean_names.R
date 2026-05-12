#' Title
#'
#' @param any_words
#'
#' @returns
#' @export
#'
#' @examples
#'
bdy_clean_names <- function(any_words){
  any_words %>%
    #iconv(., to="ASCII//TRANSLIT") %>%
    tolower %>%
    stringr::str_replace_all(., " ", "_")  %>%
    stringr::str_replace_all(., "'", "_") %>%
    stringr::str_replace_all(., "-", "_") %>%
    stringr::str_replace_all(., ":", "") %>%
    stringr::str_replace_all(., "/", "") %>%
    stringr::str_replace_all(., "«", "") %>%
    stringr::str_replace_all(., "»", "") %>%
    stringr::str_replace_all(., "\\.", "_") %>%
    stringr::str_replace_all(., "___", "_") %>%
    stringr::str_replace_all(., "__", "_")
}
