valid_url <- function(url_in, t = 300){

  # Function to check if a given URL exists or not

  # Flag      Format       Description
  # url_in   <character>   Character vector of the URL to download
  # t        <numeric>     Maximum time until timeout reached

  con <- url(url_in)
  check <- suppressWarnings(try(open.connection(con, open = "rt", timeout = t), silent = T)[1])
  suppressWarnings(try(close.connection(con), silent = T))
  ifelse(is.null(check), TRUE, FALSE)
}