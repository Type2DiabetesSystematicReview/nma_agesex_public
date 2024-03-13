PasteAnd <- function(x) {
  if(length(x) == 1) return(x)
  x1 <- x[-length(x)]
  x2 <- x[length(x)]
  x1 <- paste(x1, collapse = ", ")
  paste0(x1, " and ", x2)
}
