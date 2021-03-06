# Added '...' to some base functions. These will later be
# turned into default functions by setMethodS3().

# USED TO DO: colSums <- appendVarArgs(colSums)
colSums <- function(...) UseMethod("colSums")
setMethodS3("colSums", "default", function(...) {
  base::colSums(...)
})

# USED TO DO: colMeans <- appendVarArgs(colMeans)
colMeans <- function(...) UseMethod("colMeans")
setMethodS3("colMeans", "default", function(...) {
  base::colMeans(...)
})

write <- appendVarArgs(write)
