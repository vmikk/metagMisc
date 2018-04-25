
## Function to check the availability of dependencies within functions
package_availability <- function(pkg = "plyr"){
  # pkg = package name

  if(! requireNamespace(pkg, quietly = TRUE)){
    stop("Package ", pkg, " is not available. Please install it with:  install.packages('", pkg, "')", sep="")
  }

}
