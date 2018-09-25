
## Add leading zeroes to numeric vector
leading_zero <- function(x, z=3){
    # x = numeric vector
    # z = resulting length of number
    zz <- paste("%0", z, "d", sep="")
    sprintf(zz, x)
}
# leading_zero(1:10, z=3)
