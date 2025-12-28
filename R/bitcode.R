
# encode combinations of N and M classes into integers

# (a, b) = a*2^n + b
bitcode <- function(a,b){ a*(2^3) + b }
bitcode <- Vectorize(bitcode)

## Example
bitcode(a = rep(1:3, each = 3), b = rep(1:3, time = 3))

