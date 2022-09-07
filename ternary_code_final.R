#install.packages('Ternary')
library('Ternary')
par(mfrow = c(1, 1), mar = rep(0.2, 4))

# Figure 1: n = 4 role models and m = 3 variants
D_prime = -0.3
D = 0.9
arrow_len = 5
separation = 0.04

TernaryPlot(alab = expression('p'[1]), blab = expression('p'[2]), clab = expression('p'[3]), axis.labels = seq(0, 1, by = 0.1), grid.lines = 10, grid.lty = 'dotted', grid.minor.lines = 1, grid.minor.lty = 'dotted',)

unstable_equil <- list(
  c(1, 0, 0),
  c(0, 1, 0),
  c(0, 0, 1),
  c(1/3, 1/3, 1/3),
  c(2*D_prime/(2*D_prime-3*D), 2*D_prime/(2*D_prime-3*D), (-2*D_prime-3*D)/(2*D_prime-3*D) ),
  c(2*D_prime/(2*D_prime-3*D), (-2*D_prime-3*D)/(2*D_prime-3*D), 2*D_prime/(2*D_prime-3*D) ),
  c((-2*D_prime-3*D)/(2*D_prime-3*D), 2*D_prime/(2*D_prime-3*D), 2*D_prime/(2*D_prime-3*D) )
)
stable_equil <- list(
  c(1/2, 1/2, 0),
  c(1/2, 0, 1/2),
  c(0, 1/2, 1/2)
)
AddToTernary(points, unstable_equil, pch=1)
AddToTernary(points, stable_equil, pch=19)

next_gen <- function(vec){
  a = vec[[1]]
  b = vec[[2]]
  c = vec[[3]]
  new_a = a + D_prime*a*(1-a)*(2*a-1) + 3*a*(b^2*c + b*c^2)*(D_prime - (1/2)*D) + 3*D*a^2*b*c 
  new_b = b + D_prime*b*(1-b)*(2*b-1) + 3*b*(a^2*c + a*c^2)*(D_prime - (1/2)*D) + 3*D*b^2*a*c 
  new_c = 1-new_a-new_b
  return(c(new_a,new_b,new_c))
}

draw_arrow <- function(start, arrow_len){
  vec = start 
  for(i in 1:arrow_len){
    vec = next_gen(vec)
  }
  if(identical(start,vec) == FALSE){
    TernaryArrows(start, vec, length = 0.04, col = 'dodgerblue3')
  }
}

for(j in seq(from=0, to=1, by=separation)){
  for(k in seq(from=0, to=1, by=separation)){
    if(j+k<=1){
      draw_arrow(c(j,k,round((1-j-k),5)), arrow_len)
    }
  }
}

# Figure 2: n = 5 role models and m = 3 variants
D3 = -1.5 # D triple prime
D2 = 0.9 # D double prime
D1 = 1.9 # D prime
D = -0.5 
arrow_len = 20
separation = 0.05

under_sq_rt =  -3*(D3^2) -4*D3*(D2-D1+3*D) + 4*(D2+D1)^2
numr = sqrt(under_sq_rt) - 2*D1 + 3*D
denmr = 3*D3 - 2*D2 - 4*D1 + 3*D
p = numr / denmr

TernaryPlot(alab = expression('p'[1]), blab = expression('p'[2]), clab = expression('p'[3]), axis.labels = seq(0, 1, by = 0.1), grid.lines = 10, grid.lty = 'dotted', grid.minor.lines = 1, grid.minor.lty = 'dotted',)

unstable_equil <- list(
  c(1, 0, 0),
  c(0, 1, 0),
  c(0, 0, 1),
  c(1/3, 1/3, 1/3),
  c(p, (1-p)/2, (1-p)/2),
  c((1-p)/2, p, (1-p)/2),
  c((1-p)/2, (1-p)/2, p),
  c(1/2, 1/2, 0),
  c(1/2, 0, 1/2),
  c(0, 1/2, 1/2)
)

numr = -1*sqrt(under_sq_rt) - 2*D1 + 3*D
p = numr / denmr
stable_equil <- list(
  c((1-p)/2, (1-p)/2, p),
  c((1-p)/2, p, (1-p)/2),
  c(p, (1-p)/2, (1-p)/2)
)
AddToTernary(points, unstable_equil, pch=1)
AddToTernary(points, stable_equil, pch=19)

next_gen <- function(vec){
  a = vec[[1]]
  b = vec[[2]]
  c = vec[[3]]
  new_a = a + D3*(a*(1-a)*(2*a-1)*(a^2 - a + 1) + a*(4*b^3*c + 6*b^2*c^2 + 4*b*c^3)) +  2*D2*a*(2*a^2*b*c - b^3*c - b*c^3) + 2*D1*(a^3*(b^2+c^2) -a^2*(b^3 + c^3)) + 3*D*(a^2*(b^2*c + b*c^2) - 2*a*b^2*c^2)
  new_b = b + D3*(b*(1-b)*(2*b-1)*(b^2 - b + 1) + b*(4*a^3*c + 6*a^2*c^2 + 4*a*c^3))+  2*D2*b*(2*b^2*a*c - a^3*c - a*c^3) + 2*D1*(b^3*(a^2+c^2) -b^2*(a^3 + c^3) ) +3*D*(b^2*(a^2*c + a*c^2) - 2*b*a^2*c^2)
  new_c = 1-new_a-new_b
  return(c(new_a,new_b,new_c))
}

for(j in seq(from=0, to=1, by=separation)){
  for(k in seq(from=0, to=1, by=separation)){
    if(j+k<=1){
      
      draw_arrow(c(j,k,round((1-j-k),5)), arrow_len)
    }
  }
}

