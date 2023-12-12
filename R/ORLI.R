ORLI = function(y1, n1, y2, n2, conf.level=0.95, k, eps=1e-8) # Consider doing in log scale
{
  if (length(y1) > 1) stop("This does not support multiple strata!")
  if (any(c(y1, n1 - y1, y2, n2 - y2) < 0) | n1*n2 == 0) stop("Check the input!")
  p1 = y1/n1
  p2 = y2/n2
  o1 = y1/(n1 - y1)                 # odd of test (active) group
  o2 = y2/(n2 - y2)                 # odd of control (placebo) group
  OR0 = y1/y2*(n2 - y2)/(n1 - y1)   # point estimate of odds ratio (OR)

  if (y1 == 0 | y1 == n1 | y2 == 0 | y2 == n2) {
    maxLL = 0
  } else {
    maxLL = y1*log(p1) + (n1 - y1)*log(1 - p1) + y2*log(p2) + (n2 - y2)*log(1 - p2) # lchoose part removed
  }

  n = n1 + n2
  if (!missing(k)) {
    logk = log(k)
  } else if (n == 1) {
    logk = log(2/(1 - conf.level))
  } else {
    logk = n/2*log(1 + qf(conf.level, 1, n - 1)/(n - 1)) # - 1 regardless of number of parameters
    logk = min(logk, log(2/(1 - conf.level))) # Pawitan p240 k = 20 -> p < 0.05
  }

  Obj = function(or) {  # find or points of increased obj fx value (ofv) by v0
    A = n2*(or - 1)                         # eq 13
    B = n1*or + n2 - (y1 + y2)*(or - 1)
    C1 = -(y1 + y2)
    p2t = (-B + sqrt(B*B - 4*A*C1))/(2*A)
    p1t = p2t*or/(1 + p2t*(or - 1))

    if (p1t < eps | p1t > 1 - eps | p2t < eps | p2t > 1 - eps) {
      ll = 0
    } else {
      ll = y1*log(p1t) + (n1 - y1)*log(1 - p1t) + y2*log(p2t) + (n2 - y2)*log(1 - p2t) # lchoose part removed
    }
    return(maxLL - ll - logk)
  }

  if (y1 == 0 & y2 == 0) {           # Case 1; o1=0 & o2=0
    LL = 0
    UL = Inf
  } else if (y1 == 0 & y2 < n2) {    # Case 2: o1=0 & 0 < o2 < Inf
    LL = 0
    rTemp = try(uniroot(Obj, c(eps, 1e9)), silent=T)
    if (!inherits(rTemp, "try-error")) { UL = rTemp$root
    } else { UL = Inf }
  } else if (y1 == 0 & y2 == n2) {   # Case 3: o1=0 & o2=Inf
    LL = 0
    rTemp = try(uniroot(Obj, c(eps, 1e9)), silent=T)
    if (!inherits(rTemp, "try-error")) { UL = rTemp$root
    } else { UL = Inf }
  } else if (y1 < n1 & y2 == 0) {    # Case 4: 0 < o1 < Inf & o2=0
    rTemp = try(uniroot(Obj, c(eps, 1e9)), silent=T)
    if (!inherits(rTemp, "try-error")) { LL = rTemp$root
    } else { LL = 0 }
    UL = Inf
  } else if (y1 < n1 & y2 < n2) {    # Case 5: 0 < o1 < Inf & 0 < o2 < Inf
#    LE = 3*sqrt(1/y1 + 1/(n1 - y1) + 1/y2 + 1/(n2 - y2))
    rTemp = try(uniroot(Obj, c(1e-4, OR0)), silent=T)
    if (!inherits(rTemp, "try-error")) { LL = rTemp$root
    } else { LL = 0 }
    rTemp = try(uniroot(Obj, c(OR0, 1e4)), silent=T)
    if (!inherits(rTemp, "try-error")) { UL = rTemp$root
    } else { UL = Inf }
  } else if (y1 < n1 & y2 == n2) {   # Case 6: 0 < o1 < 1 & o2=Inf
    LL = 0
    rTemp = try(uniroot(Obj, c(eps, 1e9)), silent=T)
    if (!inherits(rTemp, "try-error")) { UL = rTemp$root
    } else { UL = Inf }
  } else if (y1 == n1 & y2 == 0) {   # Case 7: o1=Inf & o2=0
    rTemp = try(uniroot(Obj, c(eps, 1e9)), silent=T)
    if (!inherits(rTemp, "try-error")) { LL = rTemp$root
    } else { LL = 0 }
    UL = Inf
  } else if (y1 == n1 & y2 < n2) {   # Case 8: o1=Inf & 0 < o2 < Inf
    rTemp = try(uniroot(Obj, c(eps, 1e9)), silent=T)
    if (!inherits(rTemp, "try-error")) { LL = rTemp$root
    } else { LL = 0 }
    UL = Inf
  } else if (y1 == n1 & y2 == n2) {  # Case 9: o1=Inf & o2=Inf
    LL = 0
    UL = Inf
  }

  return(c(odd1 = o1, odd2 = o2, OR = OR0, lower = LL, upper = UL, k = exp(logk)))
}
