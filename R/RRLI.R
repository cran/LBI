RRLI = function(y1, n1, y2, n2, conf.level=0.95, k, eps=1e-8)
{
  if (length(y1) > 1) stop("This does not support multiple strata!")
  if (any(c(y1, n1 - y1, y2, n2 - y2) < 0) | n1*n2 == 0) stop("Check the input!")
  p1 = y1/n1                 # p of test (active) group
  p2 = y2/n2                 # p of control (placebo) group
  RR0 = p1/p2                # point estimate of relative risk (RR)

  n = n1 + n2
  if (!missing(k)) {
    logk = log(k)
  } else if (n == 1) {
    logk = log(2/(1 - conf.level))
  } else {
    logk = n/2*log(1 + qf(conf.level, 1, n - 1)/(n - 1)) # - 1 regardless of number of parameters
    logk = min(logk, log(2/(1 - conf.level))) # Pawitan p240 k = 20 -> p < 0.05
  }

  LE = 1e4*(exp(logk)/(n1 + n2))
#  logk = min(log(maxk), logk)
  p1 = min(1 - eps, max(eps, p1))
  p2 = min(1 - eps, max(eps, p2))
  if (y1 == 0 | y1 == n1 | y2 == 0 | y2 == n2) {
    maxLL = 0
  } else {
    maxLL = y1*log(p1) + (n1 - y1)*log(1 - p1) + y2*log(p2) + (n2 - y2)*log(1 - p2) # lchoose part removed
  }

  Obj = function(rr) {
    A = (n1 + n2)*rr
    B = n1*rr + y1 + n2 + y2*rr
    C1 = y1 + y2
    p2t = min(1 - eps, max(eps, (B - sqrt(B*B - 4*A*C1))/(2*A)))
    p1t = min(1 - eps, max(eps, p2t*rr))

    ll = y1*log(p1t) + (n1 - y1)*log(1 - p1t) + y2*log(p2t) + (n2 - y2)*log(1 - p2t) # lchoose part removed
    return(maxLL - ll - logk)
  }

  if (y1 == 0 & y2 == 0) {                         # Case 1: p1 = 0 & p2 = 0
    LL = 0
    UL = Inf
  } else if (y1 == 0 & y2 > 0 & y2 < n2) {         # Case 2: p1 = 0 & 0 < p2 < 1
    LL = 0
    rTemp = try(uniroot(Obj, c(eps, 1e9)), silent=T)
    if (!inherits(rTemp, "try-error")) { UL = rTemp$root
    } else { UL = 1e9 }
  } else if (y1 == 0 & y2 == n2) {                 # Case 3: p1 = 0 & p2 = 1
    LL = 0
    rTemp = try(uniroot(Obj, c(eps, 1e9)), silent=T)
    if (!inherits(rTemp, "try-error")) { UL = rTemp$root
    } else { UL = Inf }
  } else if (y1 > 0 & y1 < n1 & y2 == 0) {         # Case 4: 0 < p1 < 1 & p2 = 0
    rTemp = try(uniroot(Obj, c(eps, 1e2)), silent=T)
    if (!inherits(rTemp, "try-error")) { LL = rTemp$root
    } else { LL = 0 }
    UL = Inf
  } else if (y1 > 0 & y1 < n1 & y2 > 0 & y2 < n2) {# Case 5: 0 < p1 < 1 & 0 < p2 < 1
    rTemp = try(uniroot(Obj, c(eps, RR0 + eps)), silent=T)
    if (!inherits(rTemp, "try-error")) { LL = rTemp$root
    } else { LL = 0 }
    rTemp = try(uniroot(Obj, c(RR0 - eps, RR0 + LE)), silent=T)
    if (!inherits(rTemp, "try-error")) { UL = rTemp$root
    } else { UL = Inf }
  } else if (y1 > 0 & y1 < n1 & y2 == n2) {        # Case 6: 0 < p1 < 1 & p2 = 1
    rTemp = try(uniroot(Obj, c(eps, RR0 + eps)), silent=T)
    if (!inherits(rTemp, "try-error")) { LL = rTemp$root
    } else { LL = 0 }
    rTemp = try(uniroot(Obj, c(RR0 - eps, 1e9)), silent=T)
    if (!inherits(rTemp, "try-error")) { UL = rTemp$root
    } else { UL = Inf }
  } else if (y1 == n1 & y2 == 0) {                 # Case 7: p1 = 1 &  p2 = 0
    rTemp = try(uniroot(Obj, c(eps, 1e9)), silent=T)
    if (!inherits(rTemp, "try-error")) { LL = 1/rTemp$root
    } else { LL = 0 }
    UL = Inf
  } else if (y1 == n1 & y2 > 0 & y2 < n2) {        # Case 8: p1 = 1 & 0 < p2 < 1
    rTemp = try(uniroot(Obj, c(eps, RR0 + eps)), silent=T)
    if (!inherits(rTemp, "try-error")) { LL = rTemp$root
    } else { LL = 0 }
    rTemp = try(uniroot(Obj, c(RR0 - eps, 1e9)), silent=T)
    if (!inherits(rTemp, "try-error")) { UL = rTemp$root
    } else { UL = Inf }
  } else if (y1 == n1 & y2 == n2) {                # Case 9: p1 = 1 &  p2 = 1
    rTemp = try(uniroot(Obj, c(eps, RR0 + eps)), silent=T)
    if (!inherits(rTemp, "try-error")) { LL = rTemp$root
    } else { LL = 0 }
    rTemp = try(uniroot(Obj, c(RR0 - eps, min(1 + 1/LL, 1e9))), silent=T)
    if (!inherits(rTemp, "try-error")) { UL = rTemp$root
    } else { UL = Inf }
  }
  return(c(p1 = y1/n1, p2 = y2/n2, RR = RR0, lower = LL, upper = UL, k = exp(logk)))
}
