#' Loglikelihood functions for the different dose response models
#'
#' @param x parameter values
#' @param nvec vector containing the number of observations at each dose level
#' @param dvec unique ordered dose levels
#' @param mvec vector containing the response at each dose level
#' @param s2vec variance
#' @param qval BMR
#' @param shift value of the shift for negative geometric means
#'
#' @examples
#'
#' @return .
#'
#' @export
#'
llfE4_NI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.E4_NI(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
#' @export
llfIE4_NI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.IE4_NI(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
#' @export
llfH4_NI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.H4_NI(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
#' @export
llfLN4_NI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.LN4_NI(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
#' @export
llfG4_NI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.G4_NI(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
#' @export
llfQE4_NI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.QE4_NI(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
#' @export
llfP4_NI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.P4_NI(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
#' @export
llfSM_N=function(x,nvec,dvec,mvec,s2vec){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[length(dvec)+1]-0.5*(nvec-1)*s2vec*exp(x[length(dvec)+1])-0.5*nvec*((mvec-x[1:length(dvec)])^2)*exp(x[length(dvec)+1]))
}
#' @rdname llfE4_NI
#' @export
llfL4_NI=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.L4_NI(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
#' @export
llfE4_LNI=function(x,nvec,dvec,mvec,s2vec,qval,shift){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*(((mvec+shift)-DRM.E4_LNI(x[1:4],dvec,qval,shift))^2)*exp(x[5])) - sum((mvec+shift)*nvec)
}
#' @rdname llfE4_NI
#' @export
llfIE4_LNI=function(x,nvec,dvec,mvec,s2vec,qval,shift){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*(((mvec+shift)-DRM.IE4_LNI(x[1:4],dvec,qval,shift))^2)*exp(x[5])) - sum((mvec+shift)*nvec)
}
#' @rdname llfE4_NI
#' @export
llfH4_LNI=function(x,nvec,dvec,mvec,s2vec,qval,shift){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*(((mvec+shift)-DRM.H4_LNI(x[1:4],dvec,qval,shift))^2)*exp(x[5])) - sum((mvec+shift)*nvec)
}
#' @rdname llfE4_NI
#' @export
llfLN4_LNI=function(x,nvec,dvec,mvec,s2vec,qval,shift){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*(((mvec+shift)-DRM.LN4_LNI(x[1:4],dvec,qval,shift))^2)*exp(x[5])) - sum((mvec+shift)*nvec)
}
#' @rdname llfE4_NI
#' @export
llfG4_LNI=function(x,nvec,dvec,mvec,s2vec,qval,shift){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*(((mvec+shift)-DRM.G4_LNI(x[1:4],dvec,qval,shift))^2)*exp(x[5])) - sum((mvec+shift)*nvec)
}
#' @rdname llfE4_NI
#' @export
llfQE4_LNI=function(x,nvec,dvec,mvec,s2vec,qval,shift){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*(((mvec+shift)-DRM.QE4_LNI(x[1:4],dvec,qval,shift))^2)*exp(x[5])) - sum((mvec+shift)*nvec)
}
#' @rdname llfE4_NI
#' @export
llfP4_LNI=function(x,nvec,dvec,mvec,s2vec,qval,shift){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*(((mvec+shift)-DRM.P4_LNI(x[1:4],dvec,qval,shift))^2)*exp(x[5])) - sum((mvec+shift)*nvec)
}
#' @rdname llfE4_NI
#' @export
llfL4_LNI=function(x,nvec,dvec,mvec,s2vec,qval,shift){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*(((mvec+shift)-DRM.L4_LNI(x[1:4],dvec,qval,shift))^2)*exp(x[5])) - sum((mvec+shift)*nvec)
}
#' @rdname llfE4_NI
#' @export
llfSM_LN=function(x,nvec,dvec,mvec,s2vec,shift){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[length(dvec)+1]-0.5*(nvec-1)*s2vec*exp(x[length(dvec)+1])-0.5*nvec*(((mvec+shift)-(x[1:length(dvec)]+shift))^2)*exp(x[length(dvec)+1])) - sum((mvec+shift)*nvec)
}
#' @rdname llfE4_NI
#' @export
llfE4_ND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.E4_ND(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
#' @export
llfIE4_ND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.IE4_ND(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
#' @export
llfH4_ND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.H4_ND(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
#' @export
llfLN4_ND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.LN4_ND(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
#' @export
llfG4_ND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.G4_ND(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
#' @export
llfQE4_ND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.QE4_ND(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
#' @export
llfP4_ND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.P4_ND(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
#' @export
llfL4_ND=function(x,nvec,dvec,mvec,s2vec,qval){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*((mvec-DRM.L4_ND(x[1:4],dvec,qval))^2)*exp(x[5]))
}
#' @rdname llfE4_NI
#' @export
llfE4_LND=function(x,nvec,dvec,mvec,s2vec,qval,shift){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*(((mvec+shift)-DRM.E4_LND(x[1:4],dvec,qval,shift))^2)*exp(x[5])) - sum((mvec+shift)*nvec)
}
#' @rdname llfE4_NI
#' @export
llfIE4_LND=function(x,nvec,dvec,mvec,s2vec,qval,shift){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*(((mvec+shift)-DRM.IE4_LND(x[1:4],dvec,qval,shift))^2)*exp(x[5])) - sum((mvec+shift)*nvec)
}
#' @rdname llfE4_NI
#' @export
llfH4_LND=function(x,nvec,dvec,mvec,s2vec,qval,shift){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*(((mvec+shift)-DRM.H4_LND(x[1:4],dvec,qval,shift))^2)*exp(x[5])) - sum((mvec+shift)*nvec)
}
#' @rdname llfE4_NI
#' @export
llfLN4_LND=function(x,nvec,dvec,mvec,s2vec,qval,shift){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*(((mvec+shift)-DRM.LN4_LND(x[1:4],dvec,qval,shift))^2)*exp(x[5])) - sum((mvec+shift)*nvec)
}
#' @rdname llfE4_NI
#' @export
llfG4_LND=function(x,nvec,dvec,mvec,s2vec,qval,shift){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*(((mvec+shift)-DRM.G4_LND(x[1:4],dvec,qval,shift))^2)*exp(x[5])) - sum((mvec+shift)*nvec)
}
#' @rdname llfE4_NI
#' @export
llfQE4_LND=function(x,nvec,dvec,mvec,s2vec,qval,shift){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*(((mvec+shift)-DRM.QE4_LND(x[1:4],dvec,qval,shift))^2)*exp(x[5])) - sum((mvec+shift)*nvec)
}
#' @rdname llfE4_NI
#' @export
llfP4_LND=function(x,nvec,dvec,mvec,s2vec,qval,shift){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*(((mvec+shift)-DRM.P4_LND(x[1:4],dvec,qval,shift))^2)*exp(x[5])) - sum((mvec+shift)*nvec)
}
#' @rdname llfE4_NI
#' @export
llfL4_LND=function(x,nvec,dvec,mvec,s2vec,qval,shift){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*(((mvec+shift)-DRM.L4_LND(x[1:4],dvec,qval,shift))^2)*exp(x[5])) - sum((mvec+shift)*nvec)
}

#' Loglikelihood functions for the different dose response models - clustered
#'
#' @param x parameter values
#' @param d dose levels
#' @param n number of litters per dose level
#' @param nij dose x litter matrix with the number of observations for each combination
#' @param y response matrix
#' @param qval BMR
#' @param shift value of the shift for negative geometric means
#'
#' @examples
#'
#' @return .
#'
#' @export
#'
llfE4_NIc = function(x, d, n, nij, y, qval){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt]
      mx = DRM.E4_NI(x[1:4], d[i], qval)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfIE4_NIc = function(x, d, n, nij, y, qval){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt]
      mx = DRM.IE4_NI(x[1:4], d[i], qval)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T)

      # ym = rsp - m
      # ymt = t(ym)
      #
      # sum(-0.5*lt*log(2*pi) +
      #       0.5*lt*log(x[5]) -
      #       0.5*determinant(P, logarithm = T)$modulus[1] -
      #       0.5*(ym%*%ymt * solve(P))*x[5])

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfH4_NIc = function(x, d, n, nij, y, qval){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt]
      mx = DRM.H4_NI(x[1:4], d[i], qval)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfLN4_NIc = function(x, d, n, nij, y, qval){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt]
      mx = DRM.LN4_NI(x[1:4], d[i], qval)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfG4_NIc = function(x, d, n, nij, y, qval){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt]
      mx = DRM.G4_NI(x[1:4], d[i], qval)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfQE4_NIc = function(x, d, n, nij, y, qval){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt]
      mx = DRM.QE4_NI(x[1:4], d[i], qval)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfP4_NIc = function(x, d, n, nij, y, qval){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt]
      mx = DRM.P4_NI(x[1:4], d[i], qval)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfL4_NIc = function(x, d, n, nij, y, qval){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt]
      mx = DRM.L4_NI(x[1:4], d[i], qval)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfSM_N=function(x,nvec,dvec,mvec,s2vec){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[length(dvec)+1]-0.5*(nvec-1)*s2vec*exp(x[length(dvec)+1])-0.5*nvec*((mvec-x[1:length(dvec)])^2)*exp(x[length(dvec)+1]))
}
#' @rdname llfE4_NI
#' @export
llfE4_LNIc = function(x, d, n, nij, y, qval, shift){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt] + shift
      mx = DRM.E4_LNI(x[1:4], d[i], qval, shift)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T) - sum(rsp)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfIE4_LNIc = function(x, d, n, nij, y, qval, shift){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt] + shift
      mx = DRM.IE4_LNI(x[1:4], d[i], qval, shift)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T) - sum(rsp)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfH4_LNIc = function(x, d, n, nij, y, qval, shift){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt] + shift
      mx = DRM.H4_LNI(x[1:4], d[i], qval, shift)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T) - sum(rsp)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfLN4_LNIc = function(x, d, n, nij, y, qval, shift){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt] + shift
      mx = DRM.LN4_LNI(x[1:4], d[i], qval, shift)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T) - sum(rsp)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfG4_LNIc = function(x, d, n, nij, y, qval, shift){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt] + shift
      mx = DRM.G4_LNI(x[1:4], d[i], qval, shift)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T) - sum(rsp)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfQE4_LNIc = function(x, d, n, nij, y, qval, shift){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt] + shift
      mx = DRM.QE4_LNI(x[1:4], d[i], qval, shift)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T) - sum(rsp)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfP4_LNIc = function(x, d, n, nij, y, qval, shift){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt] + shift
      mx = DRM.P4_LNI(x[1:4], d[i], qval, shift)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T) - sum(rsp)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfL4_LNIc = function(x, d, n, nij, y, qval, shift){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt] + shift
      mx = DRM.L4_LNI(x[1:4], d[i], qval, shift)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T) - sum(rsp)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfSM_LNc = function(x, d, n, nij, y, qval, shift){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt] + shift
      mx = x[i]
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[length(d)+2]
          }
        }
      }
      Sigma = 1/exp(x[length(d)+1]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T) - sum(rsp)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfSM_Nc = function(x, d, n, nij, y, qval){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt]
      mx = x[i]
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[length(d)+2]
          }
        }
      }
      Sigma = 1/exp(x[length(d)+1]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfE4_NDc = function(x, d, n, nij, y, qval){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt]
      mx = DRM.E4_ND(x[1:4], d[i], qval)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T)

      # ym = rsp - m
      # ymt = t(ym)
      #
      # sum(-0.5*lt*log(2*pi) +
      #       0.5*lt*log(x[5]) -
      #       0.5*determinant(P, logarithm = T)$modulus[1] -
      #       0.5*(ym%*%ymt * solve(P))*x[5])

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfIE4_NDc = function(x, d, n, nij, y, qval){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt]
      mx = DRM.IE4_ND(x[1:4], d[i], qval)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfH4_NDc = function(x, d, n, nij, y, qval){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt]
      mx = DRM.H4_ND(x[1:4], d[i], qval)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfLN4_NDc = function(x, d, n, nij, y, qval){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt]
      mx = DRM.LN4_ND(x[1:4], d[i], qval)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfG4_NDc = function(x, d, n, nij, y, qval){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt]
      mx = DRM.G4_ND(x[1:4], d[i], qval)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfQE4_NDc = function(x, d, n, nij, y, qval){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt]
      mx = DRM.QE4_ND(x[1:4], d[i], qval)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfP4_NDc = function(x, d, n, nij, y, qval){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt]
      mx = DRM.P4_ND(x[1:4], d[i], qval)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfL4_NDc = function(x, d, n, nij, y, qval){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt]
      mx = DRM.L4_ND(x[1:4], d[i], qval)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, m, Sigma, log = T)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfE4_LNDc = function(x, d, n, nij, y, qval, shift){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){ ## ll per cluster
      lt = nij[i, j]
      rsp = y[cnt, 1:lt] + shift
      mx = DRM.E4_LND(x[1:4], d[i], qval, shift)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, mean = m, sigma = Sigma, log = T) - sum(rsp)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfIE4_LNDc = function(x, d, n, nij, y, qval, shift){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt] + shift
      mx = DRM.IE4_LND(x[1:4], d[i], qval, shift)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, mean = m, sigma = Sigma, log = T) - sum(rsp)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfH4_LNDc = function(x, d, n, nij, y, qval, shift){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt] + shift
      mx = DRM.H4_LND(x[1:4], d[i], qval, shift)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, mean = m, sigma = Sigma, log = T) - sum(rsp)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfLN4_LNDc = function(x, d, n, nij, y, qval, shift){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt] + shift
      mx = DRM.LN4_LND(x[1:4], d[i], qval, shift)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, mean = m, sigma = Sigma, log = T) - sum(rsp)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfG4_LNDc = function(x, d, n, nij, y, qval, shift){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt] + shift
      mx = DRM.G4_LND(x[1:4], d[i], qval, shift)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, mean = m, sigma = Sigma, log = T) - sum(rsp)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfQE4_LNDc = function(x, d, n, nij, y, qval, shift){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt] + shift
      mx = DRM.QE4_LND(x[1:4], d[i], qval, shift)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, mean = m, sigma = Sigma, log = T) - sum(rsp)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfP4_LNDc = function(x, d, n, nij, y, qval, shift){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt] + shift
      mx = DRM.P4_LND(x[1:4], d[i], qval, shift)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, mean = m, sigma = Sigma, log = T) - sum(rsp)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}
#' @rdname llfE4_NIc
#' @export
llfL4_LNDc = function(x, d, n, nij, y, qval, shift){
  ll = c()
  cnt = 1
  for(i in 1:length(d)){
    for(j in 1:n[i]){
      lt = nij[i, j]
      rsp = y[cnt, 1:lt] + shift
      mx = DRM.L4_LND(x[1:4], d[i], qval, shift)
      m = rep(mx, lt)
      P = matrix(0, nrow = lt, ncol = lt)
      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1
          }else{
            P[id1, id2] = x[6]
          }
        }
      }
      Sigma = 1/exp(x[5]) * P

      llij = mvtnorm::dmvnorm(rsp, mean = m, sigma = Sigma, log = T) - sum(rsp)

      cnt = cnt + 1

      ll = c(ll, llij)
    }
  }
  return(sum(ll))
}



#' Loglikelihood functions for the different dose response models
#'
#' @param x parameter values
#' @param nvec number of observations
#' @param dvec dose levels
#' @param yvec response
#' @param qval BMR
#'
#' @return loglikelihood value at x.
#'
#' @export
#'
### Binomial
llfE4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.E4_Q(x[1:3], dvec, qval) + .Machine$double.xmin) +
        (nvec - yvec)*log(1 - DRM.E4_Q(x, dvec, qval) + .Machine$double.xmin))
}
#' @rdname llfE4_Q
#' @export
llfIE4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.IE4_Q(x[1:3], dvec, qval) + .Machine$double.xmin) +
        (nvec - yvec)*log(1 - DRM.IE4_Q(x, dvec, qval) + .Machine$double.xmin))
}
#' @rdname llfE4_Q
#' @export
llfH4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.H4_Q(x[1:3], dvec, qval) + .Machine$double.xmin) +
        (nvec - yvec)*log(1 - DRM.H4_Q(x, dvec, qval) + .Machine$double.xmin))
}
#' @rdname llfE4_Q
#' @export
llfLN4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.LN4_Q(x[1:3], dvec, qval) + .Machine$double.xmin) +
        (nvec - yvec)*log(1 - DRM.LN4_Q(x, dvec, qval) + .Machine$double.xmin))
}
#' @rdname llfE4_Q
#' @export
llfG4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.G4_Q(x[1:3], dvec, qval) + .Machine$double.xmin) +
        (nvec - yvec)*log(1 - DRM.G4_Q(x, dvec, qval) + .Machine$double.xmin))
}
#' @rdname llfE4_Q
#' @export
llfQE4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.QE4_Q(x[1:3], dvec, qval) + .Machine$double.xmin) +
        (nvec - yvec)*log(1 - DRM.QE4_Q(x, dvec, qval) + .Machine$double.xmin))
}
#' @rdname llfE4_Q
#' @export
llfP4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.P4_Q(x[1:3], dvec, qval) + .Machine$double.xmin) +
        (nvec - yvec)*log(1 - DRM.P4_Q(x, dvec, qval)))
}
#' @rdname llfE4_Q
#' @export
llfL4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.L4_Q(x[1:3], dvec, qval) + .Machine$double.xmin) +
        (nvec - yvec)*log(1 - DRM.L4_Q(x, dvec, qval) + .Machine$double.xmin))
}
#' @rdname llfE4_Q
#' @export
llfSM_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(x[1:length(dvec)] + .Machine$double.xmin) +
        (nvec - yvec)*log(1 - x[1:length(dvec)] + .Machine$double.xmin))
}

#' Loglikelihood functions for the different dose response models (Beta-Binomial)
#'
#' @param x parameter values
#' @param nvec number of observations
#' @param dvec dose levels
#' @param yvec response
#' @param qval BMR
#' @param rho intra-cluster correlation parameter
#'
#' @return loglikelihood value at x.
#'
#' @export
#'
llfE42_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = DRM.E4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+.Machine$double.xmin) + lgamma(bbet+nvec-yvec+.Machine$double.xmin) -
        lgamma(abet+bbet+nvec+.Machine$double.xmin) - lgamma(abet+.Machine$double.xmin) -
        lgamma(bbet+.Machine$double.xmin) +
        lgamma(abet+bbet+.Machine$double.xmin))
}
#' @rdname llfE42_Q
#' @export
llfIE42_Q=function(x,nvec,dvec,yvec,qval,rho){

  m = DRM.IE4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+.Machine$double.xmin) + lgamma(bbet+nvec-yvec+.Machine$double.xmin) -
        lgamma(abet+bbet+nvec+.Machine$double.xmin) - lgamma(abet+.Machine$double.xmin) -
        lgamma(bbet+.Machine$double.xmin) +
        lgamma(abet+bbet+.Machine$double.xmin))
}
#' @rdname llfE42_Q
#' @export
llfH42_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = DRM.H4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+.Machine$double.xmin) + lgamma(bbet+nvec-yvec+.Machine$double.xmin) -
        lgamma(abet+bbet+nvec+.Machine$double.xmin) - lgamma(abet+.Machine$double.xmin) -
        lgamma(bbet+.Machine$double.xmin) +
        lgamma(abet+bbet+.Machine$double.xmin))
}
#' @rdname llfE42_Q
#' @export
llfLN42_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = DRM.LN4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+.Machine$double.xmin) + lgamma(bbet+nvec-yvec+.Machine$double.xmin) -
        lgamma(abet+bbet+nvec+.Machine$double.xmin) - lgamma(abet+.Machine$double.xmin) -
        lgamma(bbet+.Machine$double.xmin) +
        lgamma(abet+bbet+.Machine$double.xmin))
}
#' @rdname llfE42_Q
#' @export
llfG42_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = DRM.G4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+.Machine$double.xmin) + lgamma(bbet+nvec-yvec+.Machine$double.xmin) -
        lgamma(abet+bbet+nvec+.Machine$double.xmin) - lgamma(abet+.Machine$double.xmin) -
        lgamma(bbet+.Machine$double.xmin) +
        lgamma(abet+bbet+.Machine$double.xmin))
}
#' @rdname llfE42_Q
#' @export
llfQE42_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = DRM.QE4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+.Machine$double.xmin) + lgamma(bbet+nvec-yvec+.Machine$double.xmin) -
        lgamma(abet+bbet+nvec+.Machine$double.xmin) - lgamma(abet+.Machine$double.xmin) -
        lgamma(bbet+.Machine$double.xmin) +
        lgamma(abet+bbet+.Machine$double.xmin))
}
#' @rdname llfE42_Q
#' @export
llfP42_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = DRM.P4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+.Machine$double.xmin) + lgamma(bbet+nvec-yvec+.Machine$double.xmin) -
        lgamma(abet+bbet+nvec+.Machine$double.xmin) - lgamma(abet+.Machine$double.xmin) -
        lgamma(bbet+.Machine$double.xmin) +
        lgamma(abet+bbet+.Machine$double.xmin))
}
#' @rdname llfE42_Q
#' @export
llfL42_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = DRM.IE4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+.Machine$double.xmin) + lgamma(bbet+nvec-yvec+.Machine$double.xmin) -
        lgamma(abet+bbet+nvec+.Machine$double.xmin) - lgamma(abet+.Machine$double.xmin) -
        lgamma(bbet+.Machine$double.xmin) +
        lgamma(abet+bbet+.Machine$double.xmin))
}
#' @rdname llfE42_Q
#' @export
llfSM2_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = x[1:length(dvec)]
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+.Machine$double.xmin) + lgamma(bbet+nvec-yvec+.Machine$double.xmin) -
        lgamma(abet+bbet+nvec+.Machine$double.xmin) - lgamma(abet+.Machine$double.xmin) -
        lgamma(bbet+.Machine$double.xmin) +
        lgamma(abet+bbet+.Machine$double.xmin))
}

