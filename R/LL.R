
#' Loglikelihood functions for the different dose response models
#'
#' @param x parameter values
#' @param nvec vector containing the number of observations at each dose level
#' @param dvec vector containing the unique ordered dose levels
#' @param mvec vector containing the mean response at each dose level
#' @param s2vec vector containing the variance at each dose level
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
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[5]-0.5*(nvec-1)*s2vec*exp(x[5])-0.5*nvec*
        ((mvec-DRM.E4_NI(x[1:4],dvec,qval))^2)*exp(x[5]))
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


#' Loglikelihood functions for the different dose response models
#'
#' @param x parameter values
#' @param d dose levels
#' @param n observations per dose level
#' @param nij observations per dose x litter combination
#' @param y observed responses
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
#' @rdname llfE4_NI
#' @export
llfSM_N=function(x,nvec,dvec,mvec,s2vec){
  sum(-0.5*nvec*log(2*pi)+0.5*nvec*x[length(dvec)+1]-0.5*(nvec-1)*s2vec*exp(x[length(dvec)+1])-0.5*nvec*((mvec-x[1:length(dvec)])^2)*exp(x[length(dvec)+1]))
}
#' @rdname llfE4_NIc
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
#' @param nvec vector containing the number of observations at each dose level
#' @param dvec vector containing the unique ordered dose levels
#' @param yvec vector containing the number of adverse events at each dose level
#' @param qval BMR
#'
#' @examples
#'
#' @return .
#'
#' @export
#'
### Binomial
llfE4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.E4_Q(x[1:3], dvec, qval) + 1.0E-05) +
        (nvec - yvec)*log(1 - DRM.E4_Q(x, dvec, qval) + 1.0E-05))
}
#' @rdname llfE4_Q
#' @export
llfIE4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.IE4_Q(x[1:3], dvec, qval) + 1.0E-05) +
        (nvec - yvec)*log(1 - DRM.IE4_Q(x, dvec, qval) + 1.0E-05))
}
#' @rdname llfE4_Q
#' @export
llfH4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.H4_Q(x[1:3], dvec, qval) + 1.0E-05) +
        (nvec - yvec)*log(1 - DRM.H4_Q(x, dvec, qval) + 1.0E-05))
}
#' @rdname llfE4_Q
#' @export
llfLN4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.LN4_Q(x[1:3], dvec, qval) + 1.0E-05) +
        (nvec - yvec)*log(1 - DRM.LN4_Q(x, dvec, qval) + 1.0E-05))
}
#' @rdname llfE4_Q
#' @export
llfG4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.G4_Q(x[1:3], dvec, qval) + 1.0E-05) +
        (nvec - yvec)*log(1 - DRM.G4_Q(x, dvec, qval) + 1.0E-05))
}
#' @rdname llfE4_Q
#' @export
llfQE4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.QE4_Q(x[1:3], dvec, qval) + 1.0E-05) +
        (nvec - yvec)*log(1 - DRM.QE4_Q(x, dvec, qval) + 1.0E-05))
}
#' @rdname llfE4_Q
#' @export
llfP4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.P4_Q(x[1:3], dvec, qval) + 1.0E-05) +
        (nvec - yvec)*log(1 - DRM.P4_Q(x, dvec, qval) + 1.0E-05))
}
#' @rdname llfE4_Q
#' @export
llfL4_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(DRM.L4_Q(x[1:3], dvec, qval) + 1.0E-05) +
        (nvec - yvec)*log(1 - DRM.L4_Q(x, dvec, qval) + 1.0E-05))
}
#' @rdname llfE4_Q
#' @export
llfSM_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(x[1:length(dvec)] + 1.0E-05) +
        (nvec - yvec)*log(1 - x[1:length(dvec)] + 1.0E-05))
}
#' @rdname llfE4_Q
#' @export
llfH0_Q=function(x,nvec,dvec,yvec,qval){
  sum(lchoose(nvec, yvec) + yvec*log(x[1] + 1.0E-05) + (nvec-yvec)*log(1 - x[1]+ 1.0E-05))
}


### Beta-Binomial

#' Loglikelihood functions for the different dose response models
#'
#' @param x parameter values
#' @param nvec vector containing the number of observations at each dose level
#' @param dvec vector containing the unique ordered dose levels
#' @param yvec vector containing the number of adverse events at each dose level
#' @param qval BMR
#' @param rho intra-cluster correlation parameter
#'
#' @examples
#'
#' @return .
#'
#' @export
#'
llfE42_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = DRM.E4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+1.0E-05) + lgamma(bbet+nvec-yvec+1.0E-05) -
        lgamma(abet+bbet+nvec+1.0E-05) - lgamma(abet+1.0E-05) -
        lgamma(bbet+1.0E-05) +
        lgamma(abet+bbet+1.0E-05))
}
#' @rdname llfE42_Q
#' @export
llfIE42_Q=function(x,nvec,dvec,yvec,qval,rho){

  m = DRM.IE4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+1.0E-05) + lgamma(bbet+nvec-yvec+1.0E-05) -
        lgamma(abet+bbet+nvec+1.0E-05) - lgamma(abet+1.0E-05) -
        lgamma(bbet+1.0E-05) +
        lgamma(abet+bbet+1.0E-05))
}


#' @rdname llfE42_Q
#' @export
llfH42_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = DRM.H4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+1.0E-05) + lgamma(bbet+nvec-yvec+1.0E-05) -
        lgamma(abet+bbet+nvec+1.0E-05) - lgamma(abet+1.0E-05) -
        lgamma(bbet+1.0E-05) +
        lgamma(abet+bbet+1.0E-05))
}


#' @rdname llfE42_Q
#' @export
llfLN42_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = DRM.LN4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+1.0E-05) + lgamma(bbet+nvec-yvec+1.0E-05) -
        lgamma(abet+bbet+nvec+1.0E-05) - lgamma(abet+1.0E-05) -
        lgamma(bbet+1.0E-05) +
        lgamma(abet+bbet+1.0E-05))
}


#' @rdname llfE42_Q
#' @export
llfG42_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = DRM.G4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+1.0E-05) + lgamma(bbet+nvec-yvec+1.0E-05) -
        lgamma(abet+bbet+nvec+1.0E-05) - lgamma(abet+1.0E-05) -
        lgamma(bbet+1.0E-05) +
        lgamma(abet+bbet+1.0E-05))
}


#' @rdname llfE42_Q
#' @export
llfQE42_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = DRM.QE4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+1.0E-05) + lgamma(bbet+nvec-yvec+1.0E-05) -
        lgamma(abet+bbet+nvec+1.0E-05) - lgamma(abet+1.0E-05) -
        lgamma(bbet+1.0E-05) +
        lgamma(abet+bbet+1.0E-05))
}


#' @rdname llfE42_Q
#' @export
llfP42_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = DRM.P4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+1.0E-05) + lgamma(bbet+nvec-yvec+1.0E-05) -
        lgamma(abet+bbet+nvec+1.0E-05) - lgamma(abet+1.0E-05) -
        lgamma(bbet+1.0E-05) +
        lgamma(abet+bbet+1.0E-05))
}


#' @rdname llfE42_Q
#' @export
llfL42_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = DRM.IE4_Q(x[1:3], dvec, qval)
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+1.0E-05) + lgamma(bbet+nvec-yvec+1.0E-05) -
        lgamma(abet+bbet+nvec+1.0E-05) - lgamma(abet+1.0E-05) -
        lgamma(bbet+1.0E-05) +
        lgamma(abet+bbet+1.0E-05))
}



#' @rdname llfE42_Q
#' @export
llfSM2_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = yvec/nvec #x[1:length(dvec)]
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet+yvec+1.0E-05) + lgamma(bbet+nvec-yvec+1.0E-05) -
        lgamma(abet+bbet+nvec+1.0E-05) - lgamma(abet+1.0E-05) -
        lgamma(bbet+1.0E-05) +
        lgamma(abet+bbet+1.0E-05))
}
#' @rdname llfE42_Q
#' @export
llfH02_Q=function(x,nvec,dvec,yvec,qval,rho){
  m = yvec/nvec
  abet = m*((1/rho)-1)
  bbet = (1.0 - m)*((1/rho)-1)

  sum(lchoose(nvec, yvec) + lgamma(abet + yvec + 1.0E-05) + lgamma(bbet + nvec - yvec+1.0E-05) -
        lgamma(abet + bbet + nvec + 1.0E-05) - lgamma(abet+1.0E-05) -
        lgamma(bbet+1.0E-05) +
        lgamma(abet + bbet + 1.0E-05))
}

#' Loglikelihood functions for the different dose response models
#'
#' @param pars parameter values
#' @param x vector containing the unique ordered dose levels
#' @param n vector containing the number of observations at each dose level
#' @param m vector containing the mean response at each dose level
#' @param s2 vector containing the variance at each dose level
#' @param qval BMR
#' @param shift value of the shift for negative geometric means
#' @param covar which parameter includes a covariate effect
#' @param nlevels number of covariate levels
#' @param trt_ind matrix indicating which reponse corresponds to which covariate level
#'
#' @examples
#'
#' @return .
#'
#' @export
#'
llfE4_NI_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                        nlevels, trt_ind){

  covar = match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    # c <- pars[grep("c\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }

  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m - DRM.E4_NI(c(a[1], bmd[mn], c[1], d[mn]), x, qval))^2)*
                             exp(s[1]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m - DRM.E4_NI(c(a[mn], bmd[1], c[1], d[1]), x, qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m - DRM.E4_NI(c(a[mn], bmd[mn], c[1], d[mn]), x, qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*
                n*((m - DRM.E4_NI(c(a[1], bmd[1], c[1], d[1]),x,qval))^2)*exp(s[1]))

  }


  return(ll)
}

#' @rdname llfE4_NI_Cov
#' @export
llfIE4_NI_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                         nlevels, trt_ind){

  covar <- match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    # c <- pars[grep("c\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }

  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m - DRM.IE4_NI(c(a[1], bmd[mn], c[1], d[mn]), x, qval))^2)*
                             exp(s[1]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m - DRM.IE4_NI(c(a[mn], bmd[1], c[1], d[1]), x, qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m - DRM.IE4_NI(c(a[mn], bmd[mn], c[1], d[mn]), x, qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*n*
                ((m - DRM.IE4_NI(c(a[1], bmd[1], c[1], d[1]),x,qval))^2)*exp(s[1]))
  }


  return(ll)

}

#' @rdname llfE4_NI_Cov
#' @export
llfH4_NI_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                        nlevels, trt_ind){

  covar <- match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    # c <- pars[grep("c\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }

  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m-DRM.H4_NI(c(a[1], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[1]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m-DRM.H4_NI(c(a[mn], bmd[1], c[1], d[1]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m-DRM.H4_NI(c(a[mn], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*n*
                ((m - DRM.H4_NI(c(a[1], bmd[1], c[1], d[1]),x,qval))^2)*exp(s[1]))
  }

  return(ll)

}


#' @rdname llfE4_NI_Cov
#' @export
llfLN4_NI_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                         nlevels, trt_ind){

  covar <- match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    # c <- pars[grep("c\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }

  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m-DRM.LN4_NI(c(a[1], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[1]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m - DRM.LN4_NI(c(a[mn], bmd[1], c[1], d[1]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m-DRM.LN4_NI(c(a[mn], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*n*
                ((m - DRM.LN4_NI(c(a[1], bmd[1], c[1], d[1]),x,qval))^2)*exp(s[1]))
  }


  return(ll)

}

#' @rdname llfE4_NI_Cov
#' @export
llfG4_NI_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                        nlevels, trt_ind){

  covar <- match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    # c <- pars[grep("c\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }

  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m-DRM.G4_NI(c(a[1], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[1]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m - DRM.G4_NI(c(a[mn], bmd[1], c[1], d[1]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m- DRM.G4_NI(c(a[mn], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*n*
                ((m - DRM.G4_NI(c(a[1], bmd[1], c[1], d[1]),x,qval))^2)*exp(s[1]))
  }


  return(ll)

}

#' @rdname llfE4_NI_Cov
#' @export
llfQE4_NI_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                         nlevels, trt_ind){

  covar <- match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    # c <- pars[grep("c\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }

  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m-DRM.QE4_NI(c(a[1], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[1]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m - DRM.QE4_NI(c(a[mn], bmd[1], c[1], d[1]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m- DRM.QE4_NI(c(a[mn], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*n*
                ((m - DRM.QE4_NI(c(a[1], bmd[1], c[1], d[1]),x,qval))^2)*exp(s[1]))
  }


  return(ll)

}

#' @rdname llfE4_NI_Cov
#' @export
llfP4_NI_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                        nlevels, trt_ind){

  covar <- match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    # c <- pars[grep("c\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }

  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m-DRM.P4_NI(c(a[1], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[1]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m - DRM.P4_NI(c(a[mn], bmd[1], c[1], d[1]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m- DRM.P4_NI(c(a[mn], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*n*
                ((m - DRM.P4_NI(c(a[1], bmd[1], c[1], d[1]),x,qval))^2)*exp(s[1]))
  }

  return(ll)

}

#' @rdname llfE4_NI_Cov
#' @export
llfL4_NI_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                        nlevels, trt_ind){

  covar <- match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    # c <- pars[grep("c\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }

  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m-DRM.L4_NI(c(a[1], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[1]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m - DRM.L4_NI(c(a[mn], bmd[1], c[1], d[1]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m- DRM.L4_NI(c(a[mn], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*n*
                ((m - DRM.L4_NI(c(a[1], bmd[1], c[1], d[1]),x,qval))^2)*exp(s[1]))
  }


  return(ll)

}

#' @rdname llfE4_NI_Cov
#' @export
llfE4_LNI_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                         nlevels, trt_ind){

  covar = match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    # c <- pars[grep("c\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }

  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m+shift - DRM.E4_LNI(c(a[1], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[1]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.E4_LNI(c(a[mn], bmd[1], c[1], d[1]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.E4_LNI(c(a[mn], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*
                n*((m+shift - DRM.E4_LNI(c(a[1], bmd[1], c[1], d[1]),x,qval, shift))^2)*exp(s[1])) - sum((m+shift)*n)

  }


  return(ll)
}


#' @rdname llfE4_NI_Cov
#' @export
llfIE4_LNI_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                          nlevels, trt_ind){

  covar = match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }


  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m+shift - DRM.IE4_LNI(c(a[1], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[1]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.IE4_LNI(c(a[mn], bmd[1], c[1], d[1]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.IE4_LNI(c(a[mn], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*
                n*((m+shift - DRM.IE4_LNI(c(a[1], bmd[1], c[1], d[1]),x,qval, shift))^2)*exp(s[1])) - sum((m+shift)*n)

  }

  return(ll)
}

#' @rdname llfE4_NI_Cov
#' @export
llfH4_LNI_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                         nlevels, trt_ind){

  covar = match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }


  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m+shift - DRM.H4_LNI(c(a[1], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[1]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.H4_LNI(c(a[mn], bmd[1], c[1], d[1]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.H4_LNI(c(a[mn], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*
                n*((m+shift - DRM.H4_LNI(c(a[1], bmd[1], c[1], d[1]),x,qval, shift))^2)*exp(s[1])) - sum((m+shift)*n)

  }

  return(ll)
}


#' @rdname llfE4_NI_Cov
#' @export
llfLN4_LNI_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                          nlevels, trt_ind){

  covar = match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }


  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m+shift - DRM.LN4_LNI(c(a[1], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[1]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.LN4_LNI(c(a[mn], bmd[1], c[1], d[1]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.LN4_LNI(c(a[mn], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*
                n*((m+shift - DRM.LN4_LNI(c(a[1], bmd[1], c[1], d[1]),x,qval,shift))^2)*exp(s[1])) - sum((m+shift)*n)

  }

  return(ll)
}

#' @rdname llfE4_NI_Cov
#' @export
llfG4_LNI_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                         nlevels, trt_ind){

  covar = match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }


  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m+shift - DRM.G4_LNI(c(a[1], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[1]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.G4_LNI(c(a[mn], bmd[1], c[1], d[1]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.G4_LNI(c(a[mn], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*
                n*((m+shift - DRM.G4_LNI(c(a[1], bmd[1], c[1], d[1]),x,qval,shift))^2)*exp(s[1])) - sum((m+shift)*n)

  }

  return(ll)
}

#' @rdname llfE4_NI_Cov
#' @export
llfQE4_LNI_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                          nlevels, trt_ind){

  covar = match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }


  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m+shift - DRM.QE4_LNI(c(a[1], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[1]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.QE4_LNI(c(a[mn], bmd[1], c[1], d[1]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.QE4_LNI(c(a[mn], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*
                n*((m+shift - DRM.QE4_LNI(c(a[1], bmd[1], c[1], d[1]),x,qval,shift))^2)*exp(s[1])) - sum((m+shift)*n)

  }

  return(ll)
}


#' @rdname llfE4_NI_Cov
#' @export
llfP4_LNI_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                         nlevels, trt_ind){

  covar = match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }


  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m+shift - DRM.P4_LNI(c(a[1], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[1]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.P4_LNI(c(a[mn], bmd[1], c[1], d[1]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.P4_LNI(c(a[mn], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*
                n*((m+shift - DRM.P4_LNI(c(a[1], bmd[1], c[1], d[1]),x,qval,shift))^2)*exp(s[1])) - sum((m+shift)*n)

  }

  return(ll)
}


#' @rdname llfE4_NI_Cov
#' @export
llfL4_LNI_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                         nlevels, trt_ind){

  covar = match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }


  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m+shift - DRM.L4_LNI(c(a[1], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[1]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.L4_LNI(c(a[mn], bmd[1], c[1], d[1]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.L4_LNI(c(a[mn], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*
                n*((m+shift - DRM.L4_LNI(c(a[1], bmd[1], c[1], d[1]),x,qval,shift))^2)*exp(s[1])) - sum((m+shift)*n)

  }

  return(ll)
}


#' @rdname llfE4_NI_Cov
#' @export
llfE4_ND_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                        nlevels, trt_ind){

  covar = match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    # c <- pars[grep("c\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }

  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m - DRM.E4_ND(c(a[1], bmd[mn], c[1], d[mn]), x, qval))^2)*
                             exp(s[1]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m - DRM.E4_ND(c(a[mn], bmd[1], c[1], d[1]), x, qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m - DRM.E4_ND(c(a[mn], bmd[mn], c[1], d[mn]), x, qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*
                n*((m - DRM.E4_ND(c(a[1], bmd[1], c[1], d[1]),x,qval))^2)*exp(s[1]))

  }


  return(ll)
}

#' @rdname llfE4_NI_Cov
#' @export
llfIE4_ND_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                         nlevels, trt_ind){

  covar <- match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    # c <- pars[grep("c\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }

  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m - DRM.IE4_ND(c(a[1], bmd[mn], c[1], d[mn]), x, qval))^2)*
                             exp(s[1]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m - DRM.IE4_ND(c(a[mn], bmd[1], c[1], d[1]), x, qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m - DRM.IE4_ND(c(a[mn], bmd[mn], c[1], d[mn]), x, qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*n*
                ((m - DRM.IE4_ND(c(a[1], bmd[1], c[1], d[1]),x,qval))^2)*exp(s[1]))
  }


  return(ll)

}

#' @rdname llfE4_NI_Cov
#' @export
llfH4_ND_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                        nlevels, trt_ind){

  covar <- match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    # c <- pars[grep("c\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }

  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m-DRM.H4_ND(c(a[1], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[1]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m-DRM.H4_ND(c(a[mn], bmd[1], c[1], d[1]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m-DRM.H4_ND(c(a[mn], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*n*
                ((m - DRM.H4_ND(c(a[1], bmd[1], c[1], d[1]),x,qval))^2)*exp(s[1]))
  }

  return(ll)

}


#' @rdname llfE4_NI_Cov
#' @export
llfLN4_ND_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                         nlevels, trt_ind){

  covar <- match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    # c <- pars[grep("c\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }

  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m-DRM.LN4_ND(c(a[1], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[1]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m - DRM.LN4_ND(c(a[mn], bmd[1], c[1], d[1]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m-DRM.LN4_ND(c(a[mn], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*n*
                ((m - DRM.LN4_ND(c(a[1], bmd[1], c[1], d[1]),x,qval))^2)*exp(s[1]))
  }


  return(ll)

}

#' @rdname llfE4_NI_Cov
#' @export
llfG4_ND_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                        nlevels, trt_ind){

  covar <- match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    # c <- pars[grep("c\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }

  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m-DRM.G4_ND(c(a[1], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[1]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m - DRM.G4_ND(c(a[mn], bmd[1], c[1], d[1]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m- DRM.G4_ND(c(a[mn], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*n*
                ((m - DRM.G4_ND(c(a[1], bmd[1], c[1], d[1]),x,qval))^2)*exp(s[1]))
  }


  return(ll)

}

#' @rdname llfE4_NI_Cov
#' @export
llfQE4_ND_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                         nlevels, trt_ind){

  covar <- match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    # c <- pars[grep("c\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }

  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m-DRM.QE4_ND(c(a[1], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[1]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m - DRM.QE4_ND(c(a[mn], bmd[1], c[1], d[1]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m- DRM.QE4_ND(c(a[mn], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*n*
                ((m - DRM.QE4_ND(c(a[1], bmd[1], c[1], d[1]),x,qval))^2)*exp(s[1]))
  }


  return(ll)

}

#' @rdname llfE4_NI_Cov
#' @export
llfP4_ND_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                        nlevels, trt_ind){

  covar <- match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    # c <- pars[grep("c\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }

  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m-DRM.P4_ND(c(a[1], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[1]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m - DRM.P4_ND(c(a[mn], bmd[1], c[1], d[1]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m- DRM.P4_ND(c(a[mn], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*n*
                ((m - DRM.P4_ND(c(a[1], bmd[1], c[1], d[1]),x,qval))^2)*exp(s[1]))
  }

  return(ll)

}

#' @rdname llfE4_NI_Cov
#' @export
llfL4_ND_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                        nlevels, trt_ind){

  covar <- match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    # c <- pars[grep("c\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }

  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m-DRM.L4_ND(c(a[1], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[1]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m - DRM.L4_ND(c(a[mn], bmd[1], c[1], d[1]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m- DRM.L4_ND(c(a[mn], bmd[mn], c[1], d[mn]),x,qval))^2)*
                             exp(s[mn]))*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*n*
                ((m - DRM.L4_ND(c(a[1], bmd[1], c[1], d[1]),x,qval))^2)*exp(s[1]))
  }


  return(ll)

}

#' @rdname llfE4_NI_Cov
#' @export
llfE4_LND_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                         nlevels, trt_ind){

  covar = match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    # c <- pars[grep("c\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }

  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m+shift - DRM.E4_LND(c(a[1], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[1]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.E4_LND(c(a[mn], bmd[1], c[1], d[1]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.E4_LND(c(a[mn], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*
                n*((m+shift - DRM.E4_LND(c(a[1], bmd[1], c[1], d[1]),x,qval, shift))^2)*exp(s[1])) - sum((m+shift)*n)

  }


  return(ll)
}


#' @rdname llfE4_NI_Cov
#' @export
llfIE4_LND_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                          nlevels, trt_ind){

  covar = match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }


  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m+shift - DRM.IE4_LND(c(a[1], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[1]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.IE4_LND(c(a[mn], bmd[1], c[1], d[1]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.IE4_LND(c(a[mn], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*
                n*((m+shift - DRM.IE4_LND(c(a[1], bmd[1], c[1], d[1]),x,qval, shift))^2)*exp(s[1])) - sum((m+shift)*n)

  }

  return(ll)
}

#' @rdname llfE4_NI_Cov
#' @export
llfH4_LND_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                         nlevels, trt_ind){

  covar = match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }


  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m+shift - DRM.H4_LND(c(a[1], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[1]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.H4_LND(c(a[mn], bmd[1], c[1], d[1]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.H4_LND(c(a[mn], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*
                n*((m+shift - DRM.H4_LND(c(a[1], bmd[1], c[1], d[1]),x,qval, shift))^2)*exp(s[1])) - sum((m+shift)*n)

  }

  return(ll)
}


#' @rdname llfE4_NI_Cov
#' @export
llfLN4_LND_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                          nlevels, trt_ind){

  covar = match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }


  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m+shift - DRM.LN4_LND(c(a[1], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[1]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.LN4_LND(c(a[mn], bmd[1], c[1], d[1]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.LN4_LND(c(a[mn], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*
                n*((m+shift - DRM.LN4_LND(c(a[1], bmd[1], c[1], d[1]),x,qval,shift))^2)*exp(s[1])) - sum((m+shift)*n)

  }

  return(ll)
}

#' @rdname llfE4_NI_Cov
#' @export
llfG4_LND_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                         nlevels, trt_ind){

  covar = match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }


  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m+shift - DRM.G4_LND(c(a[1], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[1]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.G4_LND(c(a[mn], bmd[1], c[1], d[1]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.G4_LND(c(a[mn], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*
                n*((m+shift - DRM.G4_LND(c(a[1], bmd[1], c[1], d[1]),x,qval,shift))^2)*exp(s[1])) - sum((m+shift)*n)

  }

  return(ll)
}

#' @rdname llfE4_NI_Cov
#' @export
llfQE4_LND_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                          nlevels, trt_ind){

  covar = match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }


  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m+shift - DRM.QE4_LND(c(a[1], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[1]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.QE4_LND(c(a[mn], bmd[1], c[1], d[1]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.QE4_LND(c(a[mn], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*
                n*((m+shift - DRM.QE4_LND(c(a[1], bmd[1], c[1], d[1]),x,qval,shift))^2)*exp(s[1])) - sum((m+shift)*n)

  }

  return(ll)
}


#' @rdname llfE4_NI_Cov
#' @export
llfP4_LND_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                         nlevels, trt_ind){

  covar = match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }


  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m+shift - DRM.P4_LND(c(a[1], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[1]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.P4_LND(c(a[mn], bmd[1], c[1], d[1]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.P4_LND(c(a[mn], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*
                n*((m+shift - DRM.P4_LND(c(a[1], bmd[1], c[1], d[1]),x,qval,shift))^2)*exp(s[1])) - sum((m+shift)*n)

  }

  return(ll)
}


#' @rdname llfE4_NI_Cov
#' @export
llfL4_LND_Cov = function(pars, x, n, m, s2, qval, shift, covar = c('a_sigma2', 'BMD_d', 'all', 'none'),
                         nlevels, trt_ind){

  covar = match.arg(covar, c('a_sigma2', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    #b <- pars[grep("b\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    c <- pars[nmpar == "par3"]
    d <- pars[grep("par4\\[", nmpar)]
    s <- pars[grep("par5\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    #b <- pars[nmpar == "b"]
    bmd <- pars[nmpar == "par2"]
    c <- pars[nmpar == "par3"]
    d <- pars[nmpar == "par4"]
    s <- pars[nmpar == "par5"]
  }


  ll.level <- numeric(nlevels)
  if(covar == 'BMD_d'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[1]-
                             0.5*(n-1)*s2*exp(s[1])-
                             0.5*n*((m+shift - DRM.L4_LND(c(a[1], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[1]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else if(covar == 'a_sigma2'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.L4_LND(c(a[mn], bmd[1], c[1], d[1]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])

    }
    ll <- sum(ll.level)

  } else if(covar == 'all'){

    for(mn in 1:nlevels){
      ll.level[mn] <- sum((-0.5*n*log(2*pi)+0.5*n*s[mn]-
                             0.5*(n-1)*s2*exp(s[mn])-
                             0.5*n*((m+shift - DRM.L4_LND(c(a[mn], bmd[mn], c[1], d[mn]), x, qval, shift))^2)*
                             exp(s[mn]))*trt_ind[,mn]) - sum((m+shift)*n*trt_ind[,mn])
    }
    ll <- sum(ll.level)

  } else{

    ll <- sum(-0.5*n*log(2*pi)+0.5*n*s[1]-0.5*(n-1)*s2*exp(s[1])-0.5*
                n*((m+shift - DRM.L4_LND(c(a[1], bmd[1], c[1], d[1]),x,qval,shift))^2)*exp(s[1])) - sum((m+shift)*n)

  }

  return(ll)
}


#' Loglikelihood functions for the different dose response models
#'
#' @param pars parameter values
#' @param x vector containing the unique ordered dose levels
#' @param n vector containing the number of observations at each dose level
#' @param y vector containing the number of adverse events at each dose level
#' @param qval BMR
#' @param covar which parameter includes a covariate effect
#' @param nlevels number of covariate levels
#' @param trt_ind matrix indicating which reponse corresponds to which covariate level
#'
#' @examples
#'
#' @return .
#'
#' @export
#'
llfE4_Q_Cov = function(pars, x, n, y, qval, covar = c('background', 'BMD_d', 'all', 'none'),
                       nlevels, trt_ind){

  covar = match.arg(covar, c('background', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    d <- pars[grep("par3\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    bmd <- pars[nmpar == "par2"]
    d <- pars[nmpar == "par3"]
  }

  ll.level <- numeric(nlevels)

  if(covar == 'BMD_d'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.E4_Q(c(a[1], bmd[mn], d[mn]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.E4_Q(c(a[1], bmd[mn], d[mn]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else if(covar == 'background'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.E4_Q(c(a[mn], bmd[1], d[1]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.E4_Q(c(a[mn], bmd[1], d[1]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else if(covar == 'all'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.E4_Q(c(a[mn], bmd[mn], d[mn]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.E4_Q(c(a[mn], bmd[mn], d[mn]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else{
    ll <- sum(lchoose(n, y) + y*log(DRM.E4_Q(c(a[1], bmd[1], d[1]), x, qval) + 1.0E-05) +
                (n - y)*log(1 - DRM.E4_Q(c(a[1], bmd[1], d[1]), x, qval) + 1.0E-05))
  }

  return(ll)
}

#' @rdname llfE4_Q_Cov
#' @export
llfIE4_Q_Cov = function(pars, x, n, y, qval, covar = c('background', 'BMD_d', 'all', 'none'),
                        nlevels, trt_ind){

  covar = match.arg(covar, c('background', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    d <- pars[grep("par3\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    bmd <- pars[nmpar == "par2"]
    d <- pars[nmpar == "par3"]
  }

  ll.level <- numeric(nlevels)

  if(covar == 'BMD_d'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.IE4_Q(c(a[1], bmd[mn], d[mn]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.IE4_Q(c(a[1], bmd[mn], d[mn]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else if(covar == 'background'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.IE4_Q(c(a[mn], bmd[1], d[1]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.IE4_Q(c(a[mn], bmd[1], d[1]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else if(covar == 'all'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.IE4_Q(c(a[mn], bmd[mn], d[mn]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.IE4_Q(c(a[mn], bmd[mn], d[mn]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else{
    ll <- sum(lchoose(n, y) + y*log(DRM.IE4_Q(c(a[1], bmd[1], d[1]), x, qval) + 1.0E-05) +
                (n - y)*log(1 - DRM.IE4_Q(c(a[1], bmd[1], d[1]), x, qval) + 1.0E-05))
  }

  return(ll)
}

#' @rdname llfE4_Q_Cov
#' @export
llfH4_Q_Cov = function(pars, x, n, y, qval, covar = c('background', 'BMD_d', 'all', 'none'),
                       nlevels, trt_ind){

  covar = match.arg(covar, c('background', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    d <- pars[grep("par3\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    bmd <- pars[nmpar == "par2"]
    d <- pars[nmpar == "par3"]
  }

  ll.level <- numeric(nlevels)

  if(covar == 'BMD_d'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.H4_Q(c(a[1], bmd[mn], d[mn]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.H4_Q(c(a[1], bmd[mn], d[mn]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else if(covar == 'background'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.H4_Q(c(a[mn], bmd[1], d[1]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.H4_Q(c(a[mn], bmd[1], d[1]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else if(covar == 'all'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.H4_Q(c(a[mn], bmd[mn], d[mn]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.H4_Q(c(a[mn], bmd[mn], d[mn]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else{
    ll <- sum(lchoose(n, y) + y*log(DRM.H4_Q(c(a[1], bmd[1], d[1]), x, qval) + 1.0E-05) +
                (n - y)*log(1 - DRM.H4_Q(c(a[1], bmd[1], d[1]), x, qval) + 1.0E-05))
  }

  return(ll)
}

#' @rdname llfE4_Q_Cov
#' @export
llfLN4_Q_Cov = function(pars, x, n, y, qval, covar = c('background', 'BMD_d', 'all', 'none'),
                        nlevels, trt_ind){

  covar = match.arg(covar, c('background', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    d <- pars[grep("par3\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    bmd <- pars[nmpar == "par2"]
    d <- pars[nmpar == "par3"]
  }

  ll.level <- numeric(nlevels)

  if(covar == 'BMD_d'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.LN4_Q(c(a[1], bmd[mn], d[mn]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.LN4_Q(c(a[1], bmd[mn], d[mn]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else if(covar == 'background'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.LN4_Q(c(a[mn], bmd[1], d[1]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.LN4_Q(c(a[mn], bmd[1], d[1]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else if(covar == 'all'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.LN4_Q(c(a[mn], bmd[mn], d[mn]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.LN4_Q(c(a[mn], bmd[mn], d[mn]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else{
    ll <- sum(lchoose(n, y) + y*log(DRM.LN4_Q(c(a[1], bmd[1], d[1]), x, qval) + 1.0E-05) +
                (n - y)*log(1 - DRM.LN4_Q(c(a[1], bmd[1], d[1]), x, qval) + 1.0E-05))
  }

  return(ll)
}

#' @rdname llfE4_Q_Cov
#' @export
llfG4_Q_Cov = function(pars, x, n, y, qval, covar = c('background', 'BMD_d', 'all', 'none'),
                       nlevels, trt_ind){

  covar = match.arg(covar, c('background', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    d <- pars[grep("par3\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    bmd <- pars[nmpar == "par2"]
    d <- pars[nmpar == "par3"]
  }

  ll.level <- numeric(nlevels)

  if(covar == 'BMD_d'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.G4_Q(c(a[1], bmd[mn], d[mn]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.G4_Q(c(a[1], bmd[mn], d[mn]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else if(covar == 'background'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.G4_Q(c(a[mn], bmd[1], d[1]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.G4_Q(c(a[mn], bmd[1], d[1]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else if(covar == 'all'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.G4_Q(c(a[mn], bmd[mn], d[mn]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.G4_Q(c(a[mn], bmd[mn], d[mn]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else{
    ll <- sum(lchoose(n, y) + y*log(DRM.G4_Q(c(a[1], bmd[1], d[1]), x, qval) + 1.0E-05) +
                (n - y)*log(1 - DRM.G4_Q(c(a[1], bmd[1], d[1]), x, qval) + 1.0E-05))
  }

  return(ll)
}

#' @rdname llfE4_Q_Cov
#' @export
llfQE4_Q_Cov = function(pars, x, n, y, qval, covar = c('background', 'BMD_d', 'all', 'none'),
                        nlevels, trt_ind){

  covar = match.arg(covar, c('background', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    d <- pars[grep("par3\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    bmd <- pars[nmpar == "par2"]
    d <- pars[nmpar == "par3"]
  }

  ll.level <- numeric(nlevels)

  if(covar == 'BMD_d'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.QE4_Q(c(a[1], bmd[mn], d[mn]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.QE4_Q(c(a[1], bmd[mn], d[mn]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else if(covar == 'background'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.QE4_Q(c(a[mn], bmd[1], d[1]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.QE4_Q(c(a[mn], bmd[1], d[1]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else if(covar == 'all'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.QE4_Q(c(a[mn], bmd[mn], d[mn]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.QE4_Q(c(a[mn], bmd[mn], d[mn]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else{
    ll <- sum(lchoose(n, y) + y*log(DRM.QE4_Q(c(a[1], bmd[1], d[1]), x, qval) + 1.0E-05) +
                (n - y)*log(1 - DRM.QE4_Q(c(a[1], bmd[1], d[1]), x, qval) + 1.0E-05))
  }

  return(ll)
}

#' @rdname llfE4_Q_Cov
#' @export
llfP4_Q_Cov = function(pars, x, n, y, qval, covar = c('background', 'BMD_d', 'all', 'none'),
                       nlevels, trt_ind){

  covar = match.arg(covar, c('background', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    d <- pars[grep("par3\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    bmd <- pars[nmpar == "par2"]
    d <- pars[nmpar == "par3"]
  }

  ll.level <- numeric(nlevels)

  if(covar == 'BMD_d'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.P4_Q(c(a[1], bmd[mn], d[mn]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.P4_Q(c(a[1], bmd[mn], d[mn]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else if(covar == 'background'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.P4_Q(c(a[mn], bmd[1], d[1]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.P4_Q(c(a[mn], bmd[1], d[1]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else if(covar == 'all'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.P4_Q(c(a[mn], bmd[mn], d[mn]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.P4_Q(c(a[mn], bmd[mn], d[mn]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else{
    ll <- sum(lchoose(n, y) + y*log(DRM.P4_Q(c(a[1], bmd[1], d[1]), x, qval) + 1.0E-05) +
                (n - y)*log(1 - DRM.P4_Q(c(a[1], bmd[1], d[1]), x, qval) + 1.0E-05))
  }

  return(ll)
}

#' @rdname llfE4_Q_Cov
#' @export
llfL4_Q_Cov = function(pars, x, n, y, qval, covar = c('background', 'BMD_d', 'all', 'none'),
                       nlevels, trt_ind){

  covar = match.arg(covar, c('background', 'BMD_d', 'all', 'none'))
  nmpar <- names(pars)

  if(covar != 'none'){

    a <- pars[grep("par1\\[", nmpar)]
    bmd <- pars[grep("par2\\[", nmpar)]
    d <- pars[grep("par3\\[", nmpar)]

  } else {
    a <- pars[nmpar == "par1"]
    bmd <- pars[nmpar == "par2"]
    d <- pars[nmpar == "par3"]
  }

  ll.level <- numeric(nlevels)

  if(covar == 'BMD_d'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.L4_Q(c(a[1], bmd[mn], d[mn]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.L4_Q(c(a[1], bmd[mn], d[mn]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else if(covar == 'background'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.L4_Q(c(a[mn], bmd[1], d[1]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.L4_Q(c(a[mn], bmd[1], d[1]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else if(covar == 'all'){
    for(mn in 1:nlevels){
      ll.level[mn] <- sum(
        (lchoose(n, y) + y*log(DRM.L4_Q(c(a[mn], bmd[mn], d[mn]), x, qval) + 1.0E-05) +
           (n - y)*log(1 - DRM.L4_Q(c(a[mn], bmd[mn], d[mn]), x, qval) + 1.0E-05)) * trt_ind[,mn]
      )
    }
    ll <- sum(ll.level)
  }else{
    ll <- sum(lchoose(n, y) + y*log(DRM.L4_Q(c(a[1], bmd[1], d[1]), x, qval) + 1.0E-05) +
                (n - y)*log(1 - DRM.L4_Q(c(a[1], bmd[1], d[1]), x, qval) + 1.0E-05))
  }

  return(ll)
}
