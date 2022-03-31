#' Function to convert arithmetic to geometric mean/standard deviation
#'
#'
#'
#' @param am Arithmetic means per odered dose level, on original scale
#' @param asd Arithmetic standard deviations per odered dose level, on original scale
#'
#' @description Converst arithmetic mean (sd) on original scale to geometric mean (sd) to be used in LNI method
#'
#' @examples
#'  data("immunotoxicityData.rda")
#'  NtoLN(am = immunotoxicityData$Mean, asd = immunotoxicityData$SD)
#'
#' @return Vector containing the geometric means and standard deviations per ordered dose level.
#'
NtoLN=function(am,asd){
  gm=am/sqrt(asd^2/am^2+1)
  gsd=exp(sqrt(log(asd^2/am^2+1)))
  return(c(gm,gsd))
}

#' Function to convert geometric to arithmetic mean/standard deviation
#'
#' Converst geometric mean (sd) to arithmetic mean (sd)
#'
#' @param gm Geometric means per odered dose level, on original scale
#' @param gsd Geometric standard deviations per odered dose level, on original scale
#'
#' @description Converst arithmetic mean (sd) on original scale to geometric mean (sd) to be used in LNI method
#'
#' @examples
#'  data("immunotoxicityData.rda")
#'  LNtoN(gm = immunotoxicityData$Mean, gsd = immunotoxicityData$SD)
#'
#' @return Vector containing the arithmetic means and standard deviations per ordered dose level.

LNtoN=function(gm,gsd){
  am=exp(log(gm)+log(gsd)^2/2)
  asd=sqrt((exp(log(gsd)^2)-1)*exp(2*log(gm)+log(gsd)^2))
  return(c(am,asd))
}


