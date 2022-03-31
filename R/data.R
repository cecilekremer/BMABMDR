

#' immunotoxicity dataset. Effects on Th17 cells
#'
#' Th17 cell frequency in the spleen in offspring mice (%) was analysed. In the paper, the data are
#' represented in a graph and the actual numbers were provided to EFSA by T. Shen. Based on the
#' provided information, EFSA calculated the standard deviation (SD).
#'
#' @format A data frame with 20 rows and 5 variables:
#' \describe{
#'   \item{Dose}{dosage}
#'   \item{Mean}{average number of cells per dose}
#'   \item{SD}{standard deviation per dose}
#'   \item{n}{number of observations in the dose group}
#'   \item{group}{animal group (sex)}
#' }
#' @source Luo SM, Li Y, Li YP, Zhu QX, Jiang JH, Wu CH and Shen T, 2016. Gestational and lactational exposure
#'         to low-dose bisphenol A increases Th17 cells in mice offspring. Environmental Toxicology and
#'         Pharmacology, 47, 149–158. doi: 10.1016/j.etap.2016.09.017. RefID: 4679.
"immunotoxicityData"



#' Learning and memory during the growth phase/young age exposure period
#'
#'
#'
#' @format A data frame with 4 rows and 4 variables:
#' \describe{
#'   \item{Dose}{dosage}
#'   \item{Mean}{average per dose}
#'   \item{SEM}{standard error of the mean per dose}
#'   \item{n}{number of observations in the dose group}
#' }
#' @source Chen Z, Li T, Zhang L, Wang H and Hu F, 2018. Bisphenol A exposure remodels cognition of male rats
#'         attributable to excitatory alterations in the hippocampus and visual cortex. Toxicology, 410, 132–
#'         doi: 10.1016/j.tox.2018.10.002. RefID: 11734.
"LearningMemory"
