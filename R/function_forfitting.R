




#' The likelihood function
#'
#' @description
#' This function (i.e., \code{overalllikelihood()}) calculates the likelihood value with a list of pre-defined epidemiological parameters and a given structured contact tracing data.
#'
#' @param epi.para A list (\code{list}) of pre-defined epidemiological parameters for offspring distribution, in the format of \code{list(mean = ?, disp = ?, shift = ?)},
#' where the three parameters accept non-negative values.
#' Each parameter must be a scalar.
#' For Delaporte distribution, the value of \code{mean} should be larger than the value of \code{shift}.
#' @param offspring.type
#' A character label (\code{character}) indicating the type of distribution used to describe the offspring distribution.
#' It only accepts one of the following values:
#' \itemize{
#'   \item{\code{"D"}}{ indicates the *Delaporte* distribution, }
#'   \item{\code{"NB"}}{ indicates the *negative binomial* distribution, }
#'   \item{\code{"G"}}{ indicates the *geometric* distribution, or }
#'   \item{\code{"P"}}{ indicates the *Poisson* distribution. }
#' }
#' By default, \code{offspring.type = 'D'}.
#' @param is.log A logical variable, under which the likelihood would be taken natural logarithm, i.e., log-likelihood, if \code{is.log = TRUE}.
#' By default, \code{is.log = TRUE}.
#' @param data A data frame (\code{data.frame}), or a vector (only when \code{obs.type.lab = "offspring"}) that contains the structured contact tracing data.
#' @param var.name A list (\code{list}), or a character of variable name for the column names of dataset given in \code{data}.
#' For a list of variable names, it should be in the format of \code{list(obssize = ?, seedsize = ?, typelab = ?)}.
#' Please see the details section for more information.
#' By default, \code{var.name = list(obssize = NULL, seedsize = NULL, typelab = NULL)}.
#' @param obs.type.lab A list (\code{list}), or a character of labels (i.e., "offspring", "nextgen", or "outbreak") for the type of observations.
#' For a list of labels, it should be in the format of \code{list(offspring = ?, nextgen = ?, outbreak = ?)}.
#' Please see the details section for more information.
#' By default, \code{obs.type.lab = list(offspring = NULL, nextgen = NULL, outbreak = NULL)}.
#'
#'
#' @return
#' The log-likelihood (by default), or likelihood value from contact tracing data, with pre-defined epidemiological parameters.
#'
#' @details
#' When \code{obs.type.lab} is a character, it should be either \code{"offspring"}, \code{"nextgen"}, or \code{"outbreak"} for type of observations.
#'
#' When \code{obs.type.lab} is a list, this occurs when the contact tracing data has more than one types of observations.
#' See example 4 in the Examples section.
#'
#' When the contact tracing dataset is offspring case observations, the function arguments \code{data} could be either a vector, or a data frame.
#' If \code{data} is a vector, it is not necessary to assign any value to \code{var.name}.
#' If \code{data} is a data frame, it is necessary to identify the variable name of offspring observations in \code{var.name}.
#' See example 1 in the Examples section.
#'
#' When the contact tracing dataset is next-generation cluster size, or final outbreak size observations, the variable names of both observations and seed case size should be identified in \code{var.name} with the format of \code{list(obssize = ?, seedsize = ?)}.
#' See example 2 and example 3 in the Examples section.
#'
#' When the contact tracing dataset has more than one types of observations, the variable names of observations, seed case size, and observation type should be identified in \code{var.name} with the format of \code{list(obssize = ?, seedsize = ?, typelab = ?)}.
#' See example 4 in the Examples section.
#'
#'
#' @import Delaporte
#'
#' @export
#' overalllikelihood
#'
#' @note
#' Each parameter in \code{epi.para = list(mean = ?, disp = ?, shift = ?)} should be a scalar, which means vector is not allowed here.
#'
#' For the contact tracing data in \code{data}, unknown observations (i.e., \code{NA}) is not allowed.
#'
#'
#' @references
#' Lloyd-Smith JO, Schreiber SJ, Kopp PE, Getz WM. Superspreading and the effect of individual variation on disease emergence. *Nature*. 2005;438(7066):355-359.
#' \doi{10.1038/nature04153}
#'
#' Nishiura H, Yan P, Sleeman CK, Mode CJ. Estimating the transmission potential of supercritical processes based on the final size distribution of minor outbreaks. *Journal of Theoretical Biology*. 2012;294:48-55.
#' \doi{10.1016/j.jtbi.2011.10.039}
#'
#' Blumberg S, Funk S, Pulliam JR. Detecting differential transmissibilities that affect the size of self-limited outbreaks. *PLoS Pathogens*. 2014;10(10):e1004452.
#' \doi{10.1371/journal.ppat.1004452}
#'
#' Kucharski AJ, Althaus CL. The role of superspreading in Middle East respiratory syndrome coronavirus (MERS-CoV) transmission. *Eurosurveillance*. 2015;20(25):21167.
#' \doi{10.2807/1560-7917.ES2015.20.25.21167}
#'
#' Endo A, Abbott S, Kucharski AJ, Funk S. Estimating the overdispersion in COVID-19 transmission using outbreak sizes outside China. *Wellcome Open Research*. 2020;5:67.
#' \doi{10.12688/wellcomeopenres.15842.3}
#'
#' Adam DC, Wu P, Wong JY, Lau EH, Tsang TK, Cauchemez S, Leung GM, Cowling BJ. Clustering and superspreading potential of SARS-CoV-2 infections in Hong Kong. *Nature Medicine*. 2020;26(11):1714-1719.
#' \doi{10.1038/s41591-020-1092-0}
#'
#' Zhao S, Chong MK, Ryu S, Guo Z, He M, Chen B, Musa SS, Wang J, Wu Y, He D, Wang MH. Characterizing superspreading potential of infectious disease: Decomposition of individual transmissibility. *PLoS Computational Biology*. 2022;18(6):e1010281.
#' \doi{10.1371/journal.pcbi.1010281}
#'
#'
#' @examples
#'
#' # example 1 #
#' ## likelihood for the offspring observations
#' data(COVID19_JanApr2020_HongKong)
#' overalllikelihood(
#'   epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
#'   offspring.type = "D",
#'   data = COVID19_JanApr2020_HongKong,
#'   var.name = list(obssize = 'obs'),
#'   obs.type.lab = 'offspring'
#' )
#' overalllikelihood(
#'   epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
#'   offspring.type = "D",
#'   data = COVID19_JanApr2020_HongKong$obs,
#'   obs.type.lab = 'offspring'
#' )
#'
#'
#' # example 2 #
#' ## likelihood for the next-generation cluster size observations
#' data(smallpox_19581973_Europe)
#' overalllikelihood(
#'   epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
#'   offspring.type = 'D',
#'   data = smallpox_19581973_Europe,
#'   var.name = list(obssize = 'obs.clustersize', seedsize = 'obs.seed'),
#'   obs.type.lab = 'nextgen'
#' )
#'
#'
#' # example 3 #
#' ## likelihood for the final outbreak size observations
#' data(MERS_2013_MEregion)
#' overalllikelihood(
#'   epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
#'   offspring.type = 'D',
#'   data = MERS_2013_MEregion,
#'   var.name = list(obssize = 'obs.finalsize', seedsize = 'obs.seed'),
#'   obs.type.lab = 'outbreak'
#' )
#'
#'
#' # example 4 #
#' ## likelihood for more than one types of observations
#' data(mpox_19801984_DRC)
#' overalllikelihood(
#'   epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
#'   offspring.type = 'D',
#'   data = mpox_19801984_DRC,
#'   var.name = list(obssize = 'obs.size', seedsize = 'obs.seed', typelab = 'type'),
#'   obs.type.lab = list(offspring = 'offspring', nextgen = 'nextgen', outbreak = 'outbreak')
#' )
#'
#'
#' # example 5 #
#' ## reproducing the AIC results in Adam, et al. (2020) https://doi.org/10.1038/s41591-020-1092-0,
#' ## (see Supplementary Table 4),
#' ## where the AIC scores were calculated for NB, Geometric, and Poisson models from top to bottom.
#' ## Here, the AIC is defined as: AIC = -2 * log-likelihood + 2 * number of unknown model parameters.
#' data(COVID19_JanApr2020_HongKong)
#' overalllikelihood(
#'   epi.para = list(mean = 0.58, disp = 0.43, shift = 0.2),
#'   offspring.type = "NB",
#'   data = COVID19_JanApr2020_HongKong$obs,
#'   obs.type.lab = 'offspring'
#' ) * (-2) + 2*2
#' overalllikelihood(
#'   epi.para = list(mean = 0.63, disp = 0.43, shift = 0.2),
#'   offspring.type = "G",
#'   data = COVID19_JanApr2020_HongKong$obs,
#'   obs.type.lab = 'offspring'
#' ) * (-2) + 1*2
#' overalllikelihood(
#'   epi.para = list(mean = 0.58, disp = 0.43, shift = 0.2),
#'   offspring.type = "P",
#'   data = COVID19_JanApr2020_HongKong$obs,
#'   obs.type.lab = 'offspring'
#' ) * (-2) + 1*2
#'
overalllikelihood = function(
  epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
  offspring.type = 'D', is.log = TRUE,
  data = NULL,
  var.name = list(obssize = NULL, seedsize = NULL, typelab = NULL),
  obs.type.lab = list(offspring = NULL, nextgen = NULL, outbreak = NULL)
){
  raw.data = data
  offspring.likelihood.array = 0
  nextgen.likelihood.array = 0
  outbreak.likelihood.array = 0

  if(is.null(dim(raw.data))){
    if(identical(obs.type.lab, 'offspring')){
      offspring.likelihood.array = d_offspringdistn(
        x = c(raw.data), epi.para = epi.para,
        offspring.type = offspring.type, is.log = TRUE
      )
    } else {
      stop('for observation data in a vector form, the function argument  \"obs.type.lab\"  must be \"offspring\" for offspring type of observations.')
    }
  } else {
    colname.index.array = match(c(unlist(var.name)), colnames(raw.data))
    #colnames(raw.data)[colname.index.array] = c('obssize', 'seedsize', 'typelab')
    colnames(raw.data)[colname.index.array] = names(c(unlist(var.name)))


    if(identical(obs.type.lab, 'offspring')){
      offspring.likelihood.array = d_offspringdistn(
        x = raw.data$obssize, epi.para = epi.para,
        offspring.type = offspring.type, is.log = TRUE
      )
    } else if(identical(obs.type.lab, 'nextgen')){
      nextgen.likelihood.array = d_nextgenclusterdistn(
        x = raw.data$obssize, seed.size = raw.data$seedsize,
        epi.para = epi.para,
        offspring.type = offspring.type, is.log = TRUE
      )
    } else if(identical(obs.type.lab, 'outbreak')){
      outbreak.likelihood.array = d_outbreakdistn(
        x = raw.data$obssize, seed.size = raw.data$seedsize,
        epi.para = epi.para,
        offspring.type = offspring.type, is.log = TRUE
      )
    } else if(length(obs.type.lab) > 1){
      offspring.data = subset(
        x = raw.data, subset = (raw.data$typelab %in% c(unlist(obs.type.lab$offspring)))
      )
      nextgen.data = subset(
        x = raw.data, subset = (raw.data$typelab %in% c(unlist(obs.type.lab$nextgen)))
      )
      outbreak.data = subset(
        x = raw.data, subset = (raw.data$typelab %in% c(unlist(obs.type.lab$outbreak)))
      )

      offspring.likelihood.array = d_offspringdistn(
        x = offspring.data$obssize,
        epi.para = epi.para,
        offspring.type = offspring.type, is.log = TRUE
      )
      nextgen.likelihood.array = d_nextgenclusterdistn(
        x = nextgen.data$obssize, seed.size = nextgen.data$seedsize,
        epi.para = epi.para,
        offspring.type = offspring.type, is.log = TRUE
      )
      outbreak.likelihood.array = d_outbreakdistn(
        x = outbreak.data$obssize, seed.size = outbreak.data$seedsize,
        epi.para = epi.para,
        offspring.type = offspring.type, is.log = TRUE
      )
    } else {
      stop(
        'the value of the function argument  \"obs.type.lab\"  must be:
          -- a list as in the pre-defined format,
          -- \"offspring\" for offspring type of observations,
          -- \"nextgen\" for next-generation cluster size type of observations, or
          -- \"outbreak\" for final outbreak size type of observations.'
      )
    }
  }

  overall.loglikelihood = sum(offspring.likelihood.array, nextgen.likelihood.array, outbreak.likelihood.array)
  overall.likelihood = exp(overall.loglikelihood)
  if(is.log){
    return(overall.loglikelihood)
  } else {
    return(overall.likelihood)
  }
}
















#' To estimate model parameters using maximum likelihood approach
#'
#' @description
#' This function (i.e., \code{paraest.ML()}) performs model parameter estimation using **maximum likelihood** (ML) approach with given structured contact tracing data.
#'
#' @param can.epi.para.range A list (\code{list}) of ranges, or fixed values for unknown epidemiological parameters for offspring distribution.
#' For the ranges of unknown epidemiological parameters, the list should be in the format of \code{list(mean = c(?, ?), disp = c(?, ?), shift = c(?, ?))}.
#' For the fixed values of unknown epidemiological parameters, the list should be in the format of \code{list(mean = ?, disp = ?, shift = ?)}.
#' Each parameter must be a scalar, and only accept non-negative values.
#' The default setting is given in the code Usage section.
#' For Delaporte distribution, the value of \code{mean} should be larger than the value of \code{shift}.
#' @param offspring.type
#' A character label (\code{character}) indicating the type of distribution used to describe the offspring distribution.
#' It only accepts one of the following values:
#' \itemize{
#'   \item{\code{"D"}}{ indicates the *Delaporte* distribution, }
#'   \item{\code{"NB"}}{ indicates the *negative binomial* distribution, }
#'   \item{\code{"G"}}{ indicates the *geometric* distribution, or }
#'   \item{\code{"P"}}{ indicates the *Poisson* distribution. }
#' }
#' By default, \code{offspring.type = 'D'}.
#' @param para.comb.num A positive integer for the number of parameter combinations used to construct log-likelihood profile.
#' By default, \code{para.comb.num = 1000}, and no need to change the default setting here unless for special reasons.
#' @param can.epi.para.set A data frame (\code{data.frame}) of different parameter combinations.
#' The data frame must have three variables with names \code{"epi.para.mean"}, \code{"epi.para.disp"}, and \code{"epi.para.shift"} for the three parameters.
#' By default, \code{can.epi.para.set = NULL}.
#' Note that the function argument \code{can.epi.para.set} is usually used internally, and thus no need to change the default setting here unless for special reasons
#' @param data A data frame (\code{data.frame}), or a vector (only when \code{obs.type.lab = "offspring"}) that contains the structured contact tracing data.
#' @param var.name A list (\code{list}), or a character of variable name for the column names of dataset given in \code{data}.
#' For a list of variable names, it should be in the format of \code{list(obssize = ?, seedsize = ?, typelab = ?)}.
#' Please see the details section for more information.
#' By default, \code{var.name = list(obssize = NULL, seedsize = NULL, typelab = NULL)}.
#' @param obs.type.lab A list (\code{list}), or a character of labels (i.e., "offspring", "nextgen", or "outbreak") for the type of observations.
#' For a list of labels, it should be in the format of \code{list(offspring = ?, nextgen = ?, outbreak = ?)}.
#' Please see the details section for more information.
#' By default, \code{obs.type.lab = list(offspring = NULL, nextgen = NULL, outbreak = NULL)}.
#'
#'
#' @return
#' A list (i.e., \code{list}) contains the following three items:
#' \itemize{
#'   \item{}{a data frame (\code{data.frame}) of the maximum likelihood estimate and 95% confidence interval (CI) of each unknown parameters,}
#'   \item{}{the maximum log-likelihood value, and}
#'   \item{}{a data frame (\code{data.frame}) of different parameter combinations and their corresponding log-likelihood values.}
#' }
#'
#'
#' @details
#' For the ranges of parameters given in \code{can.epi.para.range},
#' they are some rough ranges, which are not necessarily to be precise (but have to be within a reasonable range), and
#' they will used as a "start" status to find the maximum likelihood estimate.
#'
#' When \code{obs.type.lab} is a character, it should be either \code{"offspring"}, \code{"nextgen"}, or \code{"outbreak"} for type of observations.
#' When \code{obs.type.lab} is a list, this occurs when the contact tracing data has more than one types of observations.
#'
#' When the contact tracing dataset is offspring case observations, the function arguments \code{data} could be either a vector, or a data frame.
#' If \code{data} is a vector, it is not necessary to assign any value to \code{var.name}.
#' If \code{data} is a data frame, it is necessary to identify the variable name of offspring observations in \code{var.name}.
#'
#' When the contact tracing dataset is next-generation cluster size, or final outbreak size observations, the variable names of both observations and seed case size should be identified in \code{var.name} with the format of \code{list(obssize = ?, seedsize = ?)}.
#'
#' When the contact tracing dataset has more than one types of observations, the variable names of observations, seed case size, and observation type should be identified in \code{var.name} with the format of \code{list(obssize = ?, seedsize = ?, typelab = ?)}.
#'
#'
#'
#' @import Delaporte
#' @importFrom
#' stats qchisq runif
#'
#' @export
#' paraest.ML
#'
#' @note
#' For the contact tracing data in \code{data}, unknown observations (i.e., \code{NA}) is not allowed.
#'
#' When \code{para.comb.num} is large, e.g., \code{para.comb.num} > 10000, the function \code{paraest.ML()} could take few seconds, or even minutes to complete, depending on the sample size, and model settings, etc.
#' Thus, we do not recommend the users to change the default setting of \code{para.comb.num} unless for special reasons.
#'
#'
#' @seealso
#' \code{\link[modelSSE:overalllikelihood]{modelSSE::overalllikelihood()}}
#'
#'
#' @references
#' Blumberg S, Funk S, Pulliam JR. Detecting differential transmissibilities that affect the size of self-limited outbreaks. *PLoS Pathogens*. 2014;10(10):e1004452.
#' \doi{10.1371/journal.ppat.1004452}
#'
#' Kucharski AJ, Althaus CL. The role of superspreading in Middle East respiratory syndrome coronavirus (MERS-CoV) transmission. *Eurosurveillance*. 2015;20(25):21167.
#' \doi{10.2807/1560-7917.ES2015.20.25.21167}
#'
#' Adam DC, Wu P, Wong JY, Lau EH, Tsang TK, Cauchemez S, Leung GM, Cowling BJ. Clustering and superspreading potential of SARS-CoV-2 infections in Hong Kong. *Nature Medicine*. 2020;26(11):1714-1719.
#' \doi{10.1038/s41591-020-1092-0}
#'
#' Zhao S, Chong MK, Ryu S, Guo Z, He M, Chen B, Musa SS, Wang J, Wu Y, He D, Wang MH. Characterizing superspreading potential of infectious disease: Decomposition of individual transmissibility. *PLoS Computational Biology*. 2022;18(6):e1010281.
#' \doi{10.1371/journal.pcbi.1010281}
#'
#'
#' @examples
#'
#' ## try to estimate the parameter (which is already known),
#' ## using random samples generated from a geometric distribution with mean of 1.
#' set.seed(2020)
#' paraest.ML(
#'   can.epi.para.range = list(mean = c(0.1, 2.0), disp = c(0.01, 2.5), shift = c(0.01,0.5)),
#'   offspring.type = "NB", para.comb.num = 100,
#'   data = r_offspringdistn(
#'     n = 99, epi.para = list(mean = 1, disp = 0.5, shift = 0.2), offspring.type = "G"
#'   ),
#'   obs.type.lab = 'offspring'
#' )$epi.para.est.output
#'
#'
#' \donttest{
#'
#' # example 1: for offspring observations #
#' ## reproducing the parameter estimation results in Adam, et al. (2020)
#' ## paper doi link: https://doi.org/10.1038/s41591-020-1092-0,
#' ## (see the first row in Supplementary Table 4),
#' ## where R of 0.58 (95% CI: 0.45, 0.72), and k of 0.43 (95% CI: 0.29, 0.67).
#' data(COVID19_JanApr2020_HongKong)
#' set.seed(2020)
#' paraest.ML(
#'   can.epi.para.range = list(mean = c(0.1, 2.0), disp = c(0.01, 2.5), shift = c(0.01,0.5)),
#'   offspring.type = "NB",
#'   data = COVID19_JanApr2020_HongKong$obs,
#'   obs.type.lab = 'offspring'
#' )$epi.para.est.output
#'
#'
#' # example 2: for offspring observations #
#' ## reproducing the parameter estimation results in Zhao, et al. (2020)
#' ## paper doi link: https://doi.org/10.1371/journal.pcbi.1010281,
#' ## (see the results of dataset #3 using Delaporte distribution in Table 1), where
#' ## R of 0.59 (95% CI: 0.46, 0.78),
#' ## k of 0.16 (95% CI: 0.06, 0.40), and
#' ## shift of 0.17 (95% CI: 0.04, 0.30).
#' data(COVID19_JanApr2020_HongKong)
#' set.seed(2020)
#' paraest.ML(
#'   can.epi.para.range = list(mean = c(0.1, 2.0), disp = c(0.01, 2.5), shift = c(0.01,0.5)),
#'   offspring.type = "D",
#'   data = COVID19_JanApr2020_HongKong$obs,
#'   obs.type.lab = 'offspring'
#' )$epi.para.est.output
#'
#'
#' # example 3: for next-generation cluster size observations #
#' ## reproducing the parameter estimation results in Blumberg, et al, (2014)
#' ## paper doi link: https://doi.org/10.1371/journal.ppat.1004452,
#' ## (see the last row in Table 3, and Fig 4A),
#' ## where R of 3.14 (95% CI: 2, >6), and k of 0.37 (95% CI: not reported).
#' data(smallpox_19581973_Europe)
#' set.seed(2020)
#' paraest.ML(
#'   can.epi.para.range = list(mean = c(0.1, 10.0), disp = c(0.01, 2.5), shift = c(0.01,0.5)),
#'   offspring.type = "NB",
#'   data = smallpox_19581973_Europe,
#'   var.name = list(obssize = 'obs.clustersize', seedsize = 'obs.seed'),
#'   obs.type.lab = 'nextgen'
#' )$epi.para.est.output
#'
#'
#' # example 4: final outbreak size observations #
#' ## reproducing the parameter estimation results in Kucharski, Althaus. (2015)
#' ## paper doi link: https://doi.org/10.2807/1560-7917.ES2015.20.25.21167,
#' ## (see Fig 1, and Finding section),
#' ## where R of 0.47 (95% CI: 0.29, 0.80), and k of 0.26 (95% CI: 0.09, 1.24).
#' data(MERS_2013_MEregion)
#' set.seed(2020)
#' paraest.ML(
#'   can.epi.para.range = list(mean = c(0.1, 2.0), disp = c(0.01, 2.5), shift = c(0.01,0.5)),
#'   offspring.type = "NB",
#'   data = MERS_2013_MEregion,
#'   var.name = list(obssize = 'obs.finalsize', seedsize = 'obs.seed'),
#'   obs.type.lab = 'outbreak'
#' )$epi.para.est.output
#'
#'
#' # example 5: for more than one types of observations #
#' ## reproducing the parameter estimation results in Blumberg, et al, (2014)
#' ## paper doi link: https://doi.org/10.1371/journal.ppat.1004452,
#' ## (see the last row in Table 5, and Fig 6A),
#' ## where R of 0.3 (95% CI: 0.2, 0.5), and k of 0.4 (95% CI: not reported).
#' data(mpox_19801984_DRC)
#' set.seed(2020)
#' paraest.ML(
#'   can.epi.para.range = list(mean = c(0.1, 2.0), disp = c(0.01, 2.5), shift = c(0.01,0.5)),
#'   offspring.type = "NB",
#'   data = mpox_19801984_DRC,
#'   var.name = list(obssize = 'obs.size', seedsize = 'obs.seed', typelab = 'type'),
#'   obs.type.lab = list(offspring = 'offspring', nextgen = 'nextgen', outbreak = 'outbreak')
#' )$epi.para.est.output
#'
#' }
#'
paraest.ML = function(
  can.epi.para.range = list(mean = c(0.1, 2.0), disp = c(0.01, 2.5), shift = c(0.01,0.5)),
  offspring.type = "D",
  para.comb.num = 1000, can.epi.para.set = NULL,
  data = NULL,
  var.name = list(obssize = NULL, seedsize = NULL, typelab = NULL),
  obs.type.lab = list(offspring = NULL, nextgen = NULL, outbreak = NULL)
){

  if(is.null(can.epi.para.set)){
    ## if no pre-defined data frame for parameter was given in "can.epi.para.set",
    ## try to run the code with ranges of parameter given in "can.epi.para.range", and
    ## try to narrow down the ranges of parameter for the next step.

    #para.seed.num = ceiling(para.comb.num / 44)
    para.seed.num = ceiling(para.comb.num ^(1/3) *2.5)
    #para.seed.num = 22

    if(length(can.epi.para.range$mean) == 2){
      epi.para.mean.array = exp(seq(log(can.epi.para.range$mean[1]), log(can.epi.para.range$mean[2]), length.out = para.seed.num))
      isfixed.epi.para.mean = FALSE
    } else if(length(can.epi.para.range$mean) == 1){
      epi.para.mean.array = unlist(can.epi.para.range$mean)
      isfixed.epi.para.mean = TRUE
    } else if(length(can.epi.para.range$mean) >= 3 | is.null(can.epi.para.range$mean)){
      stop("for the function argument  \"can.epi.para.range\" , the variable  \"mean\"  must be of length 1 or 2, but not larger than 2.")
    }

    if(offspring.type %in% c('G', 'P')){
      epi.para.disp.array = 1
      isfixed.epi.para.disp = TRUE
    } else {
      if(length(can.epi.para.range$disp) == 2){
        epi.para.disp.array = exp(seq(log(can.epi.para.range$disp[1]), log(can.epi.para.range$disp[2]), length.out = para.seed.num))
        isfixed.epi.para.disp = FALSE
      } else if(length(can.epi.para.range$disp) == 1){
        epi.para.disp.array = unlist(can.epi.para.range$disp)
        isfixed.epi.para.disp = TRUE
      } else if(length(can.epi.para.range$disp) >= 3 | is.null(can.epi.para.range$disp)){
        stop("for the function argument  \"can.epi.para.range\" , the variable  \"disp\"  must be of length 1 or 2, but not larger than 2.")
      }
    }

    if(offspring.type %in% c("NB", 'G', 'P')){
      epi.para.shift.array = 1e-09
      isfixed.epi.para.shift = TRUE
    } else {
      if(length(can.epi.para.range$shift) == 2){
        epi.para.shift.array = exp(seq(log(can.epi.para.range$shift[1]), log(can.epi.para.range$shift[2]), length.out = para.seed.num))
        isfixed.epi.para.shift = FALSE
      } else if(length(can.epi.para.range$shift) == 1){
        epi.para.shift.array = unlist(can.epi.para.range$shift)
        isfixed.epi.para.shift = TRUE
      } else if(length(can.epi.para.range$shift) >= 3 | is.null(can.epi.para.range$shift)){
        stop("for the function argument  \"can.epi.para.range\" , the variable  \"shift\"  must be of length 1 or 2, but not larger than 2.")
      }
    }


    record.mat = NULL
    for (epi.para.mean.jj in 1:length(epi.para.mean.array)) {#         epi.para.mean.jj = 1
      here.epi.para.mean = epi.para.mean.array[epi.para.mean.jj]

      for (epi.para.disp.jj in 1:length(epi.para.disp.array)) {#         epi.para.disp.jj = 1
        here.epi.para.disp = epi.para.disp.array[epi.para.disp.jj]

        for (epi.para.shift.jj in 1:length(epi.para.shift.array)) {#         epi.para.shift.jj = 1
          here.epi.para.shift = epi.para.shift.array[epi.para.shift.jj]
          if(here.epi.para.mean <= here.epi.para.shift){
            next
          }

          here.ll = overalllikelihood(
            epi.para = list(mean = here.epi.para.mean, disp = here.epi.para.disp, shift = here.epi.para.shift),
            offspring.type = offspring.type,
            data = data,
            var.name = var.name,
            obs.type.lab = obs.type.lab
          )
          here.record.array = c(here.epi.para.mean, here.epi.para.disp, here.epi.para.shift, here.ll)

          record.mat = rbind(record.mat, c(here.record.array))
        }
      }
    }
    record.mat = as.data.frame(record.mat)
    colnames(record.mat) = c('epi.para.mean', 'epi.para.disp', 'epi.para.shift', 'll')

    max.ll = max(record.mat$ll)
    best.record.array = record.mat[which.max(record.mat$ll),]
    best.epi.para.mean = best.record.array$epi.para.mean[1]
    best.epi.para.disp = best.record.array$epi.para.disp[1]
    best.epi.para.shift = best.record.array$epi.para.shift[1]

    sel.record.mat = subset(record.mat, subset = record.mat$ll >= (max.ll - qchisq(p = 0.95, df = 1) / 1))
    # if(!is.null(can.epi.para.set)){
    #   sel.record.mat = sel.record.mat[, c(!c(isfixed.epi.para.mean, isfixed.epi.para.disp, isfixed.epi.para.shift), T)]
    #   return(sel.record.mat)
    # }


    if(best.epi.para.mean > min(sel.record.mat$epi.para.mean) & best.epi.para.mean < max(sel.record.mat$epi.para.mean)){
      next.range.value = diff(range(sel.record.mat$epi.para.mean))
    } else {
      next.range.value = max(c(0.22,diff(range(sel.record.mat$epi.para.mean))))
    }
    next.epi.para.mean.lwr = best.epi.para.mean - next.range.value *0.66; next.epi.para.mean.lwr = ifelse(next.epi.para.mean.lwr < 0, 1e-02, next.epi.para.mean.lwr)
    next.epi.para.mean.upr = best.epi.para.mean + next.range.value *0.66

    if(best.epi.para.disp > min(sel.record.mat$epi.para.disp) & best.epi.para.disp < max(sel.record.mat$epi.para.disp)){
      next.range.value = diff(range(sel.record.mat$epi.para.disp))
    } else {
      next.range.value = max(c(0.22,diff(range(sel.record.mat$epi.para.disp))))
    }
    next.epi.para.disp.lwr = best.epi.para.disp - next.range.value *0.66; next.epi.para.disp.lwr = ifelse(next.epi.para.disp.lwr < 0, 1e-02, next.epi.para.disp.lwr)
    next.epi.para.disp.upr = best.epi.para.disp + next.range.value *0.66

    if(best.epi.para.shift > min(sel.record.mat$epi.para.shift) & best.epi.para.shift < max(sel.record.mat$epi.para.shift)){
      next.range.value = diff(range(sel.record.mat$epi.para.shift))
    } else {
      next.range.value = max(c(0.22,diff(range(sel.record.mat$epi.para.shift))))
    }
    next.epi.para.shift.lwr = best.epi.para.shift - next.range.value *0.66; next.epi.para.shift.lwr = ifelse(next.epi.para.shift.lwr < 0, 1e-04, next.epi.para.shift.lwr)
    next.epi.para.shift.upr = best.epi.para.shift + next.range.value *0.66

    next.epi.para.set = data.frame(
      epi.para.mean = ifelse(rep(isfixed.epi.para.mean, para.comb.num +1), best.epi.para.mean, c(best.epi.para.mean, exp(runif(n = para.comb.num, min = log(next.epi.para.mean.lwr), max = log(next.epi.para.mean.upr))))),
      epi.para.disp = ifelse(rep(isfixed.epi.para.disp, para.comb.num +1), best.epi.para.disp, c(best.epi.para.disp, exp(runif(n = para.comb.num, min = log(next.epi.para.disp.lwr), max = log(next.epi.para.disp.upr))))),
      epi.para.shift = ifelse(rep(isfixed.epi.para.shift, para.comb.num +1), best.epi.para.shift, c(best.epi.para.shift, exp(runif(n = para.comb.num, min = log(next.epi.para.shift.lwr), max = log(next.epi.para.shift.upr)))))
    )

    paraest.ML(
      can.epi.para.range = list(mean = NULL, disp = NULL, shift = NULL),
      offspring.type = offspring.type, #is.logtranspara = TRUE,
      can.epi.para.set = next.epi.para.set,
      data = data, var.name = var.name, obs.type.lab = obs.type.lab
    )

  } else {
    ## if a pre-defined data frame for parameter was given in "can.epi.para.set",
    ## try to run the code with ranges of parameter given, and return the likelihood profile

    epi.para.mean.array = c(can.epi.para.set$epi.para.mean)
    isfixed.epi.para.mean = FALSE

    if(offspring.type %in% c('G', 'P')){
      epi.para.disp.array = 1
      isfixed.epi.para.disp = TRUE
    } else {
      epi.para.disp.array = c(can.epi.para.set$epi.para.disp)
      isfixed.epi.para.disp = FALSE
    }

    if(offspring.type %in% c("NB", 'G', 'P')){
      epi.para.shift.array = 1e-09
      isfixed.epi.para.shift = TRUE
    } else {
      epi.para.shift.array = c(can.epi.para.set$epi.para.shift)
      isfixed.epi.para.shift = FALSE
    }


    record.mat = NULL
    for(epi.para.jj in 1:nrow(can.epi.para.set)) {#         epi.para.jj = 1
      here.epi.para.array = unlist(can.epi.para.set[epi.para.jj,])
      here.epi.para.mean = here.epi.para.array[1]
      here.epi.para.disp = here.epi.para.array[2]
      here.epi.para.shift = here.epi.para.array[3]

      if(here.epi.para.mean <= here.epi.para.shift){
        next
      }

      here.ll = overalllikelihood(
        epi.para = list(mean = here.epi.para.mean, disp = here.epi.para.disp, shift = here.epi.para.shift),
        offspring.type = offspring.type,
        data = data,
        var.name = var.name,
        obs.type.lab = obs.type.lab
      )
      here.record.array = c(here.epi.para.array, here.ll)

      record.mat = rbind(record.mat, c(here.record.array))
    }
    record.mat = as.data.frame(record.mat)
    colnames(record.mat) = c('epi.para.mean', 'epi.para.disp', 'epi.para.shift', 'll')

    max.ll = max(record.mat$ll)
    best.record.array = record.mat[which.max(record.mat$ll),]
    best.epi.para.mean = best.record.array$epi.para.mean[1]
    best.epi.para.disp = best.record.array$epi.para.disp[1]
    best.epi.para.shift = best.record.array$epi.para.shift[1]

    sel.record.mat = subset(record.mat, subset = record.mat$ll >= (max.ll - qchisq(p = 0.95, df = 1) / 1.5))
    # if(!is.null(can.epi.para.set)){
    #   sel.record.mat = sel.record.mat[, c(!c(isfixed.epi.para.mean, isfixed.epi.para.disp, isfixed.epi.para.shift), T)]
    #   return(sel.record.mat)
    # }

    #summary(sel.record.mat)
    epi.para.est.record = data.frame(
      epi.para.mean = c(best.epi.para.mean, range(sel.record.mat$epi.para.mean)),
      epi.para.disp = c(best.epi.para.disp, range(sel.record.mat$epi.para.disp)),
      epi.para.shift = c(best.epi.para.shift, range(sel.record.mat$epi.para.shift))
    )
    row.names(epi.para.est.record) = c('mle', 'ci.lwr', 'ci.upr')

    return(list(
      epi.para.est.output = epi.para.est.record[, !c(isfixed.epi.para.mean, isfixed.epi.para.disp, isfixed.epi.para.shift)],
      max.ll = max.ll,
      est.record.mat = record.mat[, c(!c(isfixed.epi.para.mean, isfixed.epi.para.disp, isfixed.epi.para.shift), T)]
    ))

  }

}

















#' To estimate model parameters using Markov chain Monte Carlo approach
#'
#' @description
#' This function (i.e., \code{paraest.MCMC()}) performs model parameter estimation using random walk **Markov chain Monte Carlo** (MCMC) approach with given structured contact tracing data.
#'
#' @param can.epi.para.start A list (\code{list}) of the starting values for unknown epidemiological parameters for offspring distribution.
#' The list should be in the format of \code{list(mean = ?, disp = ?, shift = ?)}.
#' Each parameter must be a scalar, and only accept non-negative values.
#' The default setting is given in the code Usage section.
#' For Delaporte distribution, the value of \code{mean} should be larger than the value of \code{shift}.
#' @param isfix.epi.para A list (\code{list}) of logistic values that indicate whether each parameter is fixed.
#' By default, \code{isfix.epi.para = list(mean = FALSE, disp = FALSE, shift = FALSE)}, where all three parameters are not fixed and to be estimated.
#' @param offspring.type
#' A character label (\code{character}) indicating the type of distribution used to describe the offspring distribution.
#' It only accepts one of the following values:
#' \itemize{
#'   \item{\code{"D"}}{ indicates the *Delaporte* distribution, }
#'   \item{\code{"NB"}}{ indicates the *negative binomial* distribution, }
#'   \item{\code{"G"}}{ indicates the *geometric* distribution, or }
#'   \item{\code{"P"}}{ indicates the *Poisson* distribution. }
#' }
#' By default, \code{offspring.type = 'D'}.
#' @param para.comb.num A positive integer for the number of iterations used for MCMC runs.
#' By default, \code{para.comb.num = 10000}, and no need to change the default setting here unless for special reasons.
#' @param burnin.frac A number with range strictly between 0 and 1 for the fraction of MCMC chain that will be discarded as the "burn-in" stage.
#' By default, \code{burnin.frac = 0.33}, and no need to change the default setting here unless for special reasons.
#' @param data A data frame (\code{data.frame}), or a vector (only when \code{obs.type.lab = "offspring"}) that contains the structured contact tracing data.
#' @param var.name A list (\code{list}), or a character of variable name for the column names of dataset given in \code{data}.
#' For a list of variable names, it should be in the format of \code{list(obssize = ?, seedsize = ?, typelab = ?)}.
#' Please see the details section for more information.
#' By default, \code{var.name = list(obssize = NULL, seedsize = NULL, typelab = NULL)}.
#' @param obs.type.lab A list (\code{list}), or a character of labels (i.e., "offspring", "nextgen", or "outbreak") for the type of observations.
#' For a list of labels, it should be in the format of \code{list(offspring = ?, nextgen = ?, outbreak = ?)}.
#' Please see the details section for more information.
#' By default, \code{obs.type.lab = list(offspring = NULL, nextgen = NULL, outbreak = NULL)}.
#'
#'
#' @return
#' A list (i.e., \code{list}) contains the following three items:
#' \itemize{
#'   \item{}{a data frame (\code{data.frame}) summaries the median, and 2.5% and 97.5% percentiles (95% credible interval, CrI) of the posterior MCMC samples for each unknown parameters,}
#'   \item{}{the maximum log-likelihood value, and}
#'   \item{}{a data frame (\code{data.frame}) of the posterior MCMC samples and their corresponding log-likelihood values.}
#' }
#'
#'
#' @details
#' For the values of parameters given in \code{can.epi.para.start},
#' they are some rough starting points for the MCMC to run, which are not necessarily to be precise (but have to be within a reasonable range), and
#' they will used as a "start" status to find the posterior MCMC samples.
#'
#' When \code{obs.type.lab} is a character, it should be either \code{"offspring"}, \code{"nextgen"}, or \code{"outbreak"} for type of observations.
#' When \code{obs.type.lab} is a list, this occurs when the contact tracing data has more than one types of observations.
#'
#' When the contact tracing dataset is offspring case observations, the function arguments \code{data} could be either a vector, or a data frame.
#' If \code{data} is a vector, it is not necessary to assign any value to \code{var.name}.
#' If \code{data} is a data frame, it is necessary to identify the variable name of offspring observations in \code{var.name}.
#'
#' When the contact tracing dataset is next-generation cluster size, or final outbreak size observations, the variable names of both observations and seed case size should be identified in \code{var.name} with the format of \code{list(obssize = ?, seedsize = ?)}.
#'
#' When the contact tracing dataset has more than one types of observations, the variable names of observations, seed case size, and observation type should be identified in \code{var.name} with the format of \code{list(obssize = ?, seedsize = ?, typelab = ?)}.
#'
#'
#'
#' @import Delaporte
#' @importFrom
#' stats rnorm runif quantile
#'
#' @export
#' paraest.MCMC
#'
#' @note
#' Each parameter in \code{can.epi.para.start = list(mean = ?, disp = ?, shift = ?)} should be a scalar, which means vector is not allowed here.
#'
#' For the contact tracing data in \code{data}, unknown observations (i.e., \code{NA}) is not allowed.
#'
#' As \code{para.comb.num} is the number of iterations for the MCMC runs, when \code{para.comb.num} is large, e.g., \code{para.comb.num} > 100000, the function \code{paraest.MCMC()} could take few seconds, or even minutes to complete, depending on the sample size, etc.
#' Thus, we do not recommend the users to change the default setting of \code{para.comb.num} unless for special reasons.
#'
#'
#' @seealso
#' \code{\link[modelSSE:overalllikelihood]{modelSSE::overalllikelihood()}}
#'
#'
#' @references
#' Blumberg S, Funk S, Pulliam JR. Detecting differential transmissibilities that affect the size of self-limited outbreaks. *PLoS Pathogens*. 2014;10(10):e1004452.
#' \doi{10.1371/journal.ppat.1004452}
#'
#' Kucharski AJ, Althaus CL. The role of superspreading in Middle East respiratory syndrome coronavirus (MERS-CoV) transmission. *Eurosurveillance*. 2015;20(25):21167.
#' \doi{10.2807/1560-7917.ES2015.20.25.21167}
#'
#' Adam DC, Wu P, Wong JY, Lau EH, Tsang TK, Cauchemez S, Leung GM, Cowling BJ. Clustering and superspreading potential of SARS-CoV-2 infections in Hong Kong. *Nature Medicine*. 2020;26(11):1714-1719.
#' \doi{10.1038/s41591-020-1092-0}
#'
#' Zhao S, Chong MK, Ryu S, Guo Z, He M, Chen B, Musa SS, Wang J, Wu Y, He D, Wang MH. Characterizing superspreading potential of infectious disease: Decomposition of individual transmissibility. *PLoS Computational Biology*. 2022;18(6):e1010281.
#' \doi{10.1371/journal.pcbi.1010281}
#'
#'
#' @examples
#' \donttest{
#'
#' # example 1: for offspring observations #
#' ## reproducing the parameter estimation results in Adam, et al. (2020)
#' ## paper doi link: https://doi.org/10.1038/s41591-020-1092-0,
#' ## (see the first row in Supplementary Table 4),
#' ## where R of 0.58 (95% CI: 0.45, 0.72), and k of 0.43 (95% CI: 0.29, 0.67).
#' data(COVID19_JanApr2020_HongKong)
#' set.seed(2023)
#' MCMC.output = paraest.MCMC(
#'   #can.epi.para.start = list(mean = 0.60, disp = 0.50, shift = 0.20),
#'   offspring.type = "NB", para.comb.num = 10000,
#'   data = COVID19_JanApr2020_HongKong$obs,
#'   obs.type.lab = 'offspring'
#' )
#' print(MCMC.output$epi.para.est.output)
#' ## Then, plot the posterior fitting results using 100 randomly-selected MCMC samples.
#' hist(
#'   COVID19_JanApr2020_HongKong$obs, breaks = c(0:100) -0.5, xlim = c(0,12),
#'   freq = FALSE, xlab = 'secondary cases', ylab = 'rel. freq.', main = ''
#' )
#' for(jj in 1:100){
#'   temp.random.col.index = sample.int(n = nrow(MCMC.output$est.record.mat), size = 1)
#'   est.record.array = MCMC.output$est.record.mat[temp.random.col.index,]
#'   lines(0:12, d_offspringdistn(
#'     x = 0:12,
#'     epi.para = list(
#'       mean = est.record.array$epi.para.mean,
#'       disp = est.record.array$epi.para.disp, shift = 0.1
#'     ),
#'     offspring.type = "NB"
#'   ), col = '#0066FF33', type = 'l', lty = 3)
#' }
#'
#'
#' # example 2: for offspring observations #
#' ## reproducing the parameter estimation results in Zhao, et al. (2020)
#' ## paper doi link: https://doi.org/10.1371/journal.pcbi.1010281,
#' ## (see the results of dataset #3 using Delaporte distribution in Table 1), where
#' ## R of 0.59 (95% CI: 0.46, 0.78),
#' ## k of 0.16 (95% CI: 0.06, 0.40), and
#' ## shift of 0.17 (95% CI: 0.04, 0.30).
#' data(COVID19_JanApr2020_HongKong)
#' set.seed(2023)
#' paraest.MCMC(
#'   can.epi.para.start = list(mean = 1, disp = 0.5, shift = 0.2),
#'   offspring.type = "D",
#'   data = COVID19_JanApr2020_HongKong$obs,
#'   obs.type.lab = 'offspring'
#' )$epi.para.est.output
#'
#'
#' # example 3: for next-generation cluster size observations #
#' ## reproducing the parameter estimation results in Blumberg, et al, (2014)
#' ## paper doi link: https://doi.org/10.1371/journal.ppat.1004452,
#' ## (see the last row in Table 3, and Fig 4A),
#' ## where R of 3.14 (95% CI: 2, >6), and k of 0.37 (95% CI: not reported).
#' data(smallpox_19581973_Europe)
#' set.seed(2023)
#' paraest.MCMC(
#'   can.epi.para.start = list(mean = 1, disp = 0.5, shift = 0.2),
#'   offspring.type = "NB",
#'   data = smallpox_19581973_Europe,
#'   var.name = list(obssize = 'obs.clustersize', seedsize = 'obs.seed'),
#'   obs.type.lab = 'nextgen'
#' )$epi.para.est.output
#'
#'
#' # example 4: final outbreak size observations #
#' ## reproducing the parameter estimation results in Kucharski, Althaus. (2015)
#' ## paper doi link: https://doi.org/10.2807/1560-7917.ES2015.20.25.21167,
#' ## (see Fig 1, and Finding section),
#' ## where R of 0.47 (95% CI: 0.29, 0.80), and k of 0.26 (95% CI: 0.09, 1.24).
#' data(MERS_2013_MEregion)
#' set.seed(2023)
#' paraest.MCMC(
#'   can.epi.para.start = list(mean = 1, disp = 0.5, shift = 0.2),
#'   offspring.type = "NB",
#'   data = MERS_2013_MEregion,
#'   var.name = list(obssize = 'obs.finalsize', seedsize = 'obs.seed'),
#'   obs.type.lab = 'outbreak'
#' )$epi.para.est.output
#'
#'
#' # example 5: for more than one types of observations #
#' ## reproducing the parameter estimation results in Blumberg, et al, (2014)
#' ## paper doi link: https://doi.org/10.1371/journal.ppat.1004452,
#' ## (see the last row in Table 5, and Fig 6A),
#' ## where R of 0.3 (95% CI: 0.2, 0.5), and k of 0.4 (95% CI: not reported).
#' data(mpox_19801984_DRC)
#' set.seed(2023)
#' paraest.MCMC(
#'   can.epi.para.start = list(mean = 1, disp = 0.5, shift = 0.2),
#'   offspring.type = "NB", para.comb.num = 30000,
#'   data = mpox_19801984_DRC,
#'   var.name = list(obssize = 'obs.size', seedsize = 'obs.seed', typelab = 'type'),
#'   obs.type.lab = list(offspring = 'offspring', nextgen = 'nextgen', outbreak = 'outbreak')
#' )$epi.para.est.output
#'
#' }
#'
paraest.MCMC = function(
  can.epi.para.start = list(mean = 1, disp = 0.5, shift = 0.2),
  isfix.epi.para = list(mean = FALSE, disp = FALSE, shift = FALSE),
  offspring.type = "D",
  para.comb.num = 10000, burnin.frac = 0.33,
  data = NULL,
  var.name = list(obssize = NULL, seedsize = NULL, typelab = NULL),
  obs.type.lab = list(offspring = NULL, nextgen = NULL, outbreak = NULL)
){
  if(burnin.frac <= 0 | burnin.frac >= 1){
    stop('for the function argument  \"burnin.frac\" , it should be strictly between 0 and 1, and a range from 0.1 to 0.5 is recommended.')
  }
  mcmc.runs = para.comb.num
  now.ll = -Inf

  now.epi.para.mean = unlist(can.epi.para.start$mean)[1]
  now.epi.para.disp = unlist(can.epi.para.start$disp)[1]
  now.epi.para.shift = unlist(can.epi.para.start$shift)[1]

  isfixed.epi.para.mean = unlist(isfix.epi.para$mean)[1]
  if(offspring.type %in% c('G', 'P')){
    isfixed.epi.para.disp = TRUE
    now.epi.para.mean = 1
  } else {
    isfixed.epi.para.disp = unlist(isfix.epi.para$disp)[1]
  }
  if(offspring.type %in% c("NB", 'G', 'P')){
    isfixed.epi.para.shift = TRUE
    now.epi.para.shift = 1e-09
  } else {
    isfixed.epi.para.shift = unlist(isfix.epi.para$shift)[1]
  }

  single.record.mat = NULL
  for(run.i in 1:mcmc.runs){  #           run.i = run.i +1
    cooling.term = as.integer((run.i / mcmc.runs)*7 +2) #  2 ~ 9
    disturb.factor = ifelse((run.i %% cooling.term) == 0, 0.1, 0.01)

    if(isfixed.epi.para.mean){
      can.epi.para.mean = now.epi.para.mean
    } else {
      can.epi.para.mean = exp(c(rnorm(n = 1, mean = log(now.epi.para.mean), sd = 1*disturb.factor)))
    }
    #
    if(isfixed.epi.para.disp){
      can.epi.para.disp = now.epi.para.mean
    } else {
      can.epi.para.disp = exp(c(rnorm(n = 1, mean = log(now.epi.para.disp), sd = 1*disturb.factor)))
    }
    #
    if(isfixed.epi.para.shift){
      can.epi.para.shift = now.epi.para.shift
    } else {
      can.epi.para.shift = exp(c(rnorm(n = 1, mean = log(now.epi.para.shift), sd = 1*disturb.factor)))
    }

    can.sampling.ll = 1; now.sampling.ll = 1

    can.ll = overalllikelihood(
      epi.para = list(mean = can.epi.para.mean, disp = can.epi.para.disp, shift = can.epi.para.shift),
      offspring.type = offspring.type,
      data = data,
      var.name = var.name,
      obs.type.lab = obs.type.lab
    )

    uptake.prob = exp(can.ll - now.ll + now.sampling.ll - can.sampling.ll); uptake.prob = ifelse(uptake.prob >1, 1, uptake.prob)
    random.prob = runif(1); is.update = (random.prob <= uptake.prob)
    if(is.update){
      now.epi.para.mean = can.epi.para.mean
      now.epi.para.disp = can.epi.para.disp
      now.epi.para.shift = can.epi.para.shift
      now.ll = can.ll
    }

    single.record.array = c(now.epi.para.mean, now.epi.para.disp, now.epi.para.shift, now.ll)
    single.record.mat = rbind(single.record.mat, c(single.record.array))
    #print(paste0(round(run.i /mcmc.runs*100, 2), '%, ', ifelse(is.update, 'updated', 'unchanged'), '.'))
  }
  single.record.mat = as.data.frame(single.record.mat)
  colnames(single.record.mat) = c(
    'epi.para.mean', 'epi.para.disp', 'epi.para.shift',
    'll'
  )
  #    plot(single.record.mat$epi.para.mean, type = 'l')
  mature.record.mat = single.record.mat[-c(1:ceiling(nrow(single.record.mat) *burnin.frac)),]
  mature.record.mat = as.data.frame(mature.record.mat)


  max.ll = max(mature.record.mat$ll)
  # best.record.array = mature.record.mat[which.max(mature.record.mat$ll),]
  # best.epi.para.mean = best.record.array$epi.para.mean[1]
  # best.epi.para.disp = best.record.array$epi.para.disp[1]
  # best.epi.para.shift = best.record.array$epi.para.shift[1]

  epi.para.est.record = data.frame(
    epi.para.mean = quantile(mature.record.mat$epi.para.mean, c(0.500, 0.025, 0.975)),
    epi.para.disp = quantile(mature.record.mat$epi.para.disp, c(0.500, 0.025, 0.975)),
    epi.para.shift = quantile(mature.record.mat$epi.para.shift, c(0.500, 0.025, 0.975))
  )
  row.names(epi.para.est.record) = c('med.est', 'cri.lwr', 'cri.upr')

  return(list(
    epi.para.est.output = epi.para.est.record[, !c(isfixed.epi.para.mean, isfixed.epi.para.disp, isfixed.epi.para.shift)],
    max.ll = max.ll,
    est.record.mat = mature.record.mat[, c(!c(isfixed.epi.para.mean, isfixed.epi.para.disp, isfixed.epi.para.shift), T)]
  ))

}



















