



#' To convert model parameters
#'
#' @description
#' This function (i.e., \code{convert.epipara.to.delappara()}) converts a list of pre-defined epidemiological parameters into a list of statistical parameters,
#' such that a probability function of Delaporte distribution may recognize, which is for **internal use only**.
#'
#' @param para A list (\code{list}) of epidemiological parameters for offspring distribution, in the format of \code{list(mean = ?, disp = ?, shift = ?)},
#' where the three parameters accept non-negative values.
#' Each parameter can be either a scalar, or a vector.
#' For the parameters being assigned with values of vectors, the vectors should be of the same length.
#' For Delaporte distribution, the value of \code{mean} should be larger than the value of \code{shift}.
#' @param offspring.type
#' A character label (\code{character}) indicating the type of distribution used to describe the offspring distribution.
#' It only accepts one of the following values:
#' \itemize{
#'   \item{\code{"D"}}{ indicates the *Delaporte* distribution, }
#'   \item{\code{"NB"}}{ indicates the *negative binomial* distribution, }
#'   \item{\code{"G"}}{ indicates the *geometric* distribution, or }
#'   \item{\code{"P"}}{ indicates the *Poisson* distribution.}
#' }
#' By default, \code{offspring.type = 'D'}.
#'
#' @details
#' For different values of \code{offspring.type},
#' \itemize{
#'   \item{When \code{offspring.type = "D"},}{ no action to the pre-defined epidemiological parameters; }
#'   \item{When \code{offspring.type = "NB"},}{ we set parameter \code{shift = 0} internally; }
#'   \item{When \code{offspring.type = "G"},}{ we set parameters \code{disp = 1} and \code{shift = 0} internally; and }
#'   \item{When \code{offspring.type = "P"},}{ we set parameters \code{disp = +Inf} and \code{shift = mean} internally. }
#' }
#'
#' @return
#' A list of statistical parameters in the format of \code{list(alpha = ?, beta = ?, lambda = ?)}
#' that can be recognized by Delaporte distribution.
#'
#' @export
#'
#'
#'
#' @note
#' It would be difficult to interpret the converted statistical parameters,
#' and thus we set this function \code{convert.epipara.to.delappara()} as a internal function.
#' We do not recommend the users to use this function externally unless for special reasons.
#'
#'
#' @references
#' Vose D. Risk analysis: a quantitative guide. John Wiley & Sons. 2008; pp. 618-619. ISBN: 978-0-470-51284-5
#'
#'
#' @seealso
#' \code{\link[Delaporte:ddelap]{Delaporte}} for the parameterization of Delaporte distribution.
#'
#' @examples
#' convert.epipara.to.delappara(
#'   para = list(mean = 1, disp = 0.5, shift = 0.2),
#'   offspring.type = 'D'
#' )
#' convert.epipara.to.delappara(
#'   para = list(mean = 1, disp = 1, shift = 0),
#'   offspring.type = 'G'
#' )
#'
#' convert.epipara.to.delappara(
#'   para = list(mean = c(0.5, 1, 2), disp = c(0.1, 0.1, 0.3), shift = 0.2)
#' )
#'
#' @keywords internal
#'
#@keywords internal
#@noRd
convert.epipara.to.delappara = function(
  para = list(mean = 1, disp = 0.5, shift = 0.2),
  offspring.type = 'D'
){

  para.mean = para$mean
  #     para.mean = c(0.5, 1, 2)
  if(!is.numeric(para.mean)){
    stop("the value of function argument  \"mean\"  must be numeric.")
  }

  if(offspring.type == 'D'){
    para.disp = para$disp; para.shift = para$shift
  } else if (offspring.type == 'P'){
    #para.shift = para.mean - 1e-99
    para.shift = 1e-9
    #para.shift = para$shift
    para.disp = 1e+9
  } else {
    para.shift = 1e-9
    if(offspring.type == 'NB'){
      para.disp = para$disp
    } else if(offspring.type == 'G'){
      para.disp = 1
    } else {
      stop(
        "the value of the function argument  \"offspring.type\"  must be set as one of following four labels, including:
          -- \"D\" for Delaporte distribution,
          -- \"NB\" for negative binomial distribution,
          -- \"G\" for geometric distribution, and
          -- \"P\" for Poisson distribution."
      )
    }
  }

  variable.R = para.mean - para.shift
  delap.para.alpha = para.disp; delap.para.alpha = ifelse(delap.para.alpha == 0, 1e-99, delap.para.alpha)
  delap.para.beta = variable.R / para.disp; delap.para.beta = ifelse(delap.para.beta == 0, 1e-99, delap.para.beta)
  delap.para.lambda = para.shift; delap.para.lambda = ifelse(delap.para.lambda == 0, 1e-99, delap.para.lambda)

  para = list(
    alpha = delap.para.alpha, beta = delap.para.beta, lambda = delap.para.lambda
  )
  return(para)
}

















#' The distribution of individual reproduction number
#'
#' @description
#' This function (i.e., \code{d_reproductiondistn()}) is the probability density function (PDF) of **individual reproduction number** that was modelled as a shifted gamma distribution.
#'
#' @param x A scalar, or a vector of non-negative integer.
#' @param epi.para A list (\code{list}) of pre-defined epidemiological parameters for offspring distribution, in the format of \code{list(mean = ?, disp = ?, shift = ?)},
#' where the three parameters accept non-negative values.
#' Each parameter can be either a scalar, or a vector.
#' For the parameters being assigned with values of vectors, the vectors should be of the same length.
#' For Delaporte distribution, the value of \code{mean} should be larger than the value of \code{shift}.
#' @param offspring.type
#' A character label (\code{character}) indicating the type of distribution used to describe the offspring distribution.
#' It only accepts one of the following values:
#' \itemize{
#'   \item{\code{"D"}}{ indicates the Delaporte distribution for offspring cases, where reproduction number follows a *shifted gamma* distribution; }
#'   \item{\code{"NB"}}{ indicates the negative binomial distribution for offspring cases, where reproduction number follows a (non-shifted, or standard) *gamma* distribution; }
#'   \item{\code{"G"}}{ indicates the geometric distribution for offspring cases, where reproduction number follows an *exponential* distribution; or }
#'   \item{\code{"P"}}{ indicates the Poisson distribution for offspring cases, where reproduction number follows a *Dirac delta* distribution. }
#' }
#' By default, \code{offspring.type = 'D'}.
#' @param is.log A logical variable, under which probability would be taken natural logarithm, if \code{is.log = TRUE}.
#' By default, \code{is.log = FALSE}.
#'
#'
#' @return
#' \code{d_reproductiondistn()} is the probability density function (PDF), and it returns value of probability density (non-negative value).
#'
#' @export
#' d_reproductiondistn
#'
#' @importFrom
#' stats dgamma
#'
#' @note
#' Only the PDF of individual reproduction number (i.e., \code{d_reproductiondistn()}) was created here (without cumulative distribution, quantile, or random variable generating functions).
#' The function \code{d_reproductiondistn()} was used mainly for data visualization purpose,
#' because the distribution of individual reproduction number was not explicitly used in model fitting to the disease contact tracing data.
#'
#' When \code{offspring.type = "P"}, individual reproduction number follows a *Dirac delta* distribution, which is difficult to return any value, or visualize the PDF, because of the nature of this pulse function.
#'
#' @references
#' Lloyd-Smith JO, Schreiber SJ, Kopp PE, Getz WM. Superspreading and the effect of individual variation on disease emergence. *Nature*. 2005;438(7066):355-359.
#' \doi{10.1038/nature04153}
#'
#' Zhao S, Chong MK, Ryu S, Guo Z, He M, Chen B, Musa SS, Wang J, Wu Y, He D, Wang MH. Characterizing superspreading potential of infectious disease: Decomposition of individual transmissibility. *PLoS Computational Biology*. 2022;18(6):e1010281.
#' \doi{10.1371/journal.pcbi.1010281}
#'
#'
#' @examples
#'
#' ## an example to visualize individual reproduction number is as follows.
#' plot(seq(0.01,9.99, length.out = 1001), d_reproductiondistn(
#'   x = seq(0,10, length.out = 1001),
#'   epi.para = list(mean = 2, disp = 1.5, shift = 0.5),
#'   offspring.type = "D",
#'   is.log = FALSE
#' ), type = 'l', xlab = 'individual reproduction number', ylab = 'density')
#'
d_reproductiondistn = function(
  x = 1,
  epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
  offspring.type = 'D', is.log = FALSE
){
  delap.para = convert.epipara.to.delappara(para = epi.para, offspring.type = offspring.type)

  PDF.array = dgamma(
    x = x - delap.para$lambda, log = is.log,
    shape = delap.para$alpha, scale = delap.para$beta#, lambda = delap.para$lambda
  )

  return(c(PDF.array))
}

















#' The offspring distribution
#'
#' @description
#' Density, cumulative distribution, quantile, and random variable generating functions for the **offspring** distribution with pre-defined epidemiological parameters.
#'
#' @param x A scalar, or a vector of non-negative integer.
#' @param q A scalar, or a vector of non-negative number (not necessarily integer).
#' @param p A scalar, or a vector of probability (i.e., ranging from 0 to 1).
#' @param n A scalar of positive integer.
#' @param epi.para A list (\code{list}) of pre-defined epidemiological parameters for offspring distribution, in the format of \code{list(mean = ?, disp = ?, shift = ?)},
#' where the three parameters accept non-negative values.
#' Each parameter can be either a scalar, or a vector.
#' For the parameters being assigned with values of vectors, the vectors should be of the same length.
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
#' @param is.log A logical variable, under which probability would be taken natural logarithm, if \code{is.log = TRUE}.
#' By default, \code{is.log = FALSE}.
#' @param lower.tail A logical variable, under which the probability is cumulative distribution function (CDF, i.e., *P(X <= x)*), if \code{is.log = TRUE}, and otherwise, 1 - CDF (i.e., *P(X > x)*).
#' By default, \code{lower.tail = TRUE}.
#'
#' @details
#' For different values of \code{offspring.type},
#' \itemize{
#'   \item{When \code{offspring.type = "D"},}{ no action to any parameter; }
#'   \item{When \code{offspring.type = "NB"},}{ we set parameter \code{shift = 0}; }
#'   \item{When \code{offspring.type = "G"},}{ we set parameters \code{disp = 1} and \code{shift = 0}; and }
#'   \item{When \code{offspring.type = "P"},}{ we set parameters \code{disp = +Inf} and \code{shift = mean}. }
#' }
#'
#' @return
#' For the values returned from the four functions,
#' \itemize{
#'   \item{\code{d_offspringdistn()}}{ is the *probability mass function* (PMF), and it returns value of probability (i.e., ranging from 0 to 1); }
#'   \item{\code{p_offspringdistn()}}{ is the *cumulative distribution function* (CDF), and it returns value of probability (i.e., ranging from 0 to 1); }
#'   \item{\code{q_offspringdistn()}}{ is the *quantile function*, and it returns value of quantile (non-negative integer); and }
#'   \item{\code{r_offspringdistn()}}{ is the *random variable generating function*, and it generates a set of random variables (non-negative integers) with size given in \code{n}. }
#' }
#'
#' @import Delaporte
#'
#' @export
#' d_offspringdistn p_offspringdistn q_offspringdistn r_offspringdistn
#'
#' @note
#' Depending on the values of parameters, the functions could take hours to complete, given the double-summation nature for the Delaporte distribution.
#'
#'
#' @references
#' Vose D. Risk analysis: a quantitative guide. John Wiley & Sons. 2008; pp. 618-619. ISBN: 978-0-470-51284-5
#'
#' Lloyd-Smith JO, Schreiber SJ, Kopp PE, Getz WM. Superspreading and the effect of individual variation on disease emergence. *Nature*. 2005;438(7066):355-359.
#' \doi{10.1038/nature04153}
#'
#' Zhao S, Chong MK, Ryu S, Guo Z, He M, Chen B, Musa SS, Wang J, Wu Y, He D, Wang MH. Characterizing superspreading potential of infectious disease: Decomposition of individual transmissibility. *PLoS Computational Biology*. 2022;18(6):e1010281.
#' \doi{10.1371/journal.pcbi.1010281}
#'
#'
#' @seealso
#' \code{\link[Delaporte:ddelap]{Delaporte}} for the parameterization of Delaporte distribution.
#'
#' @examples
#' ## Please see the "Usage" section.
#'
#'
#' ## the following returns the proportion of index cases that generated at least 1 offspring cases.
#' p_offspringdistn(
#'   q = 0,
#'   epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
#'   offspring.type = 'D', lower.tail = FALSE
#' )
#'
#' ## reproducing the results in Adam, et al. (2020)
#' ## paper doi link: https://doi.org/10.1038/s41591-020-1092-0 (see Fig 3b),
#' ## where the number of offspring cases were fitted
#' ## with parameter R of 0.58 and k of 0.43 under NB distribution.
#' data(COVID19_JanApr2020_HongKong)
#' hist(
#'   COVID19_JanApr2020_HongKong$obs, breaks = c(0:100) -0.5, xlim = c(0,12),
#'   freq = FALSE, xlab = 'secondary cases', ylab = 'rel. freq.', main = ''
#' )
#' lines(0:12, d_offspringdistn(
#'   x = 0:12,
#'   epi.para = list(mean = 0.58, disp = 0.43, shift = 0.2),
#'   offspring.type = "NB"
#' ), pch = 20, type = 'o', lty = 2)
#'
#'
#' ## an example to generate 100 rv of offspring case number
#' table(r_offspringdistn(
#'   n = 100,
#'   epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
#'   offspring.type = 'D'
#' ))
#'
d_offspringdistn = function(
  x = 1,
  epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
  offspring.type = 'D', is.log = FALSE
){
  delap.para = convert.epipara.to.delappara(para = epi.para, offspring.type = offspring.type)

  PMF.array = Delaporte::ddelap(
    x = x, log = is.log,
    alpha = delap.para$alpha, beta = delap.para$beta, lambda = delap.para$lambda
  )

  return(c(PMF.array))
}

#'
#' @rdname d_offspringdistn
p_offspringdistn = function(
  q = 1.5,
  epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
  offspring.type = 'D', is.log = FALSE, lower.tail = TRUE
){
  delap.para = convert.epipara.to.delappara(para = epi.para, offspring.type = offspring.type)

  CDF.array = Delaporte::pdelap(
    q = q, log.p = is.log, lower.tail = lower.tail,
    alpha = delap.para$alpha, beta = delap.para$beta, lambda = delap.para$lambda
  )

  return(c(CDF.array))
}

#'
#' @rdname d_offspringdistn
q_offspringdistn = function(
  p = 0.8,
  epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
  offspring.type = 'D', lower.tail = TRUE#, is.log = FALSE
){
  delap.para = convert.epipara.to.delappara(para = epi.para, offspring.type = offspring.type)

  q.array = Delaporte::qdelap(
    p = p, log.p = FALSE, lower.tail = lower.tail, exact = TRUE,
    alpha = delap.para$alpha, beta = delap.para$beta, lambda = delap.para$lambda
  )

  return(c(q.array))
}

#'
#' @rdname d_offspringdistn
r_offspringdistn = function(
  n = 10,
  epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
  offspring.type = 'D'
){
  delap.para = convert.epipara.to.delappara(para = epi.para, offspring.type = offspring.type)

  rv.array = Delaporte::rdelap(
    n = n, exact = TRUE,
    alpha = delap.para$alpha, beta = delap.para$beta, lambda = delap.para$lambda
  )

  return(c(rv.array))
}























#' The next-generation cluster size distribution
#'
#' @description
#' Density, cumulative distribution, quantile, and random variable generating functions for the **next-generation cluster size** distribution with pre-defined epidemiological parameters.
#'
#' @param x A scalar, or a vector of positive integer, for the next-generation cluster size.  The value of \code{x} must be not less than \code{seed.size}.
#' @param q A scalar, or a vector of positive number (not necessarily integer), for the next-generation cluster size. The value of \code{q} must be not less than \code{seed.size}.
#' @param p A scalar, or a vector of probability (i.e., ranging from 0 to 1).
#' @param n A scalar of positive integer.
#' @param seed.size A scalar, or a vector of positive integer.
#' For vector type of \code{seed.size}, it only applies to \code{d_nextgenclusterdistn()}, \code{p_nextgenclusterdistn()}, and \code{q_nextgenclusterdistn()}.
#' If \code{seed.size} and \code{x}, \code{q} or \code{p} are vectors, \code{seed.size} should be of the same length as \code{x}, \code{q} or \code{p}.
#' @param epi.para A list (\code{list}) of pre-defined epidemiological parameters for offspring distribution, in the format of \code{list(mean = ?, disp = ?, shift = ?)},
#' where the three parameters accept non-negative values.
#' Each parameter can be either a scalar, or a vector.
#' For the parameters being assigned with values of vectors, the vectors should be of the same length.
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
#' @param is.log A logical variable, under which probability would be taken natural logarithm, if \code{is.log = TRUE}.
#' By default, \code{is.log = FALSE}.
#' @param lower.tail A logical variable, under which the probability is cumulative distribution function (CDF, i.e., *P(X <= x)*), if \code{is.log = TRUE}, and otherwise, 1 - CDF (i.e., *P(X > x)*).
#' By default, \code{lower.tail = TRUE}.
#'
#' @details
#' Function \code{d_nextgenclusterdistn()} returns the probability of having a next-generation case cluster with size \code{x} generated by \code{seed.size} index cases, where (\code{seed.size}) is given.
#'
#' Function \code{p_nextgenclusterdistn()} returns the probability of having a next-generation case cluster with size less than or equal to, or larger than \code{q} (depending on the value of \code{lower.tail}), generated by \code{seed.size} index cases, where (\code{seed.size}) is given.
#'
#' Function \code{q_nextgenclusterdistn()} returns a value such that there is a probability of \code{p} for having a next-generation case cluster with size less than or equal to, or larger than this value (depending on the value of \code{lower.tail}) generated by \code{seed.size} index cases, where (\code{seed.size}) is given.
#'
#' Function \code{r_nextgenclusterdistn()} returns a set of random variables of next-generation cluster size given in \code{n}, given (\code{seed.size}).
#'
#'
#' @return
#' For the values returned from the four functions,
#' \itemize{
#'   \item{\code{d_nextgenclusterdistn()}}{ is the *probability mass function* (PMF), and it returns value of probability (i.e., ranging from 0 to 1); }
#'   \item{\code{p_nextgenclusterdistn()}}{ is the *cumulative distribution function* (CDF), and it returns value of probability (i.e., ranging from 0 to 1); }
#'   \item{\code{q_nextgenclusterdistn()}}{ is the *quantile function*, and it returns value of quantile (positive integer); and }
#'   \item{\code{r_nextgenclusterdistn()}}{ is the *random variable generating function*, and it generates a set of random variables (positive integers) with size given in \code{n}. }
#' }
#'
#' @import Delaporte
#'
#' @export
#' d_nextgenclusterdistn p_nextgenclusterdistn q_nextgenclusterdistn r_nextgenclusterdistn
#'
#' @note
#' Depending on the values of parameters, the functions could take hours to complete, given the double-summation nature for the Delaporte distribution.
#'
#'
#' @references
#' Blumberg S, Lloyd-Smith JO. Inference of R 0 and transmission heterogeneity from the size distribution of stuttering chains. *PLoS Computational Biology*. 2013 May 2;9(5):e1002993.
#' \doi{10.1371/journal.pcbi.1002993}
#'
#' Zhao S, Chong MK, Ryu S, Guo Z, He M, Chen B, Musa SS, Wang J, Wu Y, He D, Wang MH. Characterizing superspreading potential of infectious disease: Decomposition of individual transmissibility. *PLoS Computational Biology*. 2022;18(6):e1010281.
#' \doi{10.1371/journal.pcbi.1010281}
#'
#'
#' @seealso
#' \code{\link[Delaporte:ddelap]{Delaporte}} for the parameterization of Delaporte distribution.
#'
#' @examples
#' ## Please see the "Usage" section.
#'
d_nextgenclusterdistn = function(
  x = 5, seed.size = 1,
  epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
  offspring.type = 'D', is.log = FALSE
){
  if(sum(seed.size < 1) > 0){
    stop('the value of the function argument  \"seed.size\"  must NOT be less than  1.')
  }
  if(sum(x < seed.size) > 0){
    stop('the value of the function argument  \"x\"  must NOT be less than the value of  \"seed.size\" .')
  }

  epi.para$mean = epi.para$mean *seed.size
  epi.para$disp = epi.para$disp *seed.size
  epi.para$shift = epi.para$shift *seed.size
  delap.para = convert.epipara.to.delappara(para = epi.para, offspring.type = offspring.type)

  PMF.array = Delaporte::ddelap(
    x = (x - seed.size), log = is.log,
    alpha = delap.para$alpha, beta = delap.para$beta, lambda = delap.para$lambda
  )

  return(c(PMF.array))
}

#'
#' @rdname d_nextgenclusterdistn
p_nextgenclusterdistn = function(
  q = 10.5, seed.size = 1,
  epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
  offspring.type = 'D', lower.tail = TRUE, is.log = FALSE
){
  if(sum(seed.size < 1) > 0){
    stop('the value of the function argument  \"seed.size\"  must NOT be less than  1.')
  }
  if(sum(q < seed.size) > 0){
    stop('the value of the function argument  \"q\"  must NOT be less than the value of  \"seed.size\" .')
  }

  epi.para$mean = epi.para$mean *seed.size
  epi.para$disp = epi.para$disp *seed.size
  epi.para$shift = epi.para$shift *seed.size
  delap.para = convert.epipara.to.delappara(para = epi.para, offspring.type = offspring.type)

  CDF.array = Delaporte::pdelap(
    q = (q - seed.size), log.p = is.log, lower.tail = lower.tail,
    alpha = delap.para$alpha, beta = delap.para$beta, lambda = delap.para$lambda
  )

  return(c(CDF.array))
}

#'
#' @rdname d_nextgenclusterdistn
q_nextgenclusterdistn = function(
  p = 0.8, seed.size = 1,
  epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
  offspring.type = 'D', lower.tail = TRUE
){
  if(seed.size < 1){
    stop('the value of the function argument  \"seed.size\"  must NOT be less than  1.')
  }

  epi.para$mean = epi.para$mean *seed.size
  epi.para$disp = epi.para$disp *seed.size
  epi.para$shift = epi.para$shift *seed.size
  delap.para = convert.epipara.to.delappara(para = epi.para, offspring.type = offspring.type)

  q.array = Delaporte::qdelap(
    p = p, log.p = FALSE, lower.tail = lower.tail, exact = TRUE,
    alpha = delap.para$alpha, beta = delap.para$beta, lambda = delap.para$lambda
  ) + seed.size

  return(c(q.array))
}

#'
#' @rdname d_nextgenclusterdistn
r_nextgenclusterdistn = function(
  n = 10, seed.size = 1,
  epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
  offspring.type = 'D'
){
  if(seed.size < 1){
    stop('the value of the function argument  \"seed.size\"  must NOT be less than  1.')
  }

  epi.para$mean = epi.para$mean *seed.size
  epi.para$disp = epi.para$disp *seed.size
  epi.para$shift = epi.para$shift *seed.size
  delap.para = convert.epipara.to.delappara(para = epi.para, offspring.type = offspring.type)

  rv.array = Delaporte::rdelap(
    n = n, exact = TRUE,
    alpha = delap.para$alpha, beta = delap.para$beta, lambda = delap.para$lambda
  ) + seed.size

  return(c(rv.array))
}




















#' The final outbreak size distribution
#'
#' @description
#' Density, cumulative distribution, quantile, and random variable generating functions for the **final outbreak size** distribution with pre-defined epidemiological parameters.
#'
#' @param x A scalar, or a vector of final outbreak size, which is positive integer. The value of \code{x} must be not less than \code{seed.size}.
#' @param q A scalar, or a vector of positive number (not necessarily integer). The value of \code{q} must be not less than \code{seed.size}.
#' @param p A scalar, or a vector of probability (i.e., ranging from 0 to 1).
#' @param n A scalar of positive integer.
#' @param seed.size A scalar, or a vector of positive integer.
#' For vector type of \code{seed.size}, it only applies to \code{d_outbreakdistn()}, and \code{p_outbreakdistn()}.
#' If \code{seed.size} and \code{x} or \code{q} are vectors, \code{seed.size} should be of the same length as \code{x} or \code{q}.
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
#' @param is.log A logical variable, under which probability would be taken natural logarithm, if \code{is.log = TRUE}.
#' By default, \code{is.log = FALSE}.
#' @param lower.tail A logical variable, under which the probability is cumulative distribution function (CDF, i.e., *P(X <= x)*), if \code{is.log = TRUE}, and otherwise, 1 - CDF (i.e., *P(X > x)*).
#' By default, \code{lower.tail = TRUE}.
#' @param upr.limit A positive integer.
#' If the result final outbreak size is larger than \code{upr.limit}, the returned value will be set as "\code{> upr.limit}".
#' The value of \code{upr.limit} must be larger than \code{seed.size}.
#' By default, \code{upr.limit = 1000}, and no need to change the default setting here unless for special reasons.
#'
#'
#' @return
#' For the values returned from the four functions,
#' \itemize{
#'   \item{\code{d_outbreakdistn()}}{ is the *probability mass function* (PMF), and it returns value of probability (i.e., ranging from 0 to 1); }
#'   \item{\code{p_outbreakdistn()}}{ is the *cumulative distribution function* (CDF), and it returns value of probability (i.e., ranging from 0 to 1); }
#'   \item{\code{q_outbreakdistn()}}{ is the *quantile function*, and it returns value of quantile (positive integer); and }
#'   \item{\code{r_outbreakdistn()}}{ is the *random variable generating function*, and it generates a set of random variables (positive integers) with size given in \code{n}. }
#' }
#'
#' Specially, due to the computational consumption, an upper limit was set for the functions \code{q_outbreakdistn()} and \code{r_outbreakdistn()},
#' i.e., \code{upr.limit}, and thus, both functions here return value in string form (i.e., \code{character}).
#'
#' @details
#' Function \code{d_outbreakdistn()} returns the probability of having an outbreak with final size \code{x} generated by \code{seed.size} index cases, where (\code{seed.size}) is given.
#'
#' Function \code{p_outbreakdistn()} returns the probability of having an outbreak with final size less than or equal to, or larger than \code{q} (depending on the value of \code{lower.tail}), generated by \code{seed.size} index cases, where (\code{seed.size}) is given.
#'
#' Function \code{q_outbreakdistn()} returns a value such that there is a probability of \code{p} for having a final outbreak size less than or equal to, or larger than this value (depending on the value of \code{lower.tail}) generated by \code{seed.size} index cases, where (\code{seed.size}) is given.
#'
#' Function \code{r_outbreakdistn()} returns a set of random variables of final outbreak size given in \code{n}, given (\code{seed.size}).
#'
#'
#' @import Delaporte
#' @importFrom
#' stats runif
#'
#' @export
#' d_outbreakdistn p_outbreakdistn q_outbreakdistn r_outbreakdistn
#'
#' @note
#' Each parameter in \code{epi.para = list(mean = ?, disp = ?, shift = ?)} should be a scalar, which means vector is not allowed here.
#'
#' When \code{q} is large, e.g., \code{q} > 10000, the function \code{p_outbreakdistn()} could take few seconds, or even minutes to complete.
#'
#' When \code{upr.limit} is large, e.g., \code{upr.limit} > 10000, the functions \code{q_outbreakdistn()} and \code{r_outbreakdistn()} could take few seconds, or even minutes to complete.
#' Thus, we do not recommend the users to change the default setting of \code{upr.limit} unless for special reasons.
#'
#' @seealso
#' \code{\link[modelSSE:d_offspringdistn]{d_offspringdistn}}
#'
#' @references
#' Farrington CP, Kanaan MN, Gay NJ. Branching process models for surveillance of infectious diseases controlled by mass vaccination. *Biostatistics*. 2003;4(2):279-95.
#' \doi{10.1093/biostatistics/4.2.279}
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
#' Zhao S, Chong MK, Ryu S, Guo Z, He M, Chen B, Musa SS, Wang J, Wu Y, He D, Wang MH. Characterizing superspreading potential of infectious disease: Decomposition of individual transmissibility. *PLoS Computational Biology*. 2022;18(6):e1010281.
#' \doi{10.1371/journal.pcbi.1010281}
#'
#'
#' @examples
#'
#' ## an example to generate 1000 rv of final outbreak size
#' table(r_outbreakdistn(
#'   n = 1000,
#'   seed.size = 1,
#'   epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
#'   offspring.type = "D",
#'   upr.limit = 10
#' ))
#'
#'
#' \donttest{
#'
#' ## an attempt to reproduce the results in Guo, et al. (2022)
#' ## paper doi link: https://doi.org/10.1016/j.jinf.2022.05.041 (see Fig 1B),
#' ## where the probability of one seed case generating an outbreak with final size >= a given number,
#' ## with parameter R of 0.78 and k of 0.10 under NB distribution.
#' plot(1:100, 1 - c(0,cumsum(d_outbreakdistn(
#'   x = 1:99,
#'   seed.size = 1,
#'   epi.para = list(mean = 0.78, disp = 0.10, shift = 0.2),
#'   offspring.type = "NB",
#' ))), log = 'y', type = 'l', xlab = 'outbreak size', ylab = 'probability')
#' plot(1:100, c(1,p_outbreakdistn(
#'   q = 1:99,
#'   seed.size = 1,
#'   epi.para = list(mean = 0.78, disp = 0.10, shift = 0.2),
#'   offspring.type = "NB",
#'   lower.tail = FALSE
#' )), log = 'y', type = 'l', xlab = 'outbreak size', ylab = 'probability')
#'
#' }
#'
d_outbreakdistn = function(
  x = 10, seed.size = 1,
  epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
  offspring.type = 'D', is.log = FALSE
){
  if(sum(seed.size < 1) > 0){
    stop('the value of the function argument  \"seed.size\"  must NOT be less than  1.')
  }
  if(sum(x < seed.size) > 0){
    stop('the value of the function argument  \"x\"  must NOT be less than the value of  \"seed.size\" .')
  }

  epi.para$mean = epi.para$mean *x
  epi.para$disp = epi.para$disp *x
  epi.para$shift = epi.para$shift *x
  delap.para = convert.epipara.to.delappara(para = epi.para, offspring.type = offspring.type)
  PMF.array = Delaporte::ddelap(
    x = (x - seed.size), log = is.log,
    alpha = delap.para$alpha, beta = delap.para$beta, lambda = delap.para$lambda
  )

  if(is.log){
    PMF.array = PMF.array + log(seed.size / x)
  } else {
    PMF.array = PMF.array * (seed.size / x)
  }

  return(c(PMF.array))
}

#'
#' @rdname d_outbreakdistn
p_outbreakdistn = function(
  q = 30.5, seed.size = 1,
  epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
  offspring.type = 'D', lower.tail = TRUE, is.log = FALSE
){
  if(sum(seed.size < 1) > 0){
    stop('the value of the function argument  \"seed.size\"  must NOT be less than  1.')
  }
  if(sum(q < seed.size) > 0){
    stop('the value of the function argument  \"q\"  must NOT be less than the value of  \"seed.size\" .')
  }

  CDF.array = NULL
  for (item.j in 1:max(length(q),length(seed.size))) {
    if(length(q) > 1){
      here.q = q[item.j]
    } else {
      here.q = q
    }

    if(length(seed.size) > 1){
      here.seed.size = seed.size[item.j]
    } else {
      here.seed.size = seed.size
    }

    PMF.array = d_outbreakdistn(
      x = here.seed.size:(floor(here.q)), seed.size = here.seed.size,
      epi.para = epi.para,
      offspring.type = offspring.type, is.log = FALSE
    )
    CDF.value = sum(PMF.array)
    CDF.array = c(CDF.array, CDF.value)
  }

  if(lower.tail){
    if(is.log){
      CDF.array = log(CDF.array)
    } else {
      CDF.array = CDF.array
    }
  } else {
    if(is.log){
      CDF.array = log(1 - CDF.array)
    } else {
      CDF.array = 1 - CDF.array
    }
  }
  return(CDF.array)
}

#'
#' @rdname d_outbreakdistn
q_outbreakdistn = function(
  p = 0.80, seed.size = 1,
  epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
  offspring.type = 'D', lower.tail = TRUE, upr.limit = 1000#, is.log = FALSE
){
  if(seed.size < 1){
    stop('the value of the function argument  \"seed.size\"  must NOT be less than  1.')
  }

  can.outbreaksize.array = seed.size:floor(upr.limit)
  PMF.array = d_outbreakdistn(
    x = can.outbreaksize.array, seed.size = seed.size,
    epi.para = epi.para,
    offspring.type = offspring.type, is.log = FALSE
  )
  CDF.array = cumsum(PMF.array)

  outbreaksize.array = NULL
  if(lower.tail){
    for (p.j in 1:length(p)) {#      p.j = 1
      outbreaksize.value = can.outbreaksize.array[which(CDF.array >= p[p.j])[1]]
      if(is.na(outbreaksize.value)){
        outbreaksize.value = paste0('> ', upr.limit)
      }
      outbreaksize.array = c(outbreaksize.array, outbreaksize.value)
    }
  } else {
    CDF.array = 1 - max(CDF.array) + cumsum(rev(PMF.array))
    for (p.j in 1:length(p)) {#      p.j = 100
      outbreaksize.value = rev(can.outbreaksize.array)[which(CDF.array >= p[p.j])[1]]
      if(outbreaksize.value == upr.limit){
        outbreaksize.value = paste0('> ', upr.limit)
      }
      outbreaksize.array = c(outbreaksize.array, outbreaksize.value)
    }
  }

  return(as.character(outbreaksize.array))
}

#'
#' @rdname d_outbreakdistn
r_outbreakdistn = function(
  n = 10, seed.size = 1,
  epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
  offspring.type = 'D', upr.limit = 1000#, is.log = FALSE
){
  if(seed.size < 1){
    stop('the value of the function argument  \"seed.size\"  must NOT be less than  1.')
  }

  can.outbreaksize.array = seed.size:floor(upr.limit)
  PMF.array = d_outbreakdistn(
    x = can.outbreaksize.array, seed.size = seed.size,
    epi.para = epi.para,
    offspring.type = offspring.type, is.log = FALSE
  )
  CDF.array = cumsum(PMF.array)

  p = runif(n = n)
  outbreaksize.array = NULL
  for (p.j in 1:length(p)) {#      p.j = 1
    outbreaksize.value = can.outbreaksize.array[which(CDF.array >= p[p.j])[1]]
    if(is.na(outbreaksize.value)){
      outbreaksize.value = paste0('> ', upr.limit)
    }
    outbreaksize.array = c(outbreaksize.array, outbreaksize.value)
  }

  return(as.character(outbreaksize.array))
}






















#' The "20/80" rule
#'
#' @description
#' To calculate proportion of (\code{Q}) offspring cases generated from proportion of (\code{P}) the most infectious index cases with pre-defined epidemiological parameters for the offspring distribution.
#'
#' @param P,Q A scalar, or a vector of probability (i.e., ranging from 0 to 1).
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
#' @param n.seed A positive integer, for the number of seeds used to solve \code{P} or \code{Q} numerically.
#' By default, \code{n.seed = 1000}, and no need to change the default setting here unless for special reasons.
#'
#' @return
#' Function \code{tailoffspringQ()} returns the proportion of (\code{Q}) offspring cases generated from proportion of (\code{P}) index cases, where (\code{P}) is given.
#'
#' Function \code{mostinfectiousP()} returns the proportion of (\code{P}) index cases that generated proportion of (\code{Q}) offspring cases, where (\code{Q}) is given.
#'
#' @import Delaporte
#'
#' @export
#' tailoffspringQ mostinfectiousP
#'
#' @note
#' When \code{n.seed} is large, e.g., \code{n.seed} > 100000, the functions could take minutes to complete.
#' As such, we do not recommend the users to change the default setting of \code{n.seed} unless for special reasons.
#'
#' Each parameter in \code{epi.para = list(mean = ?, disp = ?, shift = ?)} should be a scalar, which means vector is not allowed here.
#'
#'
#' @seealso
#' \code{\link[modelSSE:d_offspringdistn]{d_offspringdistn}}
#'
#' @references
#' Lloyd-Smith JO, Schreiber SJ, Kopp PE, Getz WM. Superspreading and the effect of individual variation on disease emergence. *Nature*. 2005;438(7066):355-359.
#' \doi{10.1038/nature04153}
#'
#' Endo A, Abbott S, Kucharski AJ, Funk S. Estimating the overdispersion in COVID-19 transmission using outbreak sizes outside China. *Wellcome Open Research*. 2020;5:67.
#' \doi{10.12688/wellcomeopenres.15842.3}
#'
#' Adam DC, Wu P, Wong JY, Lau EH, Tsang TK, Cauchemez S, Leung GM, Cowling BJ. Clustering and superspreading potential of SARS-CoV-2 infections in Hong Kong. *Nature Medicine*. 2020;26(11):1714-9.
#' \doi{10.1038/s41591-020-1092-0}
#'
#' Zhao S, Chong MK, Ryu S, Guo Z, He M, Chen B, Musa SS, Wang J, Wu Y, He D, Wang MH. Characterizing superspreading potential of infectious disease: Decomposition of individual transmissibility. *PLoS Computational Biology*. 2022;18(6):e1010281.
#' \doi{10.1371/journal.pcbi.1010281}
#'
#'
#' @examples
#'
#' \donttest{
#'
#' ## reproducing the results in Endo, et al. (2020) https://doi.org/10.12688/wellcomeopenres.15842.3,
#' ## where ~80% offspring cases were generated from ~10% index cases
#' ## with parameters R of ~2.5 (ranging from 2 to 3) and
#' ## k of ~0.1 (ranging from 0.05 to 0.20) under NB distribution.
#' tailoffspringQ(
#'   P = 0.10,
#'   epi.para = list(mean = 2.5, disp = 0.10, shift = 0.2),
#'   offspring.type = "NB"
#' )
#' mostinfectiousP(
#'   Q = 0.80,
#'   epi.para = list(mean = 2.5, disp = 0.10, shift = 0.2),
#'   offspring.type = "NB"
#' )
#'
#'
#' ## reproducing the results in Adam, et al. (2020) https://doi.org/10.1038/s41591-020-1092-0,
#' ## where ~80% offspring cases were generated from ~19% index cases
#' ## with parameters R of 0.58 and k of 0.43 under NB distribution.
#' tailoffspringQ(
#'   P = 0.19,
#'   epi.para = list(mean = 0.58, disp = 0.43, shift = 0.2),
#'   offspring.type = "NB"
#' )
#' mostinfectiousP(
#'   Q = 0.80,
#'   epi.para = list(mean = 0.58, disp = 0.43, shift = 0.2),
#'   offspring.type = "NB"
#' )
#'
#' }
#'
tailoffspringQ = function(
  P = 0.20,
  epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
  offspring.type = 'D', n.seed = 1000
){
  #int.array = 0:1000
  int.array = 0:n.seed
  part.seed = ceiling(n.seed / 30)

  #prop.sourcecase.array = unique(sort(c(P, seq(0.01,0.10,length.out = 21), seq(0.10,0.90,length.out = 61), seq(0.90,0.99,length.out = 21))))
  prop.sourcecase.array = unique(sort(c(
    P,
    seq(0.001,0.010,length.out = part.seed *1 +1),
    seq(0.010,0.100,length.out = part.seed *2 +1),
    seq(0.100,0.900,length.out = part.seed *25 +1),
    seq(0.900,0.990,length.out = part.seed *2 +1),
    seq(0.990,0.999,length.out = part.seed *1 +1)
  )))

  delap.para = convert.epipara.to.delappara(para = epi.para, offspring.type = offspring.type)
  pm.array = Delaporte::ddelap(
    x = int.array, log = FALSE,
    alpha = delap.para$alpha, beta = delap.para$beta, lambda = delap.para$lambda
  )
  expectation.array = pm.array * int.array#; cumsum(expectation.array)

  target.index.array = apply(X = as.matrix(prop.sourcecase.array), MARGIN = 1, FUN = function(z){
    which.min(z > cumsum(pm.array))[1]
  })
  target.weight.array = 1 - (cumsum(pm.array)[target.index.array] - prop.sourcecase.array) / pm.array[target.index.array]

  prop.secondcase.array = 1 - rev(apply(X = as.matrix(1:length(target.index.array)), MARGIN = 1, FUN = function(z){
    y = target.index.array[z]
    sum(expectation.array[1:y]) - expectation.array[y]*(1-target.weight.array[z])
  })) / epi.para$mean
  #    plot(c(0,prop.sourcecase.array,1), c(0,prop.secondcase.array,1))#, type = 'l'

  return(prop.secondcase.array[match(P, prop.sourcecase.array)])
}

#'
#' @rdname tailoffspringQ
mostinfectiousP = function(
  Q = 0.80,
  epi.para = list(mean = 1, disp = 0.5, shift = 0.2),
  offspring.type = 'D', n.seed = 1000
){
  #int.array = 0:1000
  int.array = 0:n.seed

  delap.para = convert.epipara.to.delappara(para = epi.para, offspring.type = offspring.type)
  pm.array = Delaporte::ddelap(
    x = int.array, log = FALSE,
    alpha = delap.para$alpha, beta = delap.para$beta, lambda = delap.para$lambda
  )
  expectation.array = pm.array * int.array#; cumsum(expectation.array)

  sol.index.prop.array = NULL
  for (Q.j in 1:length(Q)) {#       Q.j = 1
    minor.infectee.num = epi.para$mean * (1 - Q[Q.j])
    sol.int.index = which(cumsum(expectation.array) >= minor.infectee.num)[1]

    excess.infectee.num = cumsum(expectation.array)[sol.int.index] - minor.infectee.num
    excess.infector.prop = (excess.infectee.num / expectation.array[sol.int.index]) * pm.array[sol.int.index]
    threshold.int = int.array[sol.int.index]
    sol.index.prop = Delaporte::pdelap(
      q = threshold.int -0, lower.tail = FALSE,
      alpha = delap.para$alpha, beta = delap.para$beta, lambda = delap.para$lambda
    ) + excess.infector.prop

    sol.index.prop.array = c(sol.index.prop.array, sol.index.prop)
  }

  return(c(sol.index.prop.array))
}





















