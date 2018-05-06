#' @name    gvc_pvar
#' @aliases gvc_pvar
#' @title   Phenotypic Variance
#' @description gvc_pvar computes phenotypic variances
#' for given traits of different gentypes from replicated data using methodology
#'  explained by Burton, G. W. & Devane, E. H. (1953) and Allard, R.W. (2010).
#'
#' @param .data  data.frame
#' @param .y      Response
#' @param .x      Covariate by default NULL
#' @param .rep    Repliction
#' @param .gen   gentypic Factor
#' @param .env    Environmental Factor
#'
#'
#' @return Phenotypic Variance
#'
#'
#' @author
#' \enumerate{
#'          \item  Sami Ullah (\email{samiullahuos@@gmail.com})
#'          \item  Muhammad Yaseen (\email{myaseen208@@gmail.com})
#'          }
#'
#' @references
#' \enumerate{
#'          \item R.K. Singh and B.D.Chaudhary
#'          \emph{Biometrical Methods in Quantitative Genetic Analysis}.
#'                Kalyani Publishers, New Delhi
#'              }
#' \enumerate{
#'          \item Williams, E.R., Matheson, A.C. and Harwood, C.E. (2002).\emph{Experimental Design and Analysis for Tree Improvement}.
#'                CSIRO Publishing.
#'              }
#'
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @import lme4
#' @import eda4treeR
#' @importFrom stats anova lm var
#'
#' @export
#'
#' @examples
#' set.seed(12345)
#' Response <- c(
#'                rnorm(48, mean = 15000, sd = 500)
#'              , rnorm(48, mean =  5000, sd = 500)
#'              , rnorm(48, mean =  1000, sd = 500)
#'              )
#' Rep      <- as.factor(rep(1:3, each = 48))
#' Variety  <- gl(n = 4, k =  4, length = 144, labels = letters[1:4])
#' Env      <- gl(n = 3, k = 16, length = 144, labels = letters[1:3])
#' df1      <- data.frame(Response, Rep, Variety, Env)
#'
#' # Penotypic Variance
#' pvar1 <-
#'   gvc_pvar(
#'           .data  = df1
#'          , .y    = Response
#'          , .rep  = Rep
#'          , .gen  = Variety
#'          , .env  = Env
#'           )
#' pvar1
#'
#' library(eda4treeR)
#' data(DataExam6.2)
#' pvar2 <-
#'   gvc_pvar(
#'           .data  = DataExam6.2
#'          , .y    = Dbh.mean
#'          , .rep  = Replication
#'          , .gen  = Family
#'          , .env  = Province
#'          )
#' pvar2

gvc_pvar <- function(.data, .y, .x = NULL, .rep, .gen, .env) {
  UseMethod("gvc_pvar")
}


#' @export
#' @rdname gvc_pvar

gvc_pvar.default <- function(.data, .y, .x = NULL, .rep, .gen, .env) {

  .y    <- deparse(substitute(.y))
  .rep  <- deparse(substitute(.rep))
  .gen  <- deparse(substitute(.gen))
  .env  <- deparse(substitute(.env))


  df1 <- tibble::as_tibble(data.frame(
    Env = factor(.data[[.env]])
    , Gen = factor(.data[[.gen]])
    , Rep = factor(.data[[.rep]])
    , Y   = .data[[.y]]
  ))


  df2 <-
    df1 %>%
    dplyr::group_by(Rep, Gen, Env)%>%
    dplyr::summarize(
      Mean  = mean(Y)
      , Var   = var(Y)
      , Count = length(Y)
    )


  fm1 <-
    lme4::lmer(
      formula = Mean ~ Rep + Env + (1|Gen)
      , REML    = TRUE
      , data    = df2
    )

  VarCor           <- as.data.frame(lme4::VarCorr(fm1))
  rownames(VarCor) <- VarCor$grp
  sigma2f          <- c(VarCor["Gen", "vcov"])
  gvar             <- sigma2f

  HM   <- function(x){ length(x)/sum(1/x) }
  w    <- HM(df2$Count)
  b    <- anova(
    lm(
        formula = Mean ~ Rep + Gen
      , data    = df2
    )
  )
  s2       <- b[["Mean Sq"]][length(b[["Mean Sq"]])]
  sigma2t  <- mean(df2$Var)
  sigma2m  <- s2-(sigma2t/w)

  pvar  <- (sigma2t + sigma2m + sigma2f)

  out <- list("pvar" = pvar)
  return(out)
}
