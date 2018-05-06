#' @name    gvc_herit
#' @aliases gvc_herit
#' @title   Heritability
#' @description gvc_herit computes model based genetic heritability
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
#' @return Heritability
#'
#'
#' @author
#' \enumerate{
#'          \item  Sami Ullah (\email{samiullahuos@@gmail.com})
#'          \item  Muhammad Yaseen (\email{myaseen208@@gmail.com})
#'          }
#'
#' @references
#'
#'\enumerate{
#'          \item Williams, E.R., Matheson, A.C. and Harwood, C.E. (2002).\emph{Experimental Design and Analysis for Tree Improvement}.
#'                CSIRO Publishing.
#'              }
#'
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @import lme4
#' @importFrom stats anova lm var
#' @import eda4treeR
#'
#' @export
#'
#' @examples
#'
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
#' # Heritability
#' herit1 <-
#'   gvc_herit(
#'           .data  = df1
#'          , .y    = Response
#'          , .rep  = Rep
#'          , .gen  = Variety
#'          , .env  = Env
#'          )
#' herit1
#'
#' library(eda4treeR)
#' data(DataExam6.2)
#' herit2 <-
#'   gvc_herit(
#'           .data  = DataExam6.2
#'          , .y    = Dbh.mean
#'          , .rep  = Replication
#'          , .gen  = Family
#'          , .env  = Province
#'          )
#' herit2

gvc_herit <- function(.data, .y, .x = NULL, .rep, .gen, .env) {
  UseMethod("gvc_herit")
}


#' @export
#' @rdname gvc_herit

gvc_herit.default <- function(.data, .y, .x = NULL, .rep, .gen, .env) {

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
  w    <- HM(df1$Count)
  b    <- anova(
               lm(
                  formula = Mean ~ Rep + Gen
                , data    = df2
                )
               )

  s2       <- b[["Mean Sq"]][length(b[["Mean Sq"]])]
  sigma2t  <- mean(df2$Var)
  sigma2m  <- s2-(sigma2t/w)
  pvar     <- (sigma2t + sigma2m + sigma2f)
  h2       <- (gvar/(0.3))/(pvar)
  out      <- list("h2" = h2)
  return(out)
}
