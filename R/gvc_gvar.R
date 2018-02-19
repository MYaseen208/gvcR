#' @name    gvc_gvar
#' @aliases gvc_gvar
#' @title   Genotypic Variance
#' @description gvc_gvar computes genotypic variances
#' for given traits of different genotypes from replicated data using methodology
#'  explained by Burton, G. W. & Devane, E. H. (1953) (<doi:10.2134/agronj1953.00021962004500100005x>) and Allard R.W. (1960).
#'
#' @param y     Response
#' @param x     Covariate by default NULL
#' @param rep   Repliction
#' @param geno  Genotypic Factor
#' @param env   Environmental Factor
#' @param data  data.frame
#'
#'
#' @return Genotypic Variance
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
#' # Genotypic Variance
#' gvar <-
#'  gvc_gvar(
#'            y    = Response
#'          , rep  = Rep
#'          , geno = Variety
#'          , env  = Env
#'          , data = df1
#'          )
#' gvar
#'
#' library(eda4treeR)
#' data(DataExam6.2)
#' gvar <-
#'   gvc_gvar(
#'            y    = Dbh.mean
#'          , rep  = Replication
#'          , geno = Family
#'          , env  = Province
#'          , data = DataExam6.2
#'          )
#' gvar

gvc_gvar <- function(y, x = NULL, rep, geno, env, data) {

  y         <- enquo(y)
  x         <- enquo(x)
  rep       <- enquo(rep)
  geno      <- enquo(geno)
  env       <- enquo(env)

  df1 <- data %>%
    dplyr::group_by(!! rep, !! geno, !! env)%>%
    dplyr::summarize(
        Mean  = mean(!! y)
      , Var   = var(!! y)
      , Count = length(!! y)
    )

  names(df1) <- c("rep", "geno", "env", "Mean", "Var", "Count")


  fm1 <-
        lme4::lmer(
                     formula = Mean ~ rep + env + (1|geno)
                   , REML    = TRUE
                   , data    = df1
                   )

  VarCor           <- as.data.frame(lme4::VarCorr(fm1))
  rownames(VarCor) <- VarCor$grp
  sigma2f          <- c(VarCor["geno", "vcov"])
  gvar             <- sigma2f
  out  <- list("gvar" = gvar)
  return(out)
}
