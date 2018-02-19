#' @name    gvc_herit
#' @aliases gvc_herit
#' @title   Heritability
#' @description gvc_herit computes model based genetic heritability
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
#' herit <-
#'   gvc_herit(
#'            y    = Response
#'          , rep  = Rep
#'          , geno = Variety
#'          , env  = Env
#'          , data = df1
#'          )
#' herit
#'
#' library(eda4treeR)
#' data(DataExam6.2)
#' herit <-
#'   gvc_herit(
#'            y    = Dbh.mean
#'          , rep  = Replication
#'          , geno = Family
#'          , env  = Province
#'          , data = DataExam6.2
#'          )
#' herit

gvc_herit <- function(y, x = NULL, rep, geno, env, data) {

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


  fm1 <- lme4::lmer(
            formula = Mean ~ rep + env + (1|geno)
          , REML    = TRUE
          , data    = df1
        )

  VarCor           <- as.data.frame(lme4::VarCorr(fm1))
  rownames(VarCor) <- VarCor$grp
  sigma2f          <- c(VarCor["geno", "vcov"])

  HM   <- function(x){ length(x)/sum(1/x) }
  w    <- HM(df1$Count)
  b    <- anova(
               lm(
                  formula = Mean ~ rep + geno
                , data    = df1
                )
               )
  s2       <- b[["Mean Sq"]][length(b[["Mean Sq"]])]
  sigma2t  <- mean(df1$Var)
  sigma2m  <- s2-(sigma2t/w)
  h2       <- (sigma2f/(0.3))/(sigma2t + sigma2m + sigma2f)
  out      <- list("h2" = h2)
  return(out)
}
