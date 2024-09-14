#' @name gvc
#' @aliases gvc
#' @title Genotypic Variance
#' @description The `gvc` class calculates genotypic variance, phenotypic variance, and broad-sense heritability from replicated data.
#'
#' @details
#' The `gvc` class uses methods explained by Burton, G. W. & Devane, E. H. (1953) and Allard, R.W. (2010). It includes methods for calculating genetic variance, phenotypic variance, and heritability.
#'
#' @field data A `tibble` containing the data for analysis.
#' @field y The name of the response variable.
#' @field x The name of the covariate (optional).
#' @field rep The name of the replicate factor.
#' @field gen The name of the genotype factor.
#' @field env The name of the environmental factor.
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @import lme4
#' @import rlang
#' @import eda4treeR
#' @importFrom stats anova lm var
#' @import supernova
#' @import tibble
#' @import R6
#'
#' @export
#' @examples
#' df1 <- data.frame(
#'            Response = c(
#'                           rnorm(48, mean = 15000, sd = 500)
#'                         , rnorm(48, mean =  5000, sd = 500)
#'                         , rnorm(48, mean =  1000, sd = 500)
#'                        )
#'         , Rep      = as.factor(rep(1:3, each = 48))
#'         , Variety  = gl(n = 4, k =  4, length = 144, labels = letters[1:4])
#'         , Env      = gl(n = 3, k = 16, length = 144, labels = letters[1:3])
#'         )
#'
#' # Create an instance of the class
#' gvc1 <- gvc$new(
#'            .data = df1
#'          , .y    = Response
#'          , .rep  = Rep
#'          , .gen  = Variety
#'          , .env  = Env
#'          )
#'
#' # Calculate genetic variance (gvar)
#' gvc1$calculate_gvar()
#'
#' # Calculate phenotypic variance (pvar)
#' gvc1$calculate_pvar()
#'
#' # Calculate heritability (h2)
#' gvc1$calculate_herit()
#'
gvc <- R6Class(
  "gvc",

  public = list(
    data = NULL,  # Original data
    y = NULL,     # Dependent variable (response)
    x = NULL,     # Independent variable (optional)
    rep = NULL,   # Replicate factor
    gen = NULL,   # Genotype factor
    env = NULL,   # Environment factor

    #' @description
    #' Initialize the `gvc` class with the data and variable names.
    #' @param .data A `data.frame` containing the data for analysis.
    #' @param .y The response variable.
    #' @param .x The covariate (optional).
    #' @param .rep The replicate factor.
    #' @param .gen The genotype factor.
    #' @param .env The environmental factor.
    #' @return An instance of the `gvc` class.
    initialize = function(.data, .y, .x = NULL, .rep, .gen, .env) {
      self$data <- as_tibble(.data)
      self$y    <- deparse(substitute(.y))
      self$x    <- deparse(substitute(.x))
      self$rep  <- deparse(substitute(.rep))
      self$gen  <- deparse(substitute(.gen))
      self$env  <- deparse(substitute(.env))

      # Validate data and factors
      private$validate_data()
    },

    #' @description
    #' Calculate genetic variance.
    #' @return A list with the genetic variance (`gvar`).
    calculate_gvar = function() {
      df2 <- private$prepare_data()

      # Linear mixed-effects model
      fm1 <- lmer(
        formula = Mean ~ Rep + Env + (1 | Gen),
        REML = TRUE,
        data = df2
      )

      # Extract variance component for Genotype
      VarCor <- as.data.frame(VarCorr(fm1))
      sigma2f <- VarCor[VarCor$grp == "Gen", "vcov"]

      return(list("gvar" = sigma2f))
    },

    #' @description
    #' Calculate phenotypic variance.
    #' @return A list with the phenotypic variance (`pvar`).
    calculate_pvar = function() {
      df2 <- private$prepare_data()

      # Linear mixed-effects model
      fm1 <- lmer(
        formula = Mean ~ Rep + Env + (1 | Gen),
        REML = TRUE,
        data = df2
      )

      VarCor <- as.data.frame(VarCorr(fm1))
      sigma2f <- VarCor[VarCor$grp == "Gen", "vcov"]

      # Harmonic mean function
      HM <- function(x) length(x) / sum(1 / x)
      w  <- HM(df2$Count)

      # ANOVA for mean squares
      anova_results <- anova(lm(Mean ~ Rep + Gen, data = df2))
      s2 <- anova_results[["Mean Sq"]][length(anova_results[["Mean Sq"]])]

      sigma2t <- mean(df2$Var)
      sigma2m <- s2 - (sigma2t / w)

      pvar <- sigma2t + sigma2m + sigma2f
      return(list("pvar" = pvar))
    },

    #' @description
    #' Calculate broad-sense heritability.
    #' @return A list with the heritability (`h2`).
    calculate_herit = function() {
      df2 <- private$prepare_data()

      # Linear mixed-effects model
      fm1 <- lmer(
        formula = Mean ~ Rep + Env + (1 | Gen),
        REML = TRUE,
        data = df2
      )

      VarCor <- as.data.frame(VarCorr(fm1))
      sigma2f <- VarCor[VarCor$grp == "Gen", "vcov"]

      # Harmonic mean function
      HM <- function(x) length(x) / sum(1 / x)
      w  <- HM(df2$Count)

      # ANOVA for mean squares
      anova_results <- anova(lm(Mean ~ Rep + Gen, data = df2))
      s2 <- anova_results[["Mean Sq"]][length(anova_results[["Mean Sq"]])]

      sigma2t <- mean(df2$Var)
      sigma2m <- s2 - (sigma2t / w)
      pvar <- sigma2t + sigma2m + sigma2f

      h2 <- (sigma2f / pvar)
      return(list("h2" = h2))
    }
  ),

  private = list(

    # Validate data and check for required columns
    validate_data = function() {
      required_columns <- c(self$y, self$rep, self$gen, self$env)
      missing_columns <- setdiff(required_columns, colnames(self$data))
      if (length(missing_columns) > 0) {
        stop(paste("Missing columns in data:", paste(missing_columns, collapse = ", ")))
      }
    },

    # Prepare data for analysis by grouping and summarizing
    prepare_data = function() {
      df1 <- self$data %>%
        mutate(
          Env = factor(.data[[self$env]]),
          Gen = factor(.data[[self$gen]]),
          Rep = factor(.data[[self$rep]])
        ) %>%
        select(Env, Gen, Rep, Y = all_of(self$y))

      df2 <- df1 %>%
        group_by(Rep, Gen, Env) %>%
        summarize(
          Mean  = mean(Y, na.rm = TRUE),
          Var   = var(Y, na.rm = TRUE),
          Count = n(),
          .groups = "drop"
        )

      return(df2)
    }
  )
)
