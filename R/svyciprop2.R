#' Survey-adjusted confidence intervals for proportions
#'
#' @param design Survey design object
#' @param var Variable for which the proportion and confidence interval should be calculated; may take values 0, 1, or NA
#' @param method Confidence interval calculation method. May be Logit, Agresti-Coull, Clopper-Pearson, Fleiss, Jeffreys, Wilson, or Wald.
#' @param level Confidence level, 100 - alpha
#' @param subset_condition Optional condition for subsetting the data, e.g. "var == 1"
#' @param adjust Adjust effective sample size for confidence interval calculations, per Dean and Pagano (2015) and Korn and Graubard (1998)
#' @param truncate If TRUE then the design effect (DEFF) is not allowed to be lower than 1
#'
#' @return A dataset with a survey-adjusted proportion estimate, confidence interval, and other parameters
#'
#' @import survey
#' @import dplyr
#' @import withr
#' @import rlang

svyciprop2 <- function(
    design,
    var,
    method,
    level = 95,
    subset_condition = NULL,
    adjust = TRUE,
    truncate = TRUE
){

  # Set survey design handling of lone PSUs (within function only)
  withr::local_options(list(survey.lonely.psu = "adjust"))

  dat <- dat_full <- design

  # Halt if not a survey design object
  if (!(inherits(dat, "survey.design") | inherits(dat, "survey.design2"))){
    stop("The input for design is not a survey design object created by survey::svydesign.")
  }

  # Halt if variable of interest not in svydf
  if (!any(names(dat$variables) == var)){
    stop(paste0("Variable ", var, " not found in ", design, "."))
  }

  # Halt if the variable of interest takes values other than 0, 1, or missing
  # outcome1 <- get(var, dat$variables)
  dat$variables$outcome <- get(var, dat$variables)

  if (any(! dat$variables$outcome %in% c(0, 1, NA))){

    stop(paste0("To use this function, the variable ", var,
                " should contain only 0s, 1s, or missing values."))
  }

  # Halt if CI method isn't valid
  if (!method %in% c(
    "Wald", "Wilson", "Clopper", "Clopper-Pearson", "Exact",
    "Jeffreys", "Agresti", "Agresti-Coull", "Logit", "Fleiss", "Wilsoncc"
  )){
   stop("Method options are: Wald, Wilson, Clopper, Clopper-Pearson, Exact, Jeffreys, Agresti, Agresti-Coull, Logit, Fleiss, or Wilsoncc")
  }

  if (!is.null(subset_condition)){
    dat <- subset(dat, eval(rlang::parse_expr(subset_condition)))
  }

  outcome2 <- get(var, dat$variables)
  outcome2 <- outcome2[!is.na(outcome2)]

  if (length(outcome2[outcome2 == 1]) == length(outcome2) |
      length(outcome2[outcome2 == 0]) == length(outcome2)){

    if (length(outcome2[outcome2 == 1]) == length(outcome2)){
      phat <- 1
      se <- 0
    }

    if (length(outcome2[outcome2 == 0]) == length(outcome2)){
      phat <- 0
      se <- 0
    }

  } else {
    # Given survey design, calculate proportion and standard error
    params <- survey::svyciprop( # gives same est/se as svymean
      formula = ~outcome, design = dat, df = degf(dat), na.rm = TRUE) |>
      suppressWarnings()

    phat <- as.numeric(coef(params))
    se <- as.numeric(SE(params))
  }

  # Weighted N
  nwtd <- sum(svytable(~outcome, dat))

  # Weighted number with outcome = 1
  nwtd_est <- data.frame(svytable(~outcome, dat))
  nwtd_est <- nwtd_est[nwtd_est$outcome == 1,]$Freq

  # If empty cell for outcome==1, set nwtd_est to 0
  if (length(nwtd_est) == 0){nwtd_est <- 0}

  # Calculate degrees of freedom
  df_n <- nrow(filter(dat$variables, outcome %in% c(0, 1)))

  # Cluster variable
  clustvar <- names(dat_full$cluster)

  # From design: are there strata and clusters?
  has_strata <- dat$has.strata
  has_clusters <- ifelse(nrow(unique(dat$cluster)) == nrow(dat), FALSE, TRUE)

  # Stratum variable
  if (dat_full$has.strata == TRUE){
    stratvar <- names(dat_full$strata)
  } else {
    stratvar <- NULL
  }

  # Count the number of strata and clusters involved in the prevalence calculation
  if (dat$has.strata == TRUE){
    df_strata <- dat$strata |> unique() |> nrow()
  } else {
    df_strata <- 0
  }

  # The df calculation should include ALL clusters in strata where there are
  # respondents who match the `if' criteria; Use tempvars to be sure to count
  # all those clusters, even though some of the clusters may NOT contain
  # respondents who match the `if' criteria.
  #
  # This point is mentioned specifically in West, B. T., Berglund, P., &
  # Heeringa, S. G. (2008). A closer examination of subpopulation analysis of
  # complex-sample survey data. Stata J, 8(4), 520-31.

  if (has_clusters == TRUE & has_strata == TRUE){

    strata_in_calculation <- unique(dat$strata)[,1]

    dat_temp <- dat_full
    dat_temp$variables$tempstrat <- get(stratvar, dat_full$variables)
    dat_temp$variables$tempclust <- get(clustvar, dat_full$variables)

    dat_temp <- dat_temp$variables |>
      filter(tempstrat %in% strata_in_calculation)

    df_cluster <- length(unique(dat_temp$tempclust))

  } else if (has_clusters == TRUE & has_strata == FALSE){
    df_cluster <- dat$cluster |> unique() |> nrow()
  } else if (has_clusters == FALSE){
    df_cluster <- 0
  }

  # Degrees of freedom ----

  # If strata and clusters are specified, then df = # clusters - # strata
  if (has_strata == TRUE & has_clusters == TRUE){
    df <- df_cluster - df_strata
  }

  # If no clusters, then df = N - # strata
  if (has_strata == TRUE & has_clusters == FALSE){
    df <- nrow(dat) - df_strata
  }

  # If not stratified, then df = # clusters - 1
  if (has_strata == FALSE & has_clusters == TRUE){
    df <- df_cluster - 1
  }

  # If no clusters or strata, then df = N - 1
  if (has_strata == FALSE & has_clusters == FALSE){
    df <- nrow(dat) - 1
  }

  # If the standard error is zero, there is no clustering effect, so
  # we set df = N - 1
  if (se == 0){
    df <- nrow(dat) - 1
  }

  # CI calculation section -----

  ci_df <- data.frame(level = level)

  ci_df <- ci_df |>
    mutate(
      phat = phat,
      pqhat = phat * (1-phat),
      se = se,
      n = df_n,
      neff = pqhat/(se^2)
    )

  pstring <- sprintf("%.7f", phat)
  sestring <- sprintf("%.7f", se)

  # If phat is 0 or 1 to seven digits or if stderr is 0 to seven digits
  # then we have homogeneity; set neff to n
  if (pstring == "1.0000000" | pstring == "0.0000000" | sestring == "0.0000000"){
    ci_df$neff <- df_n
  }

  ci_df <- ci_df |>
    mutate(
      df_N = round(n, 1),
      df = ifelse(df != -999, df, df_N - 1),
      DEFF = (df_N - 1)/neff
    )

  if (pstring == "1.0000000" | pstring == "0.0000000" | sestring == "0.0000000"){
    ci_df$DEFF <- 1
  }

  ci_df <- ci_df |>
    mutate(
      ao2 = (100 - level)/100/2,
      zao2 = qnorm(1 - ao2),
      acc = (zao2^2)/2, # Agresti-Coull c
      tdfNao2 = ifelse(df_N > 1, qt(1 - ao2, df_N - 1), NA), # Would use for K&G 1998 neff
      tddfao2 = qt(1 - ao2, df)
    )

  # Adjust effective sample size if user has specified the adjust option
  if (adjust == TRUE){
    ci_df <- ci_df |>
      mutate(neff = neff * (zao2/tddfao2)^2)

    # For K&G 1998 neff, instead would use:
    # mutate(neff = neff * (tdfNao2/tddfao2)^2)
  }

  # Replace effective sample size with actual sample size if DEFF < 1 and
  # user asked to truncate DEFF to be >= 1
  if (truncate == TRUE & ci_df$DEFF[1] < 1){
    ci_df <- ci_df |>
      mutate(
        neff = n,
        DEFF = 1
      )
  }

  ## Wald ----

  if (method == "Wald"){
    # If p is 0 or 1 or stderr = 0, skip Wald calculation and go to
    # Clopper-Pearson

    if (pstring == "1.0000000" | pstring == "0.0000000" | sestring == "0.0000000"){
      method <- "Clopper-Pearson"
      warning("Reverting to Clopper-Pearson confidence interval because Wald interval is not defined.")
    } else {
      ci_df <- ci_df |>
        mutate(lcb_2sided = phat - abs(zao2 * sqrt(pqhat/neff)),
               ucb_2sided = phat + abs(zao2 * sqrt(pqhat/neff)))
    }
  } # end Wald

  ## Logit ----
  if (method == "Logit"){

    # If p is 0 or 1 or stderr = 0, skip Logit calculation and go to
    # Clopper-Pearson

    if (pstring == "1.0000000" | pstring == "0.0000000" | sestring == "0.0000000"){
      method <- "Clopper-Pearson"
      warning("Reverting to Clopper-Pearson confidence interval because Logit interval is not defined.")

    } else {
      ci_df <- ci_df |>
        mutate(
          term1 = log(phat/(1-phat)),
          term2 = zao2/sqrt(neff*pqhat),
          combo1 = term1 - term2,
          combo2 = term1 + term2,
          lcb_2sided = exp(combo1)/(1 + exp(combo1)),
          ucb_2sided = exp(combo2)/(1 + exp(combo2))
        ) |>
        select(-c(term1, term2, combo1, combo2))
    }
  } # end Logit

  ## Wilson ----
  if (method == "Wilson"){
    ci_df <- ci_df |>
      mutate(
        term1 = phat + ((zao2)^2)/(2*neff),
        term2 = zao2 * sqrt((pqhat/neff) + ((zao2)^2)/((2*neff)^2)),
        term3 = 1 + ((zao2)^2)/neff,
        lcb_2sided = (term1 - term2)/term3,
        ucb_2sided = (term1 + term2)/term3
      ) |> select(-c(term1, term2, term3))
  } # end Wilson

  ## Jeffreys ----
  if (method == "Jeffreys"){
    ci_df <- ci_df |>
      mutate(
        x = phat*neff,
        alpha1 = x + 0.5,
        beta1 = neff - x + 0.5,
        lcb_2sided = qbeta(ao2, alpha1, beta1),
        ucb_2sided = qbeta(1 - ao2, alpha1, beta1),
        lcb_2sided = ifelse(phat == 0, 0, lcb_2sided),
        ucb_2sided = ifelse(phat == 1, 1, ucb_2sided)
      ) |> select(-c(x, alpha1, beta1))
  } # end Jeffreys

  ## Agresti-Coull ----
  if (method %in% c("Agresti", "Agresti-Coull")){
    ci_df <- ci_df |>
      mutate(
        xtilde  = phat * neff + acc,
        ntilde  = neff + 2 * acc,
        ptilde  = xtilde/ntilde,
        pqtilde = ptilde * (1 - ptilde),

        lcb_2sided = ptilde - zao2 * sqrt(pqtilde/ntilde),
        ucb_2sided = ptilde + zao2 * sqrt(pqtilde/ntilde)
      ) |> select(-c(xtilde, ntilde, ptilde, pqtilde))
  } # end Agresti-Coull

  ## Fleiss ----
  if (method %in% c("Fleiss", "Wilsoncc")){
    ci_df <- ci_df |>
      mutate(
        term1l = 2*neff*phat + zao2^2 - 1,
        term2l = zao2*sqrt(zao2^2 - ( 2 + (1/neff)) + 4*phat*(neff*(1-phat)+1)),
        term3  = 2*(neff+zao2^2),
        lcb_2sided = (term1l - term2l)/term3,
        term1u = term1l + 2,
        term2u = zao2*sqrt(zao2^2 + ( 2 - (1/neff)) + 4*phat*(neff*(1-phat)-1)),
        ucb_2sided = (term1u + term2u)/term3,
        lcb_2sided = ifelse(phat == 0, 0, lcb_2sided),
        ucb_2sided = ifelse(phat == 1, 1, ucb_2sided)
      ) |> select(-c(term1l, term2l, term1u, term2u, term3))

  } # end Fleiss

  ## Clopper-Pearson ----
  if (method %in% c("Clopper-Pearson", "Clopper", "Exact")){

    # If the sample proportion is 0 or 1,
    # or if the standard error is zero (meaning all clusters have the same
    # observed proportion) then consider the effective sample size to be equal
    # to the actual sample size
    if (round(phat, 10) == 0 | round(phat, 10) == 1 | round(se, 10) == 0){
      ci_df <- ci_df |>
        mutate(neff = n,
               DEFF = 1)
    }

    ci_df <- ci_df |>
      mutate(
        x = phat*neff,
        v1 = 2*x,
        v2 = 2*(neff-x+1),
        v3 = 2*(x+1),
        v4 = 2*(neff-x)) |>
      rowwise() |>
      mutate(
        v1 = max(v1, 2e-10),
        v1 = min(v1, 2e+17),
        v2 = max(v2, 2e-10),
        v2 = min(v2, 2e+17),
        v3 = max(v3, 2e-10),
        v3 = min(v3, 2e+17),
        v4 = max(v4, 2e-10),
        v4 = min(v4, 2e+17)) |>
      ungroup() |>
      mutate(
        fao2 = qf(ao2, v1, v2),
        lcb_2sided = (v1*fao2)/(v2 + v1*fao2),
        f1mao2  = qf(1-ao2, v3, v4),
        ucb_2sided = (v3*f1mao2)/(v4 + v3*f1mao2),
        # If v4 is very small, the UCB ratio is set to missing instead of 1
        # so check for this condition here
        ucb_2sided = ifelse(v4 <= 2e-10, 1, ucb_2sided)
      ) |>
      select(-c(x, v1, v2, v3, v4, fao2, f1mao2))
  } # end Clopper-Pearson

  # Replace infinitesimal values of lcb with 0
  ci_df <- ci_df |>
    mutate(
      lcb_2sided = ifelse(abs(lcb_2sided) <= 1e-10, 0, lcb_2sided)
    )

  out <- ci_df |>
    mutate(method = method) |>
    select(phat, method, level, lcb_2sided, ucb_2sided,
           se, n, neff, df, DEFF)

  dplyr::tibble(out)
}

