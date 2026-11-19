suppressPackageStartupMessages(library('nlme'))
suppressPackageStartupMessages(library('MASS'))
suppressPackageStartupMessages(library('future.apply'))

# I nest things in if (TRUE) just to allow them to be collapsed in the IDEs

num_cores <- max(1, parallelly::availableCores() - 2)
plan(multisession, workers = num_cores)

##############################
#### COVARIANCE FUNCTIONS ####
##############################
if (TRUE) {
  gen_ind_cov <- function(pars) {
      n <- pars$n_trt * pars$n_blk * pars$n_rep
      diag(rep(pars$s2e, n))
  }
  gen_comp_symm_cov <- function(pars) {
      # Initial n_rep x n_rep block matrix with just s2b
      unit_cov <- matrix(pars$s2b, nrow = pars$n_rep, ncol = pars$n_rep)
      # Add s2e to the diagonals of the block matrix
      unit_cov <- unit_cov + diag(rep(pars$s2e, pars$n_rep))
      kronecker(diag(pars$n_trt * pars$n_blk), unit_cov)
  }
  gen_ar1_cov <- function(pars) {
      times <- 1:pars$n_rep
      unit_cov <- pars$s2e * (pars$rho^abs(outer(times, times, "-")))
      kronecker(diag(pars$n_trt * pars$n_blk), unit_cov)
  }
  gen_unstr_cov <- function(pars) {
      kronecker(diag(pars$n_trt * pars$n_blk), pars$unstr_cov)
  }
  gen_all_cov <- function(pars) {
      list(
          Independent = gen_ind_cov(pars), 
          CompoundSymmetric = gen_comp_symm_cov(pars), 
          AR1 = gen_ar1_cov(pars), 
          Unstructured = gen_unstr_cov(pars))
  }
}

##############################
#### SIMULATION FUNCTIONS ####
##############################
if (TRUE) {
  # Generate the data for a specific covariance/parameter set
  gen_data <- function(pars, trt_means, Sigma) {
      block <- unlist(lapply(1:(pars$n_trt*pars$n_blk), 
          \(x) rep(x, times = pars$n_rep)))
      reps <- rep(1:5, pars$n_trt*pars$n_blk)
      means <- unlist(lapply(trt_means, 
          \(x) rep(x, times = pars$n_blk*pars$n_rep)))
      trt_labels <- paste0("Trt", 1:pars$n_trt)
      treatment_factor <- unlist(lapply(trt_labels, 
          \(x) rep(x, times = pars$n_blk*pars$n_rep)))
      errors <- mvrnorm(n = 1, 
          mu = rep(0, pars$n_trt*pars$n_blk*pars$n_rep), Sigma = Sigma)
      dt <- data.frame(
          Treatment = means, 
          TrtGroup = as.factor(treatment_factor), 
          Block = as.factor(block),
          Repetition = as.factor(reps),
          Response = means + errors)
  }

  # Given the data, fit all models and pick the best
  pick_best_model_fit <- function(dt) {
    # Independent
    ind_res <- tryCatch({
            mod_ind <- gls(Response ~ TrtGroup, data = dt)
            bic_ind <- BIC(mod_ind)
            list(mod = mod_ind, bic = bic_ind)
        }, error = function(e) {
            list(mod = NA, bic = NA)
        })

    # Compound Symmetric
    cs_res <- tryCatch({
            mod_cs <- gls(Response ~ TrtGroup, correlation = corCompSymm(form = ~ 1 | Block), data = dt)
            bic_cs <- BIC(mod_cs)
            list(mod = mod_cs, bic = bic_cs)
        }, error = function(e) {
            list(mod = NA, bic = NA)
        })
    
    # AR1
    ar1_res <- tryCatch({
            mod_ar1 <- gls(Response ~ TrtGroup, correlation = corAR1(form = ~ as.numeric(Repetition) | Block), data = dt)
            bic_ar1 <- BIC(mod_ar1)
            list(mod = mod_ar1, bic = bic_ar1)
    }, error = function(e) {
            list(mod = NA, bic = NA)
        })
    
    # Unstructured
    unstr_res <- tryCatch({
            mod_unstr <- gls(Response ~ TrtGroup, correlation = corSymm(form = ~ 1 | Block), weights = varIdent(form = ~ 1 | Repetition), data = dt)
            bic_unstr <- BIC(mod_unstr)
            list(mod = mod_unstr, bic = bic_unstr)
    }, error = function(e) {
            list(mod = NA, bic = NA)
        })
    
    # Create all the necessary lists separately
    mod_list <- list(
        Ind = ind_res$mod, 
        CS = cs_res$mod, 
        AR1 = ar1_res$mod, 
        Unstr = unstr_res$mod
    )
    bic_list <- c(
        Ind = ind_res$bic, 
        CS = cs_res$bic, 
        AR1 = ar1_res$bic, 
        Unstr = unstr_res$bic
    )
    hyp_test_list <- c(
        IndRej = 0,
        CSRej = 0,
        AR1Rej = 0,
        UnstrRej = 0
    )
    fail_list <- c(
        IndFail = is.na(ind_res$bic),
        CSFail = is.na(cs_res$bic),
        AR1Fail = is.na(ar1_res$bic),
        UnstrFail = is.na(unstr_res$bic)
    )

    # Identify best model
    bic_list[1:4] <- bic_list[1:4] == min(bic_list[1:4], na.rm = TRUE)

    # With the chosen best model, check the pvalue and see if we reject the Null Hypothesis
    mod_summary <- summary(mod_list[[which(bic_list == TRUE)]])
    hyp_test_list[which(bic_list == TRUE)] <- (mod_summary$tTable["TrtGroupTrt2", "p-value"] <= .05)

    # Return the joined list of all the indicators
    return(c(bic_list, hyp_test_list, fail_list))
  }

  # Do one interation of generating data and picking the best model
  run_one_sim <- function(pars, trt_means, Sigma) {
      sim_data <- gen_data(pars, trt_means, Sigma)
      pick_best_model_fit(sim_data)
  }

  # Do one interation of generating data and picking the best model
  run_one_sim_drop_20 <- function(pars, trt_means, Sigma) {
        data <- gen_data(pars, trt_means, Sigma)
        data <- data[sample.int(nrow(data), nrow(data)*.8), ] 
        pick_best_model_fit(data)
  }

  run_sim_set <- function(pars, trt_means, n_sims, drop_20, covs) {
    sim_function <- run_one_sim
    if (drop_20) {
      sim_function <- run_one_sim_drop_20
    }
    i <- 1
    dt <- lapply(covs, \(x) {
      sims <- future_sapply(1:n_sims, \(sim) sim_function(pars, trt_means, x), future.seed = TRUE)
      print(paste0("Finished Covariance: ", i))
      i <<- i + 1
      colSums(t(sims), na.rm = TRUE)
    })
    print(paste0("Finished Simulation | trt2=", trt_means[2], " | Drop20=", drop_20))
    sim_data <- as.data.frame(do.call(rbind, dt))
    sim_data$GeneratingModel <- names(covs)
    sim_data$Trt2 <- trt_means[2]
    sim_data$Drop20 <- drop_20
    sim_data <- sim_data[, c(
        "GeneratingModel", "Trt2", "Drop20", 
        "Ind", "CS", "AR1", "Unstr", 
        "IndRej", "CSRej", "AR1Rej", "UnstrRej", 
        "IndFail", "CSFail", "AR1Fail", "UnstrFail")]
    sim_data
  }
}

###############################
#### SIMULATION PARAMETERS ####
###############################
if (TRUE) {
  # Initialize simulation parameters
  sim_params <- list(
      n_trt = 2,  # Treatments
      n_blk = 10, # Blocks
      n_rep = 5,  # Repetitions per Treatment/Block combo
      s2e = 10,   # Sigma^2_e
      s2b = 12,   # Sigma^2_b
      rho = 0.8,  # rho for AR1 covariance
      unstr_cov = matrix(c(
          15, 10, 12, 12, 10, 
          10, 18, 5, 13, 8, 
          12, 5, 20, 5, 3, 
          12, 13, 5, 15, 10, 
          10, 8, 3, 10, 12), nrow = 5, byrow = TRUE)
  )

  # List of the covariance functions with their names
  cov_list <- list(
      Independent = gen_ind_cov, 
      CompoundSymmetric = gen_comp_symm_cov, 
      AR1 = gen_ar1_cov, 
      Unstructured = gen_unstr_cov
  )
}

############################
#### RUN THE SIMULATION ####
############################
set.seed(537)
n_sims <- 500
covs <- gen_all_cov(sim_params)

trt_means <- c(50, 50)
trt2_50_incl_100 <- run_sim_set(sim_params, trt_means, n_sims, drop_20 = FALSE, covs)
trt2_50_incl_80 <- run_sim_set(sim_params, trt_means, n_sims, drop_20 = TRUE, covs)
trt_means <- c(50, 55)
trt2_55_incl_100 <- run_sim_set(sim_params, trt_means, n_sims, drop_20 = FALSE, covs)
trt2_55_incl_80 <- run_sim_set(sim_params, trt_means, n_sims, drop_20 = TRUE, covs)

res <- rbind(trt2_50_incl_100, trt2_50_incl_80, trt2_55_incl_100, trt2_55_incl_80)
rownames(res) <- NULL
write.csv(res, 'sim_data.csv', row.names = FALSE)


