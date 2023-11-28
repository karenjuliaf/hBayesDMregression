#' @noRd
#'
write_stan_model <- function(model_name, regression_pars) {
  switch(
    model_name,
    "pvl_decay_regression" = write_stan_model_pvl_decay_regression(regression_pars),
    "orl_regression" = write_stan_model_orl_regression(regression_pars),
    "vpp_regression" = write_stan_model_vpp_regression(regression_pars)
  )
}

write_stan_model_pvl_decay_regression <- function(regression_pars) {
  regression_strings <- list(
    "A" = "A[i] = Phi_approx(mu_pr[1] + covariate_matrix[i] * beta[{{par_index}}]' + sigma[1] * A_pr[i]);",
    "alpha" = "alpha[i] = Phi_approx(mu_pr[2] + covariate_matrix[i] * beta[{{par_index}}]' + sigma[2] * alpha_pr[i]) * 2;",
    "cons" = "cons[i] = Phi_approx(mu_pr[3] + covariate_matrix[i] * beta[{{par_index}}]' + sigma[3] * cons_pr[i]) * 5;",
    "lambda" = "lambda[i] = Phi_approx(mu_pr[4] + covariate_matrix[i] * beta[{{par_index}}]' + sigma[4] * lambda_pr[i]) * 10;"
  )
  non_regression_strings <- list(
    "A" = "A[i] = Phi_approx(mu_pr[1] + sigma[1] * A_pr[i]);",
    "alpha" = "alpha[i] = Phi_approx(mu_pr[2] + sigma[2] * alpha_pr[i]) * 2;",
    "cons" = "cons[i] = Phi_approx(mu_pr[3] + sigma[3] * cons_pr[i]) * 5;",
    "lambda" = "lambda[i] = Phi_approx(mu_pr[4] + sigma[4] * lambda_pr[i]) * 10;"
  )
  non_regression_pars <- setdiff(names(regression_strings), regression_pars)

  main_chunk <- regression_strings[regression_pars]
  main_chunk <- mapply(
    gsub,
    x = main_chunk,
    replacement = as.character(1:length(regression_pars)),
    MoreArgs = list(pattern = "\\{\\{par_index\\}\\}"),
    SIMPLIFY = FALSE
  )
  main_chunk <- c(main_chunk, non_regression_strings[non_regression_pars])
  main_chunk <- main_chunk[sort(names(main_chunk))]

  npars <- length(regression_pars)

  model_code <- glue::glue("
  data {
    int<lower=1> N;
    int<lower=1> T;
    int<lower=1, upper=T> Tsubj[N];
    int choice[N, T];
    real outcome[N, T];

    int<lower=1> ncov;
    matrix[N, ncov] covariate_matrix;
    matrix<lower=0>[{{npars}}, ncov] covariate_precisions;
  }
  transformed data {
    vector[4] initV;
    initV  = rep_vector(0.0, 4);
  }
  parameters {
  // Declare all parameters as vectors for vectorizing
    // Hyper(group)-parameters
    vector[4] mu_pr;
    vector<lower=0>[4] sigma;
    matrix<lower=0>[{{npars}}, ncov] sigma_beta;
    matrix[{{npars}}, ncov] beta;

    // Subject-level raw parameters (for Matt trick)
    vector[N] A_pr;
    vector[N] alpha_pr;
    vector[N] cons_pr;
    vector[N] lambda_pr;
  }
  transformed parameters {
    // Transform subject-level raw parameters
    vector<lower=0, upper=1>[N]  A;
    vector<lower=0, upper=2>[N]  alpha;
    vector<lower=0, upper=5>[N]  cons;
    vector<lower=0, upper=10>[N] lambda;

    for (i in 1:N) {
      {{main_chunk$A}}
      {{main_chunk$alpha}}
      {{main_chunk$cons}}
      {{main_chunk$lambda}}
    }
  }
  model {
    // Hyperparameters
    mu_pr  ~ normal(0, 1);
    sigma  ~ normal(0, 0.2);

    // Regression parameters
    for (i in 1:{{npars}}) {
      sigma_beta[i] ~ exponential(covariate_precisions[i]);
      beta[i] ~ normal(0, sigma_beta[i]);
    }

    // individual parameters
    A_pr      ~ normal(0, 1);
    alpha_pr  ~ normal(0, 1);
    cons_pr   ~ normal(0, 1);
    lambda_pr ~ normal(0, 1);

    for (i in 1:N) {
      // Define values
      vector[4] ev;
      real curUtil;     // utility of curFb
      real theta;       // theta = 3^c - 1

      // Initialize values
      theta = pow(3, cons[i]) -1;
      ev = initV; // initial ev values

      for (t in 1:Tsubj[i]) {
        // softmax choice
        choice[i, t] ~ categorical_logit(theta * ev);

        if (outcome[i, t] >= 0) {  // x(t) >= 0
          curUtil = pow(outcome[i, t], alpha[i]);
        } else {                  // x(t) < 0
          curUtil = -1 * lambda[i] * pow(-1 * outcome[i, t], alpha[i]);
        }

        // decay-RI
        ev *= A[i];
        ev[choice[i, t]] += curUtil;
      }
    }
  }
  generated quantities {
    // For group level parameters
    real<lower=0, upper=1>  mu_A;
    real<lower=0, upper=2>  mu_alpha;
    real<lower=0, upper=5>  mu_cons;
    real<lower=0, upper=10> mu_lambda;

    // For log likelihood calculation
    real log_lik[N];

    // For posterior predictive check
    real y_pred[N, T];

    // Set all posterior predictions to 0 (avoids NULL values)
    for (i in 1:N) {
      for (t in 1:T) {
        y_pred[i, t] = -1;
      }
    }

    mu_A      = Phi_approx(mu_pr[1]);
    mu_alpha  = Phi_approx(mu_pr[2]) * 2;
    mu_cons   = Phi_approx(mu_pr[3]) * 5;
    mu_lambda = Phi_approx(mu_pr[4]) * 10;

    { // local section, this saves time and space
      for (i in 1:N) {
        // Define values
        vector[4] ev;
        real curUtil;     // utility of curFb
        real theta;       // theta = 3^c - 1

        // Initialize values
        log_lik[i] = 0;
        theta = pow(3, cons[i]) -1;
        ev = initV; // initial ev values

        for (t in 1:Tsubj[i]) {
          // softmax choice
          log_lik[i] += categorical_logit_lpmf(choice[i, t] | theta * ev);

          // generate posterior prediction for current trial
          y_pred[i, t] = categorical_rng(softmax(theta * ev));

          if (outcome[i, t] >= 0) {  // x(t) >= 0
            curUtil = pow(outcome[i, t], alpha[i]);
          } else {                  // x(t) < 0
            curUtil = -1 * lambda[i] * pow(-1 * outcome[i, t], alpha[i]);
          }

          // decay-RI
          ev *= A[i];
          ev[choice[i, t]] += curUtil;
        }
      }
    }
  }", .open="{{", .close="}}", .trim=FALSE)

  return(model_code)
}

write_stan_model_orl_regression <- function(regression_pars) {
  regression_strings <- list(
    "Arew" = "Arew[i] = Phi_approx(mu_pr[1] + covariate_matrix[i] * beta[{{par_index}}]' + sigma[1] * Arew_pr[i]);",
    "Apun" = "Apun[i] = Phi_approx(mu_pr[2] + covariate_matrix[i] * beta[{{par_index}}]' + sigma[2] * Apun_pr[i]);",
    "K" = "K[i] = Phi_approx(mu_pr[3] + covariate_matrix[i] * beta[{{par_index}}]' + sigma[3] * K_pr[i]) * 5;",
    "betaF" = "betaF = mu_pr[4] + covariate_matrix[i] * beta[{{par_index}}]' + sigma[4] * betaF_pr[i];",
    "betaP" = "betaP = mu_pr[5] + covariate_matrix[i] * beta[{{par_index}}]' + sigma[5] * betaP_pr[i];"
  )
  non_regression_strings <- list(
    "Arew" = "Arew[i] = Phi_approx(mu_pr[1] + sigma[1] * Arew_pr[i]);",
    "Apun" = "Apun[i] = Phi_approx(mu_pr[2] + sigma[2] * Apun_pr[i]);",
    "K" = "K[i] = Phi_approx(mu_pr[3] + sigma[3] * K_pr[i]) * 5;",
    "betaF" = "betaF = mu_pr[4] + sigma[4] * betaF_pr[i];",
    "betaP" = "betaP = mu_pr[5] + sigma[5] * betaP_pr[i];"
  )
  non_regression_pars <- setdiff(names(regression_strings), regression_pars)

  main_chunk <- regression_strings[regression_pars]
  main_chunk <- mapply(
    gsub,
    x = main_chunk,
    replacement = as.character(1:length(regression_pars)),
    MoreArgs = list(pattern = "\\{\\{par_index\\}\\}"),
    SIMPLIFY = FALSE
  )
  main_chunk <- c(main_chunk, non_regression_strings[non_regression_pars])
  main_chunk <- main_chunk[sort(names(main_chunk))]

  npars <- length(regression_pars)

  model_code <- glue::glue("
  data {
    int<lower=1> N;
    int<lower=1> T;
    int<lower=1, upper=T> Tsubj[N];
    int choice[N, T];
    real outcome[N, T];
    real sign_out[N, T];

    int<lower=1> ncov;
    matrix[N, ncov] covariate_matrix;
    matrix<lower=0>[{{npars}}, ncov] covariate_precisions;
  }
  transformed data {
    vector[4] initV;
    initV  = rep_vector(0.0, 4);
  }
  parameters {
  // Declare all parameters as vectors for vectorizing
    // Hyper(group)-parameters
    vector[5] mu_pr;
    vector<lower=0>[5] sigma;
    matrix<lower=0>[{{npars}}, ncov] sigma_beta;
    matrix[{{npars}}, ncov] beta;

    // Subject-level raw parameters (for Matt trick)
    vector[N] Arew_pr;
    vector[N] Apun_pr;
    vector[N] K_pr;
    vector[N] betaF_pr;
    vector[N] betaP_pr;
  }
  transformed parameters {
    // Transform subject-level raw parameters
    vector<lower=0, upper=1>[N] Arew;
    vector<lower=0, upper=1>[N] Apun;
    vector<lower=0, upper=5>[N] K;
    vector[N]                   betaF;
    vector[N]                   betaP;

    for (i in 1:N) {
      {{main_chunk$Arew}}
      {{main_chunk$Apun}}
      {{main_chunk$K}}
      {{main_chunk$betaF}}
      {{main_chunk$betaP}}
    }

  }
  model {
    // Hyperparameters
    mu_pr  ~ normal(0, 1);
    sigma[1:3] ~ normal(0, 0.2);
    sigma[4:5] ~ cauchy(0, 1.0);

    // Regression parameters
    for (i in 1:{{npars}}) {
      sigma_beta[i] ~ exponential(covariate_precisions[i]);
      beta[i] ~ normal(0, sigma_beta[i]);
    }

    // individual parameters
    Arew_pr  ~ normal(0, 1.0);
    Apun_pr  ~ normal(0, 1.0);
    K_pr     ~ normal(0, 1.0);
    betaF_pr ~ normal(0, 1.0);
    betaP_pr ~ normal(0, 1.0);

    for (i in 1:N) {
      // Define values
      vector[4] ef;
      vector[4] ev;
      vector[4] PEfreq_fic;
      vector[4] PEval_fic;
      vector[4] pers;   // perseverance
      vector[4] util;

      real PEval;
      real PEfreq;
      real efChosen;
      real evChosen;
      real K_tr;

      // Initialize values
      ef    = initV;
      ev    = initV;
      pers  = initV; // initial pers values
      util  = initV;
      K_tr = pow(3, K[i]) - 1;

      for (t in 1:Tsubj[i]) {
        // softmax choice
        choice[i, t] ~ categorical_logit( util );

        // Prediction error
        PEval  = outcome[i,t] - ev[ choice[i,t]];
        PEfreq = sign_out[i,t] - ef[ choice[i,t]];
        PEfreq_fic = -sign_out[i,t]/3 - ef;

        // store chosen deck ev
        efChosen = ef[ choice[i,t]];
        evChosen = ev[ choice[i,t]];

        if (outcome[i,t] >= 0) {
          // Update ev for all decks
          ef += Apun[i] * PEfreq_fic;
          // Update chosendeck with stored value
          ef[ choice[i,t]] = efChosen + Arew[i] * PEfreq;
          ev[ choice[i,t]] = evChosen + Arew[i] * PEval;
        } else {
          // Update ev for all decks
          ef += Arew[i] * PEfreq_fic;
          // Update chosendeck with stored value
          ef[ choice[i,t]] = efChosen + Apun[i] * PEfreq;
          ev[ choice[i,t]] = evChosen + Apun[i] * PEval;
        }

        // Perseverance updating
        pers[ choice[i,t] ] = 1;   // perseverance term
        pers /= (1 + K_tr);        // decay

        // Utility of expected value and perseverance
        util  = ev + ef * betaF[i] + pers * betaP[i];
      }
    }
  }
  generated quantities {
    // For group level parameters
    real<lower=0,upper=1> mu_Arew;
    real<lower=0,upper=1> mu_Apun;
    real<lower=0,upper=5> mu_K;
    real                  mu_betaF;
    real                  mu_betaP;

    // For log likelihood calculation
    real log_lik[N];

    // For posterior predictive check
    real y_pred[N,T];

    // Set all posterior predictions to -1 (avoids NULL values)
    for (i in 1:N) {
      for (t in 1:T) {
        y_pred[i,t] = -1;
      }
    }

    mu_Arew   = Phi_approx(mu_pr[1]);
    mu_Apun   = Phi_approx(mu_pr[2]);
    mu_K      = Phi_approx(mu_pr[3]) * 5;
    mu_betaF  = mu_pr[4];
    mu_betaP  = mu_pr[5];

    { // local section, this saves time and space
      for (i in 1:N) {
        // Define values
        vector[4] ef;
        vector[4] ev;
        vector[4] PEfreq_fic;
        vector[4] PEval_fic;
        vector[4] pers;   // perseverance
        vector[4] util;

        real PEval;
        real PEfreq;
        real efChosen;
        real evChosen;
        real K_tr;

        // Initialize values
        log_lik[i] = 0;
        ef    = initV;
        ev    = initV;
        pers  = initV; // initial pers values
        util  = initV;
        K_tr = pow(3, K[i]) - 1;

        for (t in 1:Tsubj[i]) {
          // softmax choice
          log_lik[i] += categorical_logit_lpmf( choice[i, t] | util );

          // generate posterior prediction for current trial
          y_pred[i,t] = categorical_rng(softmax(util));

          // Prediction error
          PEval  = outcome[i,t] - ev[ choice[i,t]];
          PEfreq = sign_out[i,t] - ef[ choice[i,t]];
          PEfreq_fic = -sign_out[i,t]/3 - ef;

          // store chosen deck ev
          efChosen = ef[ choice[i,t]];
          evChosen = ev[ choice[i,t]];

          if (outcome[i,t] >= 0) {
            // Update ev for all decks
            ef += Apun[i] * PEfreq_fic;
            // Update chosendeck with stored value
            ef[ choice[i,t]] = efChosen + Arew[i] * PEfreq;
            ev[ choice[i,t]] = evChosen + Arew[i] * PEval;
          } else {
            // Update ev for all decks
            ef += Arew[i] * PEfreq_fic;
            // Update chosendeck with stored value
            ef[ choice[i,t]] = efChosen + Apun[i] * PEfreq;
            ev[ choice[i,t]] = evChosen + Apun[i] * PEval;
          }

          // Perseverance updating
          pers[ choice[i,t] ] = 1;   // perseverance term
          pers /= (1 + K_tr);        // decay

          // Utility of expected value and perseverance
          util  = ev + ef * betaF[i] + pers * betaP[i];
        }
      }
    }
  }", .open="{{", .close="}}", .trim=FALSE)

  return(model_code)
}

write_stan_model_vpp_regression <- function(regression_pars) {
  regression_strings <- list(
    "Arew" = "Arew[i] = Phi_approx(mu_pr[1] + covariate_matrix[i] * beta[{{par_index}}]' + sigma[1] * Arew_pr[i]);",
    "Apun" = "Apun[i] = Phi_approx(mu_pr[2] + covariate_matrix[i] * beta[{{par_index}}]' + sigma[2] * Apun_pr[i]);",
    "K" = "K[i] = Phi_approx(mu_pr[3] + covariate_matrix[i] * beta[{{par_index}}]' + sigma[3] * K_pr[i]) * 5;",
    "betaF" = "betaF = mu_pr[4] + covariate_matrix[i] * beta[{{par_index}}]' + sigma[4] * betaF_pr[i];",
    "betaP" = "betaP = mu_pr[5] + covariate_matrix[i] * beta[{{par_index}}]' + sigma[5] * betaP_pr[i];"
  )
  non_regression_strings <- list(
    "Arew" = "Arew[i] = Phi_approx(mu_pr[1] + sigma[1] * Arew_pr[i]);",
    "Apun" = "Apun[i] = Phi_approx(mu_pr[2] + sigma[2] * Apun_pr[i]);",
    "K" = "K[i] = Phi_approx(mu_pr[3] + sigma[3] * K_pr[i]) * 5;",
    "betaF" = "betaF = mu_pr[4] + sigma[4] * betaF_pr[i];",
    "betaP" = "betaP = mu_pr[5] + sigma[5] * betaP_pr[i];"
  )
  non_regression_pars <- setdiff(names(regression_strings), regression_pars)

  main_chunk <- regression_strings[regression_pars]
  main_chunk <- mapply(
    gsub,
    x = main_chunk,
    replacement = as.character(1:length(regression_pars)),
    MoreArgs = list(pattern = "\\{\\{par_index\\}\\}"),
    SIMPLIFY = FALSE
  )
  main_chunk <- c(main_chunk, non_regression_strings[non_regression_pars])
  main_chunk <- main_chunk[sort(names(main_chunk))]

  npars <- length(regression_pars)

  model_code <- glue::glue("
  data {
    int<lower=1> N;
    int<lower=1> T;
    int<lower=1, upper=T> Tsubj[N];
    int choice[N, T];
    real outcome[N, T];

    int<lower=1> ncov;
    matrix[N, ncov] covariate_matrix;
    matrix<lower=0>[{{npars}}, ncov] covariate_precisions;
  }
  transformed data {
    vector[4] initV;
    initV  = rep_vector(0.0, 4);
  }
  parameters {
  // Declare all parameters as vectors for vectorizing
    // Hyper(group)-parameters
    vector[8] mu_pr;
    vector<lower=0>[8] sigma;
    matrix<lower=0>[{{npars}}, ncov] sigma_beta;
    matrix[{{npars}}, ncov] beta;

    // Subject-level raw parameters (for Matt trick)
    vector[N] A_pr;
    vector[N] alpha_pr;
    vector[N] cons_pr;
    vector[N] lambda_pr;
    vector[N] epP_pr;
    vector[N] epN_pr;
    vector[N] K_pr;
    vector[N] w_pr;
  }
  transformed parameters {
    // Transform subject-level raw parameters
    vector<lower=0, upper=1>[N]  A;
    vector<lower=0, upper=2>[N]  alpha;
    vector<lower=0, upper=5>[N]  cons;
    vector<lower=0, upper=10>[N] lambda;
    vector[N] epP;
    vector[N] epN;
    vector<lower=0, upper=1>[N] K;
    vector<lower=0, upper=1>[N] w;

    for (i in 1:N) {
      {{main_chunk$A}}
      {{main_chunk$alpha}}
      {{main_chunk$cons}}
      {{main_chunk$lambda}}
      {{main_chunk$K}}
      {{main_chunk$w}}
      {{main_chunk$epP}}
      {{main_chunk$epN}}
    }
  }
  model {
    // Hyperparameters
    mu_pr       ~ normal(0, 1.0);
    sigma[1:4] ~ normal(0, 0.2);
    sigma[5:6] ~ cauchy(0, 1.0);
    sigma[7:8] ~ normal(0, 0.2);

    // Regression parameters
    for (i in 1:{{npars}}) {
      sigma_beta[i] ~ exponential(covariate_precisions[i]);
      beta[i] ~ normal(0, sigma_beta[i]);
    }

    // individual parameters
    A_pr      ~ normal(0, 1.0);
    alpha_pr  ~ normal(0, 1.0);
    cons_pr   ~ normal(0, 1.0);
    lambda_pr ~ normal(0, 1.0);
    epP_pr    ~ normal(0, 1.0);
    epN_pr    ~ normal(0, 1.0);
    K_pr      ~ normal(0, 1.0);
    w_pr      ~ normal(0, 1.0);

    for (i in 1:N) {
      // Define values
      vector[4] ev;
      vector[4] p_next;
      vector[4] str;
      vector[4] pers;   // perseverance
      vector[4] V;   // weighted sum of ev and pers

      real curUtil;     // utility of curFb
      real theta;       // theta = 3^c - 1

      // Initialize values
      theta = pow(3, cons[i]) -1;
      ev    = initV; // initial ev values
      pers  = initV; // initial pers values
      V     = initV;

      for (t in 1:Tsubj[i]) {
        // softmax choice
        choice[i, t] ~ categorical_logit(theta * V);

        // perseverance decay
        pers *= K[i]; // decay

        if (outcome[i, t] >= 0) {  // x(t) >= 0
          curUtil = pow(outcome[i, t], alpha[i]);
          pers[choice[i, t]] += epP[i];  // perseverance term
        } else {                  // x(t) < 0
          curUtil = -1 * lambda[i] * pow(-1 * outcome[i, t], alpha[i]);
          pers[choice[i, t]] += epN[i];  // perseverance term
        }

        ev[choice[i, t]] += A[i] * (curUtil - ev[choice[i, t]]);
        // calculate V
        V = w[i] * ev + (1-w[i]) * pers;
      }
    }
  }
  generated quantities {
    // For group level parameters
    real<lower=0, upper=1>  mu_A;
    real<lower=0, upper=2>  mu_alpha;
    real<lower=0, upper=5>  mu_cons;
    real<lower=0, upper=10> mu_lambda;
    real mu_epP;
    real mu_epN;
    real<lower=0, upper=1> mu_K;
    real<lower=0, upper=1> mu_w;

    // For log likelihood calculation
    real log_lik[N];

    // For posterior predictive check
    real y_pred[N, T];

    // Set all posterior predictions to 0 (avoids NULL values)
    for (i in 1:N) {
      for (t in 1:T) {
        y_pred[i, t] = -1;
      }
    }

    mu_A      = Phi_approx(mu_pr[1]);
    mu_alpha  = Phi_approx(mu_pr[2]) * 2;
    mu_cons   = Phi_approx(mu_pr[3]) * 5;
    mu_lambda = Phi_approx(mu_pr[4]) * 10;
    mu_epP    = mu_pr[5];
    mu_epN    = mu_pr[6];
    mu_K      = Phi_approx(mu_pr[7]);
    mu_w      = Phi_approx(mu_pr[8]);

    { // local section, this saves time and space
      for (i in 1:N) {
        // Define values
        vector[4] ev;
        vector[4] p_next;
        vector[4] str;
        vector[4] pers;   // perseverance
        vector[4] V;   // weighted sum of ev and pers

        real curUtil;     // utility of curFb
        real theta;       // theta = 3^c - 1

        // Initialize values
        log_lik[i] = 0;
        theta      = pow(3, cons[i]) -1;
        ev         = initV; // initial ev values
        pers       = initV; // initial pers values
        V          = initV;

        for (t in 1:Tsubj[i]) {
          // softmax choice
          log_lik[i] += categorical_logit_lpmf(choice[i, t] | theta * V);

          // generate posterior prediction for current trial
          y_pred[i, t] = categorical_rng(softmax(theta * V));

          // perseverance decay
          pers *= K[i]; // decay

          if (outcome[i, t] >= 0) {  // x(t) >= 0
            curUtil = pow(outcome[i, t], alpha[i]);
            pers[choice[i, t]] += epP[i];  // perseverance term
          } else {                  // x(t) < 0
            curUtil = -1 * lambda[i] * pow(-1 * outcome[i, t], alpha[i]);
            pers[choice[i, t]] += epN[i];  // perseverance term
          }

          ev[choice[i, t]] += A[i] * (curUtil - ev[choice[i, t]]);
          // calculate V
          V = w[i] * ev + (1-w[i]) * pers;
        }
      }
    }
  }", .open="{{", .close="}}", .trim=FALSE)

  return(model_code)
}
