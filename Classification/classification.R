df_viro <- df1[!is.na(df1$vironostika), ]
df_viro <- df_viro[!is.na(df_viro$single_multiple_lineage), ]
df_viro$Y <- ifelse(df_viro$TSI <= 180, 1, 0)
df_viro$vironostika_scaled <- as.numeric(scale(df_viro$vironostika))
df_viro$M <- ifelse(df_viro$single_multiple_lineage=='multiple', 1, 0)

id <- unique(df_viro$ID)
set.seed(123)
shuffled_id <- sample(id)
test_id <- shuffled_id[1:7]
train_id <- shuffled_id[8:length(shuffled_id)]
df_viro_train <- df_viro %>% filter(ID %in% train_id)
df_viro_test  <- df_viro %>% filter(ID %in% test_id)
df_viro_test <- df_viro_test %>%
  mutate(id_id = as.integer(as.factor(ID)))
df_viro_train <- df_viro_train %>%
  mutate(id_id = as.integer(as.factor(ID)))

viro_txt <- "
data{
  int<lower=1> N; // number of observations
  array[N] int<lower=0> Y;
  vector[N] M; // binary: 1 iff multiple
  vector[N] vironostika;
  int<lower=1> N_tilde; // number of observations
  array[N_tilde] int<lower=0> Y_tilde;
  vector[N_tilde] M_tilde; // binary: 1 iff multiple
  vector[N_tilde] vironostika_tilde;
}
parameters{
  real alpha;
  real beta_M;
  real beta_vironostika;
}
transformed parameters {
  vector[N] mu;
  
  mu = alpha + beta_M * M + beta_vironostika * vironostika;
}
model {
  // Priors
  alpha ~ normal(0, 1);
  beta_M ~ normal(0, 1);
  beta_vironostika ~ normal(0, 1);

  // Likelihood (logistic regression)
  Y ~ bernoulli_logit(alpha + beta_M * M + beta_vironostika * vironostika);
}
generated quantities {
  vector[N_tilde] p_tilde;
  
  for (n in 1:N_tilde) {
    real eta = alpha + beta_M * M_tilde[n] + beta_vironostika * vironostika_tilde[n];
    p_tilde[n] = inv_logit(eta);
  }
}
"

# compile the model
viro_filename <- cmdstanr::write_stan_file(
  gsub('\t',' ', viro_txt ),
  dir = out,
  basename = NULL,
  force_overwrite = FALSE,
  hash_salt = ""
)

# compile Stan model
viro_compiled <- cmdstanr::cmdstan_model(viro_filename)

# data
stan_data <- list(
  N = nrow(df_viro_train),
  Y = df_viro_train$Y,
  M = df_viro_train$M,
  vironostika = df_viro_train$vironostika_scaled,
  N_tilde = nrow(df_viro_test),
  Y_tilde = df_viro_test$Y,
  M_tilde = df_viro_test$M,
  vironostika_tilde = df_viro_test$vironostika_scaled
)

# sample
viro_fit <- viro_compiled$sample(
  data = stan_data,
  seed = 6373,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 2000,
  refresh = 500
)

viro_fit$save_object(file = file.path(out, "HR_viro.rds"))

p_tilde <- viro_fit$draws(variables = "p_tilde", 
                                     inc_warmup = FALSE,
                                     format = "draws_df"
)
p_mean <- colMeans(p_tilde)[1:46]
Y_pred <- ifelse(p_mean > 0.5, 1, 0)
Y_true <- df_viro_test$Y

accuracy <- mean(Y_pred == Y_true)
precision <- sum(Y_pred == 1 & Y_true == 1) / sum(Y_pred == 1)
recall <- sum(Y_pred == 1 & Y_true == 1) / sum(Y_true == 1)

c(Accuracy = accuracy, Precision = precision, Recall = recall)
