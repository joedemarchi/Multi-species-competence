
data {
  
  int<lower = 1> N; // Number of data points
  int<lower = 1> M; // Missing loads
  int<lower = 1> P; // Missing sheds
  array[N] real load;
  array[N] real shedding;
  array[N] int missing_ind;
  array[N] int missing_shed_ind;
  
} parameters {
  
  // Shedding relationship parameters
  real beta0; // The shedding rate
  real beta1; // Would expect to be one for proportional shedding
  real<lower=0> sigma; // 
  
  // Missing value parameters
  array[M] real missing_logload;
  array[P] real missing_logshed;

} transformed parameters {
  
  // Hold the systematic component for the swad vs. shedding relationship
  array[N] real mu;
  
  // Set up the systematic components for missing covariate data
  for(i in 1:N){
    
    if(load[i] > 0){
    
      // Swab data present
      mu[i] = beta0 + beta1*log(load[i]);
      
    } else {
      
      // Missing swab data
      mu[i] = beta0 + beta1*missing_logload[missing_ind[i]];
    }
    
  }
    
} model {
  
  // Priors
  beta0 ~ normal(0, 5);
  beta1 ~ normal(0, 5);
  missing_logload ~ normal(0, 5);
  missing_logshed ~ normal(0, 5);
  sigma ~ exponential(1);
  
  // Loop over data points
  for(i in 1:N){
    
    if(shedding[i] != 0){
      
      // Not missing shed data
      target += normal_lpdf(log(shedding[i]) | mu[i], sigma);
      
    } else{
    
      // Account for missing shed dat 
      target += normal_lpdf(missing_logshed[missing_shed_ind[i]] | mu[i], sigma);
      
    }

  }
}

