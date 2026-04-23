
data {
  int<lower=1> N;          // number of observations
  int<lower=1> T;          // number of distinct sampling times
  array[N] int<lower=1, upper=T> t_index; // time index for each obs (vector 1: 4)
  vector[N] y;             // observed log loads
}

parameters {
  real a;                  // intercept (per time step)
  real b;                  // slope (per time step) how strongly today's load                                      //influences the next time step
  real<lower=0> sigma;     // process error - variability in growth between steps
  real m1;                 // population mean log load at first sampling time
  real<lower=0> v1;        // population variance of log-load at first sampling
}

transformed parameters {    // Compute latent pop mean and variance at each step
  vector[T] m;             // mean at each time
  vector<lower=0>[T] v;    // variance at each time
  
  m[1] = m1;
  v[1] = v1;
  
  for (t in 2:T) {
    m[t] = a + b * m[t-1];
    v[t] = square(b) * v[t-1] + square(sigma);
  }
}

model {
  // Priors
  a ~ normal(0, 2);     // can be positive or neg but centered near 0 on log scale
  b ~ normal(1, 1);        // centered near 1 meaning today is close to tomorrow
  sigma ~ exponential(1);   // favors small process noise
  m1 ~ normal(0, 5);      // wide prior on initial mean load
  v1 ~ exponential(1);  // positive initial variance
  
  // Likelihood: each observation comes from population at its time
  for (i in 1:N) {
    y[i] ~ normal(m[t_index[i]], sqrt(v[t_index[i]]));
  }
}

