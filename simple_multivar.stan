data {
  int<lower = 1> N; //number of genes
  int<lower = 1> J; //number of experiments
  matrix[N,J] y; //data;
} parameters {
  vector<lower = 0, upper = 1>[N] theta; //mixing proportions
  vector[N] mu0; //locations of mixture components
  vector[N] mu1;
  real<lower = 0> alpha;
  vector<lower = 0>[N] sigma;
  real<lower = 0, upper = 10> aleph;
  real<lower = 0, upper = 10> tau;

} model {
  sigma ~ lognormal(0, alpha);
  
  alpha ~ normal(0, .5);
  
  aleph ~ uniform(0,10);
  tau ~ uniform(0,10);
  
  mu0 ~ normal(0, .05);
  mu1 ~ normal(0, 2);
  for(n in 1:N) {

    theta[n] ~ beta(aleph, tau);


    for(j in 1:J) {
      target += log_mix(theta[n], normal_lpdf(y[n,j] | mu0[n], sigma[n]),
                        normal_lpdf(y[n,j] |  mu0[n] + mu1[n], sigma[n]));
    }
  }
}

