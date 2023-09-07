

data {
  int K;
  int D;
  int N;
  int bin_size;
  real generation_time;
  matrix[N,D] X;
  int Y[N,K];
  int Y_sum[N];
}

transformed data {
  vector[D] Zeros;
  real Max_date;
  matrix[N,D] X_norm;
  Zeros = rep_vector(0,D);
  Max_date = X[N,2];
  for (n in 1:N) {
    X_norm[n,1] = 1;
    X_norm[n,2] = X[n,2] / Max_date;
  }
}

parameters {
  matrix[D,K-1] b_raw;
}

transformed parameters {
  matrix[D,K] b;
  matrix[N,K] mu;
  b = append_col(Zeros, b_raw);
  mu = X_norm*b;
}

model {
  for (n in 1:N)
    Y[n,] ~ multinomial_logit(mu[n,]');
}

generated quantities {
  vector[K-1] growth_rate;
  matrix[N,K] theta;
  int Y_predict[N,K];
  
  for(k in 1:(K-1)){
      growth_rate[k] = exp(((b_raw[2,k] / Max_date) / bin_size)  * generation_time);
  }
  
  for(n in 1:N){
    theta[n,] = softmax(mu[n,]')';
    Y_predict[n,] = multinomial_rng(softmax(mu[n,]'),Y_sum[n]);
  }

}

