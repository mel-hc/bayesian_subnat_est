stan_script <- "

data {
int<lower=1>            N_obs;               // Number of observations (rows of data)
int<lower=1>            N_geog;              // Number of municipalities
int<lower=1>            N_time;              // Number of years
int                     y_notif[N_obs];      // Vector of notification counts 
int                     y_all_mort[N_obs];   // Vector of total death counts 
int                     x_id[N_obs];         // Municipality or state ID (as continuous intergers)
int                     x_year[N_obs];       // Vector of years (indexed from 1)
vector[N_obs]           fhs;                 // FHS teams per 4000 pop (std)
vector[N_obs]           pop_100k;            // Population in 100Ks 
vector[N_obs]           ln_gdp;              // ln(GDP) per capita in 1000s of $R (std)
vector[N_obs]           p_cov;               // Fraction of mortality system coverage 
vector[N_obs]           idc;                 // Poorly defined cause (fraction of all mortality)
vector[N_obs]           mort_treat;          // Vector of fitted treatment death probability
vector[N_obs]           aban_treat;          // Vector of fitted treatment abandonment probability
}

parameters {
// INCIDENCE
real                        b_inc_00; // intercept
vector[N_geog]              b_inc_i0; // place random effect
matrix[N_geog, N_time-1 ]   b_inc_ij; // place/time random effect
real                        b_inc_gdp; // covariate 
real                        b_inc_fhs; // covariate 
real<lower=0>               sigma_inc_i; // sigma place effect
real<lower=0>               sigma_inc_ij; // sigma place/time effect
// FRACTION TREATED
real                        b_ft_00; // intercept
vector[N_geog]              b_ft_i0; // place random effect
matrix[N_geog, N_time-1 ]   b_ft_ij; // place/time random effect
real                        b_ft_fhs; // covariate
real                        b_ft_gdp; // covariate
real<lower=0>               sigma_ft_i; // sigma place effect
real<lower=0>               sigma_ft_ij; // sigma place/time effect
// DEATH ADJUSTMENT
real                        b_adj_00; // intercept
real                        b_adj_0j; // time effect
vector[N_geog]              b_adj_i0; // place effect
real                        b_adj_3; // coeffecient for poorly defined cause
real                        y; // from expert opinion 
real                        z; // from expert opinion 
real<lower=0>               sigma_adj_i; // sigma place effect 

real<lower=0, upper=1>      p_surv_no_notif; // probability survival without treatment
real<lower=0, upper=1>      p_mort_abandon; // probability of death given treatment abandoned
}

transformed parameters {
// INCIDENCE
vector[N_obs]               inc; 
matrix[N_geog, N_time]      b_inc_tmp; // a temporary container  
matrix[N_geog, N_time]      b_inc; 
real                        b_inc_i0_mean; // for the mean   
vector[N_geog]              b_inc_ij_mean; // for the mean
// FRACTION TREATED
vector[N_obs]               ft;
matrix[N_geog, N_time]      b_ft_tmp; // a temporary container  
matrix[N_geog, N_time]      b_ft;
real                        b_ft_i0_mean; // for the mean   
vector[N_geog]              b_ft_ij_mean; // for the mean. 
// DEATH ADJUSTMENT
vector[N_obs]               death_adj; 
vector[N_geog]              place_adj; 
vector[N_time]              time_adj; 
real                        mean_adj; 
real<lower=0, upper=1>      a; 
real<lower=0, upper=1>      b;
// MODELED DEATHS
vector[N_obs]               m_deaths;
vector[N_obs]               p_mort; 

//// INCIDENCE
for(i in 1:N_geog){ 
  b_inc_tmp[i,1] = 0.0; 
  for(j in 2:N_time){ 
     b_inc_tmp[i,j] = b_inc_tmp[i,j-1] + b_inc_ij[i,j-1]; // random walk started from zero
  } 
} 

b_inc_i0_mean = mean(b_inc_i0); // mean of the place random effects (to subtract off later)
for(i in 1:N_geog){ 
  b_inc_ij_mean[i] = mean(b_inc_tmp [i,]); // mean of the place-year 
  //random effects (by place), to subtract off later
} 

for(i in 1:N_geog){ 
  for(j in 1:N_time){ 
    b_inc[i,j] = b_inc_00 + b_inc_i0[i] - b_inc_i0_mean + b_inc_tmp[i,j] - 
    b_inc_ij_mean[i]; 
  } 
} 

for(i in 1:N_obs){
inc[i]      = exp(b_inc[ x_id[i], x_year[i] ] + b_inc_fhs*fhs[i] + b_inc_gdp*ln_gdp[i]);  
}

//// FRACTION TREATED
for(i in 1:N_geog){ 
  b_ft_tmp[i,1] = 0.0; 
  for(j in 2:N_time){ 
     b_ft_tmp[i,j] = b_ft_tmp[i,j-1] + b_ft_ij[i,j-1]*sigma_ft_ij; 
  }
}

b_ft_i0_mean = mean(b_ft_i0*sigma_ft_i);
for(i in 1:N_geog){ 
  b_ft_ij_mean[i] = mean(b_ft_tmp [i,]); 
}


for(i in 1:N_geog){ 
  for(j in 1:N_time){ 
    b_ft[i,j] = b_ft_00 + b_ft_i0[i]*sigma_ft_i - b_ft_i0_mean + b_ft_tmp[i,j] - 
    b_ft_ij_mean[i]; 
  } 
} 

for(i in 1:N_obs){
ft[i]       = inv_logit((b_ft[ x_id[i], x_year[i] ]) + b_ft_fhs*fhs[i] + b_ft_gdp*ln_gdp[i]); 
}

//// DEATH ADJUSTMENT
for(i in 1:N_time){
  time_adj[i] = b_adj_0j*(x_year[i] - 10);
}
for(i in 1:N_geog){
  place_adj[i] = b_adj_00 + b_adj_i0[i]*sigma_adj_i; 
}
mean_adj = b_adj_00 + mean(b_adj_i0*sigma_adj_i);

a   = inv_logit(mean_adj + b_adj_3*y); 
b   = inv_logit(mean_adj + b_adj_3*z);

for(i in 1:N_obs){
    death_adj[i] = inv_logit(place_adj[ x_id[i] ] + time_adj[ x_year[i] ] + b_adj_3*idc[i]);  
}
//// PR DEATH GIVEN NOTIFICATION
for(i in 1:N_obs){ 
  p_mort[i] = mort_treat[i] + p_mort_abandon*aban_treat[i];  
 } 

// MODELED DEATHS
for(i in 1:N_obs){
m_deaths[i] = pop_100k[i] * inc[i] * (pn[i] * p_mort[i] + ((1-pn[i]) * (1-p_surv_no_notif))); 
}
}

model {

to_vector(b_inc_ij)     ~ normal(0, sigma_inc_ij);  
b_inc_00                ~ normal(0, 10); 
b_inc_i0                ~ normal(0, sigma_inc_i); 
b_inc_gdp               ~ normal(0, 10);
b_inc_fhs               ~ normal(0, 10); 
sigma_inc_i             ~ cauchy(0, 2); 
sigma_inc_ij            ~ cauchy(0, 2); 

to_vector(b_ft_ij)      ~ normal(0, 1);  
b_ft_00                 ~ normal(0, 10); 
b_ft_i0                 ~ normal(0, 1); 
b_ft_fhs                ~ normal(0, 10); 
b_ft_gdp                ~ normal(0, 10);
sigma_ft_i              ~ cauchy(0, 2);  
sigma_ft_ij             ~ cauchy(0, 2);  

b_adj_00                ~ normal(0, 1); 
b_adj_i0                ~ normal(0, 1);
sigma_adj_i             ~ cauchy(0, 2);  
b_adj_0j                ~ normal(0, 0.05); 
b_adj_3                 ~ normal(0, 1); 
y                       ~ normal(0.01, .001);
z                       ~ normal(0.15, .001); 
a                       ~  beta(52.9737, 451.15377); // from expert survey 
b                       ~  beta(97.82789, 285.81089); // from expert survey

p_surv_no_notif         ~ beta(25.65, 33.32); // from expert survey
p_mort_abandon          ~ beta(4.287894, 81.469979); // 0.01 - 0.1, mean 0.05 

//// LIKELIHOODS
for(i in 1:N_obs){
 y_notif[i]      ~ poisson(pop_100k[i] * inc[i] * ft[i]); 
}

for(i in 1:N_obs){
y_all_mort[i]   ~ poisson(m_deaths[i] * p_cov[i] * (1-death_adj[i]));
}
}
"

####################################################
