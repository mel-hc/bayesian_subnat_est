stan_script <- "

///INPUTS//
////////////////////////

data {
int<lower=1>            N_obs;               // Number of observations (rows of data)
int<lower=1>            N_geog;              // Number of municipalities
int<lower=1>            N_time;              // Number of years
int                     y_notif[N_obs];      // Vector of notification counts (N_notif_it)
int                     y_all_mort[N_obs];   // Vector of total death counts (N_mort_it) 
int                     x_id[N_obs];         // municipality or state ID (as continuous intergers)
int                     x_year[N_obs];       // Vector of years (indexed from 1)
vector[N_obs]           fhs;                 // FHS teams per 4000 pop, truncated (std)
vector[N_obs]           pop_100k;            // population_it, in 100Ks (std)
vector[N_obs]           ln_gdp;              // ln_gdp per capita in 1000s of $R, 2015 (std)
vector[N_obs]           p_cov;               // fraction coverage, by it 
vector[N_obs]           idc;                 // ill/poorly defined cause (fraction all)
vector[N_obs]           mort_treat;          // Vector of treatment deaths 
vector[N_obs]           aban_treat;          // Vector of treatment abandonment 
}

///////////////////////////////////////////

parameters {
real                        b_inc_00; 
vector[N_geog]              b_inc_i0;  
matrix[N_geog, N_time-1 ]   b_inc_ij; 
real                        b_inc_gdp;
real                        b_inc_fhs; 
real<lower=0>               sigma_inc_i;  
real<lower=0>               sigma_inc_ij; 

real                        b_pn_00; 
vector[N_geog]              b_pn_i0;  
matrix[N_geog, N_time-1 ]   b_pn_ij;
real                        b_pn_fhs; 
real                        b_pn_gdp;
real<lower=0>               sigma_pn_i;
real<lower=0>               sigma_pn_ij;

real                        b_adj_00; 
real                        b_adj_0j; 
vector[N_geog]              b_adj_i0;  
real                        b_adj_3;
real                        y;
real                        z;
real<lower=0>               sigma_adj_i;

real<lower=0, upper=1>      p_surv_no_notif;
real<lower=0, upper=1>      p_mort_abandon;
}
/////////////////////////////////////////////////

transformed parameters {

vector[N_obs]               inc;
matrix[N_geog, N_time]      b_inc_tmp; // a temporary container  
matrix[N_geog, N_time]      b_inc;
real                        b_inc_i0_mean; // for the mean   
vector[N_geog]              b_inc_ij_mean; // for the mean. 
vector[N_obs]               pn;
matrix[N_geog, N_time]      b_pn_tmp; // a temporary container  
matrix[N_geog, N_time]      b_pn;
real                        b_pn_i0_mean; // for the mean   
vector[N_geog]              b_pn_ij_mean; // for the mean. 

vector[N_obs]               death_adj; 
vector[N_geog]              state_adj; 
vector[N_time]              time_adj; 
real                        mean_adj; 
real<lower=0, upper=1>      a; 
real<lower=0, upper=1>      b;
vector[N_obs]               m_deaths;
vector[N_obs]               p_mort; 

vector[N_obs]               cases;
vector[N_obs]               mc_count;
vector[N_obs]               mort_notif;
vector[N_obs]               mort_no_treat;


//// INCIDENCE
for(i in 1:N_geog){ 
  b_inc_tmp[i,1] = 0.0; 
  for(j in 2:N_time){ 
     b_inc_tmp[i,j] = b_inc_tmp[i,j-1] + b_inc_ij[i,j-1]; // random walk started from zero
  } 
} 

b_inc_i0_mean = mean(b_inc_i0); // mean of the state random effects (to subtract off later)
for(i in 1:N_geog){ 
  b_inc_ij_mean[i] = mean(b_inc_tmp [i,]); // mean of the state-year 
  //random effects (by state), to subtract off later
} 

for(i in 1:N_geog){ 
  for(j in 1:N_time){ 
    b_inc[i,j] = b_inc_00 + b_inc_i0[i] - b_inc_i0_mean + b_inc_tmp[i,j] - 
    b_inc_ij_mean[i]; 
  } //~~
} //~~


for(i in 1:N_obs){
inc[i]      = exp(b_inc[ x_id[i], x_year[i] ] + b_inc_fhs*fhs[i] + 
b_inc_gdp*ln_gdp[i]);  
}

//// FRACTION TREATED
for(i in 1:N_geog){ 
  b_pn_tmp[i,1] = 0.0; 
  for(j in 2:N_time){ 
     b_pn_tmp[i,j] = b_pn_tmp[i,j-1] + b_pn_ij[i,j-1]*sigma_pn_ij; 
  }
}

b_pn_i0_mean = mean(b_pn_i0*sigma_pn_i);
for(i in 1:N_geog){ 
  b_pn_ij_mean[i] = mean(b_pn_tmp [i,]); 
}


for(i in 1:N_geog){ 
  for(j in 1:N_time){ 
    b_pn[i,j] = b_pn_00 + b_pn_i0[i]*sigma_pn_i - b_pn_i0_mean + b_pn_tmp[i,j] - 
    b_pn_ij_mean[i]; 
  } 
} 

for(i in 1:N_obs){
pn[i]       = inv_logit((b_pn[ x_id[i], x_year[i] ]) + b_pn_fhs*fhs[i] + 
b_pn_gdp*ln_gdp[i]); 
}

//// DEATH ADJUSTMENT
for(i in 1:N_time){
  time_adj[i] = b_adj_0j*(x_year[i] - 10);
}
for(i in 1:N_geog){
  state_adj[i] = b_adj_00 + b_adj_i0[i]*sigma_adj_i; 
}
mean_adj = b_adj_00 + mean(b_adj_i0*sigma_adj_i);

a   = inv_logit(mean_adj + b_adj_3*y); 
b   = inv_logit(mean_adj + b_adj_3*z);

for(i in 1:N_obs){
    death_adj[i] = inv_logit(state_adj[ x_id[i] ] + time_adj[ x_year[i] ] + 
    b_adj_3*idc[i]);  
}
//// PR DEATH GIVEN NOTIF
for(i in 1:N_obs){ 
  p_mort[i] = mort_treat[i] + p_mort_abandon*aban_treat[i];  
 } 

// MODELED DEATHS
for(i in 1:N_obs){
m_deaths[i] = pop_100k[i] * inc[i] * (pn[i] * p_mort[i] + ((1-pn[i]) * 
(1-p_surv_no_notif))); 
}

 // Outcomes to make life easier
       for(i in 1:N_obs){ 
         cases[i] = inc[i] *  pop_100k[i];
      }
       
      for(i in 1:N_obs){ 
          mc_count[i] = pop_100k[i] * inc[i] *  (1 - pn[i]);
      }
      
      for(i in 1:N_obs){
          mort_notif[i] = p_mort[x_id[i]] * y_notif[i];
      }
      
      for(i in 1:N_obs){
          mort_no_treat[i] = m_deaths[i] - (p_mort[x_id[i]] * y_notif[i]);
}
}
/////////////////////////////////////////////////

model {

to_vector(b_inc_ij)     ~ normal(0, sigma_inc_ij);  
b_inc_00                ~ normal(0, 10); 
b_inc_i0                ~ normal(0, sigma_inc_i); 
b_inc_gdp               ~ normal(0, 10);
b_inc_fhs               ~ normal(0, 10); 
sigma_inc_i             ~ cauchy(0, 2); 
sigma_inc_ij            ~ cauchy(0, 2); 

to_vector(b_pn_ij)      ~ normal(0, 1);  
b_pn_00                 ~ normal(0, 10); 
b_pn_i0                 ~ normal(0, 1); 
b_pn_fhs                ~ normal(0, 10); 
b_pn_gdp                ~ normal(0, 10);
sigma_pn_i              ~ cauchy(0, 2);  
sigma_pn_ij             ~ cauchy(0, 2);  

b_adj_00                ~ normal(0, 1); 
b_adj_i0                ~ normal(0, 1);
sigma_adj_i             ~ cauchy(0, 2);  
b_adj_0j                ~ normal(0, 0.1); //normal(0, 0.05); note much tighter prior 
b_adj_3                 ~ normal(0, 1); 
y                       ~ normal(0.01, .001); // state: 0.001 sd
z                       ~ normal(0.15, .001); // state: 0.001 sd
a                       ~  beta(52.9737, 451.15377); // from expert survey 
b                       ~  beta(97.82789, 285.81089); // from expert survey

p_surv_no_notif         ~ beta(25.65, 33.32); // from expert survey
p_mort_abandon          ~ beta(4.287894, 81.469979); // lowered to 0.01 - 0.1, mean 0.05 

//// LIKELIHOODS
for(i in 1:N_obs){
 y_notif[i]      ~ poisson(pop_100k[i] * inc[i] * pn[i]); 
}

for(i in 1:N_obs){
y_all_mort[i]   ~ poisson(m_deaths[i] * p_cov[i] * (1-death_adj[i]));
}
}
"

####################################################