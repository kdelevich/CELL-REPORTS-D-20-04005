data{
  int<lower=1> n_s;                               //number of subjects
  int<lower=1> n_g;                               //number of subjects
  int<lower=1> n_t[n_s+1];                               //number of trials per subject, plus total number of trials
  int<lower=1,upper=4> Choice[n_t[n_s+1]];           // choice options trial n_t. (All subjects stacked)
  int<lower=0,upper=1> Reward[n_t[n_s+1]];           // reward (=1, yes)? trial n_t
  int<lower=1,upper=n_s> Subject[n_t[n_s+1]];        // subject number
  int<lower=1,upper=n_g> Group[n_s];        // group number
  int<lower=0,upper=1> Init[n_t[n_s+1]];           // is this first trial of a subject? Should RL be initialized?
  int<lower=1,upper=2> Phase[n_t[n_s+1]];           // is this a phase change trial?
}// end data
  

parameters{
  // full group level mean parameters
	real<lower=0, upper=10> mua_alpha_b1; 				  //inverse gain parameter
	real<lower=0, upper=10> mub_alpha_b1; 				  //inverse gain parameter
	real<lower=0, upper=10> mua_beta_b1; 				  //inverse gain parameter
	real<lower=0, upper=10> mub_beta_b1; 				  //inverse gain parameter
	real<lower=0, upper=10> mua_alpha_b2; 				  //inverse gain parameter
	real<lower=0, upper=10> mub_alpha_b2; 				  //inverse gain parameter
	real<lower=0, upper=10> mua_beta_b2; 				  //inverse gain parameter
	real<lower=0, upper=10> mub_beta_b2; 				  //inverse gain parameter
	real<lower=0, upper=10> mua_alpha_a1; 				  //inverse gain parameter
	real<lower=0, upper=10> mub_alpha_a1; 				  //inverse gain parameter
	real<lower=0, upper=10> mua_beta_a1; 				  //inverse gain parameter
	real<lower=0, upper=10> mub_beta_a1; 				  //inverse gain parameter
	real<lower=0, upper=10> mua_alpha_a2; 				  //inverse gain parameter
	real<lower=0, upper=10> mub_alpha_a2; 				  //inverse gain parameter
	real<lower=0, upper=10> mua_beta_a2; 				  //inverse gain parameter
	real<lower=0, upper=10> mub_beta_a2; 				  //inverse gain parameter
    

  // group level parameters
  real<lower=0> alpha_b1_gr[n_g];           //inverse gain parameter
  real<lower=0> beta_b1_gr[n_g];           //inverse gain parameter
  real<lower=0> alpha_b2_gr[n_g];           //inverse gain parameter
  real<lower=0> beta_b2_gr[n_g];           //inverse gain parameter
  real<lower=0> alpha_a1_gr[n_g];           //inverse gain parameter
  real<lower=0> beta_a1_gr[n_g];           //inverse gain parameter
  real<lower=0> alpha_a2_gr[n_g];           //inverse gain parameter
  real<lower=0> beta_a2_gr[n_g];           //inverse gain parameter

  // individual level parameters
  real<lower=0, upper=40> b1_ind[n_s];   			  //inverse gain parameter
  real<lower=0, upper=40> b2_ind[n_s];   			  //inverse gain parameter
	real<lower=0, upper=1> a1_ind[n_s];   				//alphaG
	real<lower=0, upper=1> a2_ind[n_s];   				//alphaG

    // group level odor preference
    real<lower=0,upper=1> Q0[4];
}//end paramters
	
model{
  // define general variables needed for subject loop
  int si;
  vector[4] Q;
  int a;
  real alpha;
  real beta;
  vector[4] pchoice;
  real epsilon;
  epsilon <- .00001;

  // set prior on full group level mean parameters
  mua_alpha_b1 ~  uniform(0,10);   			  //inverse gain parameter
  mub_alpha_b1 ~ uniform(0,10);   				//alphaG
  mua_beta_b1 ~  uniform(0,10);   			  //inverse gain parameter
  mub_beta_b1 ~ uniform(0,10);   				//alphaG
  mua_alpha_b2 ~  uniform(0,10);   			  //inverse gain parameter
  mub_alpha_b2 ~ uniform(0,10);   				//alphaG
  mua_beta_b2 ~  uniform(0,10);   			  //inverse gain parameter
  mub_beta_b2 ~ uniform(0,10);   				//alphaG
  mua_alpha_a1 ~  uniform(0,10);   			  //inverse gain parameter
  mub_alpha_a1 ~ uniform(0,10);   				//alphaG
  mua_beta_a1 ~  uniform(0,10);   			  //inverse gain parameter
  mub_beta_a1 ~ uniform(0,10);   				//alphaG
  mua_alpha_a2 ~  uniform(0,10);   			  //inverse gain parameter
  mub_alpha_a2 ~ uniform(0,10);   				//alphaG
  mua_beta_a2 ~  uniform(0,10);   			  //inverse gain parameter
  mub_beta_a2 ~ uniform(0,10);   				//alphaG
 
  
  // set prior on group level odor preferences
  for (o in 1:4){
        Q0[o] ~ uniform(0,1);
  }

  // set prior for group level parameters
  for (g in 1:n_g)
  {
   alpha_b1_gr[g] ~ gamma(mua_alpha_b1,mub_alpha_b1);            //inverse gain parameter
    beta_b1_gr[g] ~  gamma(mua_beta_b1,mub_beta_b1);          //inverse gain parameter
   alpha_b2_gr[g] ~ gamma(mua_alpha_b2,mub_alpha_b2);            //inverse gain parameter
    beta_b2_gr[g] ~  gamma(mua_beta_b2,mub_beta_b2);          //inverse gain parameter
    alpha_a1_gr[g]~ gamma(mua_alpha_a1,mub_alpha_a1);           //alphaG
    beta_a1_gr[g] ~ gamma(mua_beta_a1,mub_beta_a1);          //alphaG
    alpha_a2_gr[g]~ gamma(mua_alpha_a2,mub_alpha_a2);           //alphaG
    beta_a2_gr[g] ~ gamma(mua_beta_a2,mub_beta_a2);          //alphaG
  }

  // set prior for individual level parameters
  for (s in 1:n_s)
  {
    b1_ind[s] ~ gamma(alpha_b1_gr[Group[s]],   beta_b1_gr[Group[s]]);      //inverse gain parameter
    b2_ind[s] ~ gamma(alpha_b2_gr[Group[s]],   beta_b2_gr[Group[s]]);      //inverse gain parameter
	  a1_ind[s]~ beta(alpha_a1_gr[Group[s]],  beta_a1_gr[Group[s]]);   				//alphaG
	  a2_ind[s]~ beta(alpha_a2_gr[Group[s]],  beta_a2_gr[Group[s]]);   				//alphaG
  }
  //print("log-posterior = ", get_lp());

  // now start looping over subjects
  for (t in 1:n_t[n_s+1])
  {

      // set initial values subject
      if (Init[t]==1){
            si<- Subject[t];
            for (v in 1:4)
              {
                Q[v] <- Q0[v];
              }// end inital values loop
       }  
       
       beta <- (2-Phase[t])*b1_ind[si]+(Phase[t]-1)*b2_ind[si];
       pchoice<-(1-epsilon)*softmax(beta*Q)+epsilon/4;
       Choice[t]~categorical(pchoice);

      // reinforcement
 
       alpha <- (2-Phase[t])*a1_ind[si]+(Phase[t]-1)*a2_ind[si];
      Q[(Choice[t])] <- Q[(Choice[t])] + alpha*(Reward[t]-Q[(Choice[t])]);

      
   }// end subject loop
}// end of model loop
