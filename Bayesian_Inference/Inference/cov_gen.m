%Function used to randomly generate covariance matrix for each genomic region.
%Mainly used to generate parameters for synthetic data or to initialize
%parameters for inference. Function takes a variable number of arguments
%and depending the scenario, it returns covariance matrix and unconstrained 
%standard deviation/correlation sample(s) stacked in tensors for compactness.
%Warning: This function accepts standard deviations and correlation as
%unconstrained parameters.

%cov_gen(): returns a single sample for each of the output

%cov_gen(N_samp): returns N_samp samples for each of the output paramters

%cov_gen(N_samp,sig_lt): returns a N_samp samples for each of the output
%parmaters assuming known standard deviations (Number of samples in sig_lt
%and N_samp need to be the same)

%cov_gen(N_samp,sig_lt,rho_lt): returns a N_samp samples for each of the output
%parmaters assuming known standard deviations and correlations (Number of
%samples in sig_lt/rho_lt and N_samp need to be the same)
function [CV,sig_lt,rho_lt] = cov_gen(varargin)
dim=2;
if nargin==0
   N_samp=1;
   sig_lt=randn(N_samp,dim);
   rho_lt=randn(N_samp,1);
elseif nargin==1
   N_samp=varargin{1};
   sig_lt=randn(N_samp,dim);
   rho_lt=randn(N_samp,1);    
elseif nargin==2
   N_samp=varargin{1};
   sig_lt=varargin{2};
   rho_lt=randn(N_samp,1);  
   if size(sig_lt,1)~=size(rho_lt,1)
     errorMessage = sprintf('Please check check size input paramters');
     uiwait(warndlg(errorMessage));
     return;  
   end
   elseif nargin==3
   N_samp=varargin{1};
   sig_lt=varargin{2};
   rho_lt=varargin{3};   
   if size(sig_lt,1)~=size(rho_lt,1)
     errorMessage = sprintf('Please check check size input paramters');
     uiwait(warndlg(errorMessage));
     return;  
   end
else
    errorMessage = sprintf('Please input correct number of arguments');
    uiwait(warndlg(errorMessage));
    return;   
end
rho=tanh(rho_lt/2);
sig=exp(sig_lt);

CV=Build_Mtrx(sig,rho);
    
    
end

