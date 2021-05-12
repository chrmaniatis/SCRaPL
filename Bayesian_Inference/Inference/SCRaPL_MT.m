clc
clear all %clear all variables in the workspace

parpool(32)%Start a pool with 32 workers. If parallel computing toolbox please comment that line.

[parentdir,~,~]=fileparts(pwd);%Get parent directory.
DataFolder=[parentdir,'/Data'];%Folder with data
SaveFolder=[parentdir,'/Results'];%Saving folder

File=strcat("MT_enh_MT_25000_0.1.mat");%Data file name
load(fullfile(DataFolder,File))%loading data

File_meta="MT_enh_MT_25000_0.1_Meta.mat";%Name of with wchich meta-data will be saved

N_genes=size(unique(y_true(:,4)),1);%Number of genomic regions

%Initializations

%Use data to initialize latent covariance,mean and state.
x_lat=[norminv(0.0000001+0.99999*y_true(:,1)./y_true(:,2)),log(1+y_true(:,3)./norm_fact_fact_all.norm_fact(y_true(:,5))),y_true(:,4)];%Initialize posterior latent states
m_gen=[accumarray(y_true(:,4),x_lat(:,1),[],@mean),accumarray(y_true(:,4),x_lat(:,2),[],@mean)];%Initialize posterior latent means
sig_gen=[accumarray(y_true(:,4),x_lat(:,1),[],@std),accumarray(y_true(:,4),x_lat(:,2),[],@std)];%Initialize posterior latent standard deviations
cor_gen= 2*atanh(prs_ft_cor_MT(y_true(:,4),y_true));%Initialize posterior latent correlations

m_gene_new=zeros(N_genes,2);%posterior latent mean initialization
x_lat=[x_lat,y_true(:,end-1)];%latent state values initialization
zero_inf=-1+2*randn(N_genes,1); %zero inflation initilization (in logit-scale) 

alpha_inf=2;beta_inf=8;%zero inflation prior parameters initalization
m_prior=2*randn(N_genes,2); Psi_0=3*eye(2); %latent mean prior parameters initialization

%Prior parameters in rest components initialization
a_gene=[3,3];b_gene=[1,1];%gene prior covariance parameters
alpha_gene=3;beta_gene=3;%prior correlation parameters

%Arrays for collected samples initialization
Col_smps=1050;%Number of collected samples
Burn_in=450;%Warm-up steps
sig_gen_col=zeros(N_genes,2,Col_smps);%Collected std for each gene
cor_gen_col=zeros(N_genes,Col_smps);%Collected correlation for each gene
zero_inf_col=zeros(N_genes,Col_smps);%Collected zero-inflation for each gene
m_gen_col=zeros(N_genes,2,Col_smps);%Collected means for each gene
x_mn=0*x_lat;%Estimated mean for each gene from latent states
lg_x_1=zeros(N_genes,2,Col_smps);%Log-likelihood pt1
lg_x_2=zeros(N_genes,Col_smps);%Log-likelihood pt2
step_times=zeros(3,1);%collect times to complete each inference step (useful for predicting remaining time)
ign=[];%Indices of ignored features

for jj=1:(Burn_in+Col_smps)
tic
    Sima_gene=Build_Mtrx(exp(sig_gen),tanh(cor_gen/2));%Gene Covarainces
    ssg=zeros(N_genes,3);
    m_gen1=[];
    zero_inf1=[];
    x_col=[]; 
    parfor ii=1:N_genes
         yy=y_true(y_true(:,end-1)==ii,:);
         xx=x_lat(x_lat(:,end)==ii,:);

        if ismember(ii,ign)==1
            zero_inf1=[zero_inf1;zero_inf(ii)];
            m_gen1=[m_gen1;m_prior(ii,:)];
            x_col=[x_col;xx];
            continue
        end
         
	     log_pdf_x1=@(prmx1) llik_ps_smp_x1_bino(prmx1,yy(:,1:2),m_gen(ii,:),sig_gen(ii,:),cor_gen(ii),xx(:,2));
         smp_x1 = hmcSampler(log_pdf_x1,xx(:,1));
         xx(:,1) = drawSamples(smp_x1,'NumSamples',1,'Burnin',10); %sample latent methylation/accessibility

         log_pdf_x2=@(prmx2) llik_ps_smp_x2_zero_inf_new(prmx2,yy(:,3),m_gen(ii,:),sig_gen(ii,:),cor_gen(ii),xx(:,1),zero_inf(ii),table2array(norm_fact_all(y_true(:,4)==ii,1)));
         smp_x2 = hmcSampler(log_pdf_x2,xx(:,2));
         xx(:,2) = drawSamples(smp_x2,'NumSamples',1,'Burnin',10); %sample latent expression
                  
         sm=sum(xx(:,1:2),1);
         xx1=xx(:,1:2)-m_gen(ii,:);
         tt=[sum(xx1(:,1).^2),sum(xx1(:,2).^2),sum(xx1(:,1).*xx1(:,2))];
         
         log_pdf_gene=@(gparam)lg_std_rho_samp(gparam,tt,a_gene,b_gene,alpha_gene,beta_gene,size(xx,1));
         smp_sg_g = hmcSampler(log_pdf_gene,[sig_gen(ii,:),cor_gen(ii)]);
         [smp_sg_g,tuneinfo_sg_g] = tuneSampler(smp_sg_g);%Tune HMC sampler
         ssg(ii,:) = drawSamples(smp_sg_g,'NumSamples',1,'Burnin',4); %sample latent covariance
         
         m_gene_new(ii,:)=((Psi_0)\(m_prior(ii,:))'+Sima_gene(:,:,ii)\(sm'))';
         Psi_hat=inv(Psi_0)+size(xx,1)*inv(Sima_gene(:,:,ii));
         m_gene_new(ii,:)=(Psi_hat\m_gene_new(ii,:)')';
         
         m_gen1=[m_gen1;mvnrnd(m_gene_new(ii,:),inv(Psi_hat))]; % sample latent mean

         
         x_col=[x_col;xx];
         log_pdf_zinf=@(zinfparam)llik_ps_smp_zero_inf(zinfparam,yy(:,3),xx(:,2),alpha_inf,beta_inf,table2array(norm_fact_all(y_true(:,4)==ii,1)));
         smp_zinf=hmcSampler(log_pdf_zinf,zero_inf(ii));%'UseNumericalGradient',true
         [smp_zinf,tuneinfo_zinf] = tuneSampler(smp_zinf);%Tune HMC sampler
         zero_inf1=[zero_inf1;drawSamples(smp_zinf,'NumSamples',1,'Burnin',4)]; % sample latent inflation    
    
    end
    
    m_gen=m_gen1; sig_gen=ssg(:,1:2); cor_gen=ssg(:,3); x_lat=x_col; zero_inf=zero_inf1;
    clear m_gen1 ssg x_col zero_inf1
      %If inference exceeds burn in steps then record parameters of interest
      if jj>Burn_in
       ttt=mod(jj-Burn_in-1,Col_smps)+1;
       sig_gen_col(:,:,ttt)=sig_gen;
       cor_gen_col(:,ttt)=cor_gen;
       m_gen_col(:,:,ttt)=m_gen;
       zero_inf_col(:,ttt)=zero_inf;
       x_mn=x_mn+x_lat/Col_smps;
       
       [llk1,llk2]=llk(y_true,x_lat,m_gen,sig_gen,cor_gen,zero_inf,norm_fact_all);
       lg_x_1(:,:,ttt)=llk1;
       lg_x_2(:,ttt)=llk2;
       
      end
    if mod(jj,20)==0    
        save(strjoin([SaveFolder,File_meta],'/')) %save current state (helpfull for recovery)
    end      
    step_times(mod(jj,3)+1)=toc;
    [ mean(tanh(cor_gen'/2)),jj/(Burn_in+Col_smps)] %print average latent correlation across genomic regions and percentage of completeness 
    remaining_time=(Burn_in+Col_smps-jj)*mean(step_times) %print remaining time 
end
