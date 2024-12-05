% 将TV-SBL推广到PC-SBL中，希望在引入元素之间的耦合性的同时，加强图像边缘的保留
function [Result]=PC_TV_SBL(y,A,beta_tilde,epi,PRINT)
% Input:
% y: measurement data
% A: measurement matrix
% lambda: noise variance
% beta_tilde: weight of TV regulaization prior, usually less than .5
% epi: little value term
%=======================Initialization=====================================
beta=0.5; % coupled coefficients
gamma=1./(abs(A'*y).^2); % prior precision
gamma0=1/var(y); % the precision of noise
Iter_max=2e2;
convergence=false;
iter=0;
[M,N]=size(A);
Ind_cor=2:2:N;
Ind_odd=1:2:N;
conver_value=1e-4;
iter_in=1;
dg_set=[];
FuncM=[];
Large=1e8;
c=1e-4;
d=1e-4;
gamma_rela=Mycouple(gamma,beta);
%========================Algorithm=========================================
while ~(iter>Iter_max||convergence)
    iter=iter+1;   
    Sigma=inv(A'*A*gamma0+diag(gamma_rela)); % covariance matrix
    mu=Sigma*A'*y*gamma0; % 1-order moment
    sigma_diag=real(diag(Sigma));
    omega=mu.*conj(mu)+sigma_diag; 
    omega_rela=Mycouple(omega,beta);
    
    gamma_left = [gamma(1);gamma(1:end-1)];
    gamma_right=[gamma(2:end);gamma(end)];
    gamma_rela = gamma+beta*(gamma_left+gamma_right);    
    gamma_rela_left=[gamma_rela(1);gamma_rela(1:end-1)];
    gamma_rela_right=[gamma_rela(2:end);gamma_rela(end)];
    Func_value=sum(omega.*gamma_rela-log(gamma_rela)+beta_tilde*norm(gradient(gamma),1));
    FuncM=[FuncM Func_value];    
    gamma_old=gamma;
    gamma_rela_old=gamma;
    Ind_cor=1:3:N;
    %================Even Indices========================================= 
    for ii=1:iter_in
    var_rela=1./gamma_rela;
    var_rela_left=[0;var_rela(1:end-1)];
    var_rela_right=[var_rela(2:end);0];
    var_sq=1./(gamma_rela-gamma);
    var_sq_left=1./(gamma_rela_left-beta*gamma);
    var_sq_right=1./(gamma_rela_right-beta*gamma);    
    omega_bar=omega_rela-var_rela-beta*(var_rela_left+var_rela_right);
    theta=var_rela.*var_sq+beta^2*(var_rela_left.*var_sq_left+var_rela_right.*var_sq_right);       
    alpha_even_hat=gamma(Ind_cor)-(omega_bar(Ind_cor)+2*beta_tilde)./theta(Ind_cor);
    alpha_even_check=gamma(Ind_cor)-omega_bar(Ind_cor)./theta(Ind_cor);
    alpha_even_bar=gamma(Ind_cor)-(omega_bar(Ind_cor)-2*beta_tilde)./theta(Ind_cor);   
    gamma_even_old=gamma;
    gamma_left=[gamma(1);gamma(1:N-1)];
    gamma_right=[gamma(2:N);gamma(N)]; 
    gamma_min=min(gamma_left,gamma_right);
    gamma_max=max(gamma_left,gamma_right);     
    
    gamma_even_hat=max(gamma_max(Ind_cor),alpha_even_hat);
    gamma_even_check=min(gamma_max(Ind_cor),max(gamma_min(Ind_cor),alpha_even_check));
    if any(alpha_even_bar<0)
        fprintf('  Gamma occurs negative element\n');
        keyboard;
    end
    gamma_even_bar=max(epi,min(gamma_min(Ind_cor),alpha_even_bar));
    Func_even_hat=Func_gamma(omega_rela,gamma_even_old,gamma_even_hat,gamma_min,gamma_max,Ind_cor,beta_tilde,beta);
    Func_even_check=Func_gamma(omega_rela,gamma_even_old,gamma_even_check,gamma_min,gamma_max,Ind_cor,beta_tilde,beta);
    Func_even_bar=Func_gamma(omega_rela,gamma_even_old,gamma_even_bar,gamma_min,gamma_max,Ind_cor,beta_tilde,beta);
    Func=[Func_even_hat,Func_even_check,Func_even_bar];
    [~,ind]=min(Func,[],2);
    GammaMtr_even=[gamma_even_hat,gamma_even_check,gamma_even_bar];
    for nn=1:length(Ind_cor)
        gamma(Ind_cor(nn))=GammaMtr_even(nn,ind(nn));
    end
    
%     gamma_left = [gamma(1);gamma(1:end-1)];
%     gamma_right=[gamma(2:end);gamma(end)];
%     gamma_rela = gamma+beta*(gamma_left+gamma_right);    
%     gamma_rela_left=[gamma_rela(1);gamma_rela(1:end-1)];
%     gamma_rela_right=[gamma_rela(2:end);gamma_rela(end)];
%     Func_value=sum(omega_rela.*gamma-log(gamma_rela)-log(gamma_rela_left)...
%         -log(gamma_rela_right)+beta_tilde*norm(gradient(gamma),1));
%     FuncM=[FuncM Func_value];
    end
    %================Odd Indices=========================================
    for ii=1:iter_in
%     gamma_left=[gamma(1);gamma(1:N-1)];
%     gamma_right=[gamma(2:N);gamma(N)]; 
%     gamma_min=min(gamma_left,gamma_right);
%     gamma_max=max(gamma_left,gamma_right);  
%     gamma_rela_odd = gamma+beta*(gamma_left+gamma_right);  
%     
%     var_rela=1./gamma_rela_odd;
%     var_rela_left=[0;var_rela(1:end-1)];
%     var_rela_right=[var_rela(2:end);0];
%     var_sq=1./(gamma_rela-gamma);
%     var_sq_left=1./(gamma_rela_left-beta*gamma);
%     var_sq_right=1./(gamma_rela_right-beta*gamma);    
%     omega_bar=omega_rela-var_rela-beta*(var_rela_left+var_rela_right);
%     theta=var_rela.*var_sq+beta^2*(var_rela_left.*var_sq_left+var_rela_right.*var_sq_right);  
    alpha_odd_hat=gamma(Ind_odd)-(omega_bar(Ind_odd)+2*beta_tilde)./theta(Ind_odd);
    alpha_odd_check=gamma(Ind_odd)-omega_bar(Ind_odd)./theta(Ind_odd);
    alpha_odd_bar=gamma(Ind_odd)-(omega_bar(Ind_odd)-2*beta_tilde)./theta(Ind_odd);    
    
    gamma_odd_hat=max(gamma_max(Ind_odd),alpha_odd_hat);
    gamma_odd_check=min(gamma_max(Ind_odd),max(gamma_min(Ind_odd),alpha_odd_check));
    gamma_odd_bar=max(0,min(gamma_min(Ind_odd),alpha_odd_bar));
    gamma_old_odd=gamma;
    Func_odd_hat=Func_gamma(omega_rela,gamma_old_odd,gamma_odd_hat,gamma_min,gamma_max,Ind_odd,beta_tilde,beta);
    Func_odd_check=Func_gamma(omega_rela,gamma_old_odd,gamma_odd_check,gamma_min,gamma_max,Ind_odd,beta_tilde,beta);
    Func_odd_bar=Func_gamma(omega_rela,gamma_old_odd,gamma_odd_bar,gamma_min,gamma_max,Ind_odd,beta_tilde,beta);
    Func=[Func_odd_hat,Func_odd_check,Func_odd_bar];
    [~,ind]=min(Func,[],2);
    GammaMtr_odd=[gamma_odd_hat,gamma_odd_check,gamma_odd_bar];
    for nn=1:length(Ind_odd)
        gamma(Ind_odd(nn))=GammaMtr_odd(nn,ind(nn));
    end  
         
    end
%     %=================================Update gamma0========================
    rho=1-sigma_diag.*gamma_rela_old;
    gamma0 =  (M+2*c)/(norm(y-A*mu)^2+sum(rho)/gamma0+2*d);
    
    gamma_left_old = [gamma_old(1);gamma_old(1:end-1)];
    gamma_right_old=[gamma_old(2:end);gamma_old(end)];
    gamma_rela = gamma+beta*(gamma_left_old+gamma_right_old);    
    gamma_rela_left=gamma_left_old+beta*gamma+beta*[gamma_old(1);gamma_left_old(1:end-1)];
    gamma_rela_right=gamma_right_old+beta*gamma+beta*[gamma_right_old(2:end);gamma_old(end)];
    Func_value=sum(omega_rela.*gamma-log(gamma_rela)-log(gamma_rela_left)...
        -log(gamma_rela_right)+beta_tilde*norm(gradient(gamma),1));
    FuncM=[FuncM Func_value];       
    %=================================Judge convergence====================
    dg=norm(gamma_old-gamma)/norm(gamma);
    dg_set=[dg_set dg];
    if dg<conver_value
        convergence=true;
        if PRINT
            fprintf("    PC+TV_SBL reach convergence at %d iteration\n", iter);
        end
    end   
end

 %==================================Output=================================
 Result.x=mu;
 Result.conver=dg_set;
 Result.gamma=gamma;
 Result.var=1/gamma0;
end

function Func_hat=Func_gamma(omega_rela,gamma_old,gamma,gamma_min,gamma_max,Ind,beta_tilde,beta)
%     gamma_old(Ind)=gamma;
    gamma_left = [gamma_old(1);gamma_old(1:end-1)];
    gamma_right=[gamma_old(2:end);gamma_old(end)];
%     gamma_rela = gamma_old+beta*(gamma_left+gamma_right); 
    gamma_dleft=[gamma_left(1);gamma_left(1:end-1)];
    gamma_dright=[gamma_right(2:end);gamma_right(end)];   
    gamma_rela=gamma+beta*(gamma_left(Ind)+gamma_right(Ind)); 
    gamma_rela_left=gamma_left(Ind)+beta*(gamma+gamma_dleft(Ind)); 
    gamma_rela_right=gamma_right(Ind)+beta*(gamma+gamma_dright(Ind)); 
    Func_hat=omega_rela(Ind).*gamma-log(gamma_rela)-log(gamma_rela_left)...
        -log(gamma_rela_right)+beta_tilde*(abs(gamma-gamma_min(Ind))+abs(gamma-gamma_max(Ind)));
end
