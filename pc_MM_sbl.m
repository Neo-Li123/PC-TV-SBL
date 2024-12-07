% 利用MM算法给出PC-SBL中M步的迭代解，从而使得有收敛性
function [Result]=pc_MM_sbl(y,A,beta,dim,PRINT)
% Input:
% y: measurement data
% A: measurement matrix
% lambda: noise variance,
% beta_tilde: weight of TV regulaization prior, usually less than .5
% epi: little value term
%=======================Initialization=====================================
b=1e-8;
a=1+1e-2;
c=1e-4;d=1e-4;
% beta=5e-1; % coupled coefficients
[M,N]=size(A);
% gamma=1./(abs(A'*y).^2); % prior precision
gamma=ones(N,1);
gamma0=1/var(y); % the precision of noise
Iter_max=2e2;
convergence=false;
iter=0;
conver_value=1e-4;
dg_set=[];
FuncM=[];
FuncM_up=[];
Large=1e8;
x_len=sqrt(N);
y_len=N/x_len;
if dim==1
    iter_in_y=1;
else
    iter_in_y=3;
end
iter_in_x=3;
%========================Algorithm=========================================
while ~(iter>Iter_max||convergence)
    iter=iter+1;
    %=======================================E Step=========================
    gamma_rela=Mycouple(gamma,beta,dim);
    Sigma=inv(A'*A*gamma0+diag(gamma_rela)); % covariance matrix
    mu=Sigma*A'*y*gamma0; % 1-order moment
    sigma_diag=real(diag(Sigma));
    omega=mu.*conj(mu)+sigma_diag;
    omega_rela=Mycouple(omega,beta,dim);
    
    gamma_rela_left=[gamma_rela(1);gamma_rela(1:end-1)];
    gamma_rela_right=[gamma_rela(2:end);gamma_rela(end)];
    %         gamma_rela_cp=gamma_rela-gamma;
    Func_value_test=omega_rela.*gamma-log(gamma_rela)-log(gamma_rela_left)-log(gamma_rela_right);
    %         Func_value=sum(-2*(a-1)*log(gamma)+2*b*gamma+omega.*gamma_rela-log(gamma_rela));
    %         FuncM=[FuncM Func_value];
    %     if iter>1
    %         if FuncM(end)>FuncM(end-1)
    %             keyboard;
    %         end
    %     end
    %=======================M Step=========================================
    gamma_old_org=gamma;
    
    gamma_rela_old=gamma_rela;
    Ind_cor_y=0:3:y_len;
    for yy=1:iter_in_y
        Ind_cor_x=0:3:x_len;
        Ind_cor_y=Ind_cor_y+1;
        Ind_cor_y(Ind_cor_y>y_len)=[];
%         disp(['  y axis=',num2str(Ind_cor_y)]);
        for nn=1:iter_in_x
            Ind_cor_x=Ind_cor_x+1;
            Ind_cor_x(Ind_cor_x>x_len)=[];
            Ind_cor=Ind_cor_x+(Ind_cor_y-1)'*x_len;
%             Ind_cor=(yy-1)*x_len+Ind_cor;
%             Ind_cor(Ind_cor>N)=[];
            
%             disp(['  x axis=',num2str(Ind_cor_x)]);
            gamma_old=gamma;
            var_rela=1./gamma_rela;
            var_sum=Mycouple(var_rela,beta,dim);
            omega_bar=omega_rela-var_sum;
            var_sq=1./(gamma_rela-gamma);
            var_sq_left=1./(gamma_rela_left-beta*gamma);
            var_sq_right=1./(gamma_rela_right-beta*gamma);
            if dim==1
                theta=var_sq./gamma_rela+beta^2*(var_sq_left./gamma_rela_left+var_sq_right./gamma_rela_right);
            elseif dim==2
                gamma_rela_up=[gamma_rela(x_len+1:end);gamma_rela(end-x_len+1:end)];
                gamma_rela_down=[gamma_rela(1:x_len);gamma_rela(1:end-x_len)];
                var_sq_up=1./(gamma_rela_up-beta*gamma);
                var_sq_down=1./(gamma_rela_down-beta*gamma);
                theta=var_sq./gamma_rela+beta^2*(var_sq_left./gamma_rela_left+var_sq_right./gamma_rela_right...
                    +var_sq_up./gamma_rela_up+var_sq_down./gamma_rela_down);
            else
                error('   The input of dimension is 1 or 2\n ');
            end
            %         alpha=gamma-omega_bar./theta;
            sig_term=omega_bar+2*b-theta.*gamma;
            judge_term=sig_term.^2+8*(a-1)*theta;
            Ind_pos=find(judge_term>0);
            %         if any(judge_term<0)
            %             keyboard;
            %         end
            %         alpha_left=(-sig_term(Ind_pos)-sqrt(judge_term(Ind_pos)))./(2*theta(Ind_pos));
            alpha=(-sig_term(Ind_pos)+sqrt(judge_term(Ind_pos)))./(2*theta(Ind_pos));
            %         if any(alpha<0)
            %             fprintf('  Gamma occurs negative element\n');
            %                     keyboard;
            %         end
            gamma(Ind_cor)=alpha(Ind_cor);
            %         gamma(alpha<0)=gamma_old(alpha<0)*0.99;
            gamma_rela = Mycouple(gamma,beta,dim);
            gamma_rela_left=[gamma_rela(1);gamma_rela(1:end-1)];
            gamma_rela_right=[gamma_rela(2:end);gamma_rela(end)];
            %============================测试上界的下降=====================
            Func_value=sum(-2*(a-1)*log(gamma)+2*b*gamma+omega.*gamma_rela-log(gamma_rela));
            FuncM=[FuncM  double(Func_value)];
            [ind,F,F_org,F_up]=Func_org(omega_rela,gamma,gamma_old,Ind_cor(2:end-1),beta,a,b);
            Err=[F(ind),F_org(ind),F_up(ind)];
            %         if iter>1
            %         if ~isempty(ind)
            %             keyboard;
            %         end
%                 if iter>1
%                     if FuncM(end)>FuncM(end-1)
%                         keyboard;
%                     end
%                 end
%             dg_in=norm(gamma_old-gamma)/norm(gamma);
%             if dg_in<conver_value
%                 break;
%             end
        end
    end
    %     %=================================Update gamma0========================
    rho=1-sigma_diag.*gamma_rela_old;
    gamma0 =  (M+2*c)/(norm(y-A*mu)^2+sum(rho)/gamma0+2*d);
    
    % prevent error of inverse
    gamma(gamma>1e8)=1e8;
    %=================================Judge convergence====================
    dg=norm(gamma_old_org-gamma)/norm(gamma);
    dg_set=[dg_set dg];
    if dg<conver_value
        convergence=true;
        if PRINT
            fprintf("    PC-MM-SBL reach convergence at %d iteration\n", iter);
        end
    end
end

%==================================Output=================================
Result.x=mu;
Result.conver=dg_set;
Result.gamma=gamma;
Result.FuncM=FuncM;
Result.var=1/gamma0;
end

function [ind,F,F_org,F_up]=Func_org(omega_rela,gamma,gamma_old,idx,beta,a,b)
gamma_l=gamma(idx-1);
gamma_ll=gamma(idx-2);
gamma_r=gamma(idx+1);
gamma_rr=gamma(idx+2);
gamma_l_old=gamma_old(idx-1);
gamma_ll_old=gamma_old(idx-2);
gamma_r_old=gamma_old(idx+1);
gamma_rr_old=gamma_old(idx+2);
gamma_c_old=gamma_old(idx)+beta*(gamma_l_old+gamma_r_old);
gamma_c_r_old=gamma_r_old+beta*(gamma_rr_old+gamma_old(idx));
gamma_c_l_old=gamma_l_old+beta*(gamma_ll_old+gamma_old(idx));
F=-2*(a-1)*log(gamma(idx))+(omega_rela(idx)+2*b).*gamma(idx)-log(gamma(idx)+beta*(gamma_l+gamma_r))-log(gamma_l+beta*(gamma_ll+gamma(idx)))...
    -log(gamma_r+beta*(gamma_rr+gamma(idx)));
F_org=-2*(a-1)*log(gamma_old(idx))+(omega_rela(idx)+2*b).*gamma_old(idx)-log(gamma_c_old)-log(gamma_c_l_old)-log(gamma_c_r_old);
% omega_bar=omega_rela(idx)-1./(gamma_c_old)-beta./(gamma_c_r_old)-beta./(gamma_c_l_old);
theta=1./gamma_c_old./(gamma_c_old-gamma_old(idx))+beta^2./gamma_c_r_old./(gamma_c_r_old-beta*gamma_old(idx))...
    +beta^2./gamma_c_l_old./(gamma_c_l_old-beta*gamma_old(idx));
F_up=-2*(a-1)*log(gamma(idx))+(omega_rela(idx)+2*b).*gamma(idx)+0.5*theta.*(gamma(idx)-gamma_old(idx)).^2-...
    log(gamma_c_old)-log(gamma_c_l_old)-log(gamma_c_r_old)+(gamma_old(idx)-gamma(idx)).*(1./gamma_c_old+beta./gamma_c_l_old+beta./gamma_c_r_old);
ind=find(F>F_org);
end
