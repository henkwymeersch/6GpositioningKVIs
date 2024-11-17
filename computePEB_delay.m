function [PEB_delay,J_delay,C_position]=computePEB_delay(par,geom,J_tau,x_index)
% PEB and FIM for delay only positioning
    U=(geom.UE-geom.BS(:,x_index))./vecnorm(geom.UE-geom.BS(:,x_index));      % unit norm vectors used in FIM
    U=U(1:geom.dim,:);
     J_delay=[U*(J_tau)*U' par.c*U*diag(J_tau) ; 
                     (par.c*U*diag(J_tau))' trace(J_tau)*(par.c)^2];   % FIM expression from Henk's GC 2023 paper    

     if (rank(J_delay)<geom.dim+1)          % location and clock bias
         PEB_delay=+Inf;
         C_position=+Inf*eye(geom.dim);
         return;
     end
     tmp=inv(J_delay);
     PEB_delay=sqrt(trace(tmp(1:geom.dim,1:geom.dim))); 
     C_position=tmp(1:geom.dim,1:geom.dim);