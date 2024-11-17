%% ------------------------------------
% This code generates Fig 3 from the paper "6G Positioning and Sensing 
% Through the Lens of Sustainability, Inclusiveness, and Trustworthiness"
% by H. Wymeersch et al, IEEE WCM, 2025. 
% (c) H. Wymeersch, November 2024, all rights reserved. 
% ------------------------------------

close all
clear all
clc
warning off 
rng(0);
% ------------------------------
% 1. Set system parameters:
% ------------------------------
par.Delta_f=40e-6;          % subcarrier spacing in Hz
par.Ns=2000;                 % number of subcarriers
par.W=par.Ns*par.Delta_f;   % system bandwidth in GHz
par.c=0.3;                  % speed of light in m/ns
par.N0_dBmHz=-174;          % noise PSD in dbm/Hz
par.NF_dB=0;                % receiver noise figure dB
par.N0=10^(0.1*(par.N0_dBmHz+par.NF_dB))*1e9; % noise PSD  [W/GHz] (290 Kelvin * Boltzmann constant in W/Hz)
geom.dim=2;                 % this is how many dimensions we estimate: 1 = X-coordinate, 2 = X&Y coordinate, 3 = not supported. I recommend setting to 1 to reduce simulation time of the estimators
par.TIR=0.1;

% UE locations where we compute the rate and PEB
Nx=100;
Ny=99;
x = linspace(0,1000,Nx);
y = linspace(0,1000,Ny);
[X,Y] = meshgrid(x,y);
geom.allUE=[X(:) Y(:) ones(Nx*Ny,1)*2]'; % all the UE locations we test

% ------------------------------
% 2. Choose variable parameter
% ------------------------------
free_parameter='base stations';
%free_parameter='power';
%free_parameter='bandwidth';
%free_parameter='frequency';



% ------------------------------
% 3. Default values for variable parameter
% ------------------------------
f_c_vec=5.9;
B_vec=par.W;
M_BS_vec=6;
P_vec=20;

% generate a grid of all the BS locations
N_BS_x=8;
x = linspace(0,1000,N_BS_x);
y = linspace(0,1000,N_BS_x);
[X,Y] = meshgrid(x,y);
geom.allBS=[X(:) Y(:) ones(N_BS_x^2,1)*10]';
figure(100)
plot3(geom.allBS(1,:),geom.allBS(2,:),geom.allBS(3,:),'rs')
geom.BS=geom.allBS(:,1:M_BS_vec);

col_num = size(geom.allBS, 2);
random_cols= randperm(col_num );
geom.allBS=geom.allBS(:,random_cols);
geom.BS=geom.allBS(:,1:M_BS_vec);
geom.M_BS=M_BS_vec;


% ------------------------------
% 4. Adapt range of variable parameters
% ------------------------------

switch free_parameter
    case {'power'}
        N_val=20;       
        P_vec=linspace(-10,15,N_val);         
        x_vector=P_vec;
    case {'bandwidth'}
        N_val=15;
        B_vec=round(logspace(-4,0,N_val)/par.Delta_f)*par.Delta_f;        
        x_vector=B_vec;
    case {'frequency'}
        N_val=10;
        f_c_vec=logspace(0.2,2.4,N_val);                   
        N_val=length(f_c_vec);
        x_vector=f_c_vec;
    case {'base stations'}              
        M_BS_vec=1:N_BS_x^2;     
        N_val=length(M_BS_vec);        
        x_vector=M_BS_vec;
    otherwise 
        error('not supported');
end


% some vectors to collect statistics
N_val=length(x_vector);
PEB_delay=zeros(N_val,Nx*Ny);
PEB_delay_robust=zeros(N_val,Nx*Ny);
Rate=zeros(N_val,Nx*Ny);
PEB_delay_cov=zeros(N_val,1);
Rate_cov=zeros(N_val,1);
PEB_delay_avg=zeros(N_val,1);
Rate_avg=zeros(N_val,1);
power_comm=zeros(N_val,1);
energy_pos=zeros(N_val,1);

Precision_avg=zeros(N_val,1);
Precision_cov=zeros(N_val,1);



% ------------------------------
% 5. Start simulation, loop over free parameter
% ------------------------------
for counter=1:N_val                                 % loop over the free parameter
    counter    
    % set the remaining parametrs
    switch free_parameter
        case {'power'}
            par.P_dBm=P_vec(counter);                % transmit power in dBm        
            par.P=10^(0.1*par.P_dBm);                % transmit power in mW              
            par.W=B_vec;
            par.Ns=round(par.W/par.Delta_f);        
            par.fc=f_c_vec;                         % carrier in GHz
            par.lambda=par.c/par.fc;                % wavelength in meter
            geom.M_BS=M_BS_vec;
            geom.BS=geom.allBS(:,1:geom.M_BS);            
        case {'bandwidth'}
            par.P_dBm=P_vec;                        % transmit power in dBm        
            par.P=10^(0.1*par.P_dBm);               % transmit power in mW  
            par.W=B_vec(counter);
            par.Ns=round(par.W/par.Delta_f);      
            par.fc=f_c_vec;                         % carrier in GHz
            par.lambda=par.c/par.fc;                % wavelength in meter
            geom.M_BS=M_BS_vec;
            geom.BS=geom.allBS(:,1:geom.M_BS);            
        case {'frequency'}
            par.P_dBm=P_vec;                % transmit power in dBm        
            par.P=10^(0.1*par.P_dBm);                % transmit power in mW  
            par.W=B_vec;
            par.Ns=round(par.W/par.Delta_f);        
            par.fc=f_c_vec(counter);               % carrier in GHz
            par.lambda=par.c/par.fc;                % wavelength in meter
             geom.M_BS=M_BS_vec;
            geom.BS=geom.allBS(:,1:geom.M_BS);            
        case {'base stations'}            
            par.P_dBm=P_vec;                % transmit power in dBm        
            par.P=10^(0.1*par.P_dBm);                % transmit power in mW  
            par.W=B_vec;
            par.Ns=round(par.W/par.Delta_f);        
            par.fc=f_c_vec;               % carrier in GHz
            par.lambda=par.c/par.fc;                % wavelength in meter           
            geom.M_BS=M_BS_vec(counter);
            geom.BS=geom.allBS(:,1:geom.M_BS);            
    end
    
    M=geom.M_BS;    
    power_comm(counter)=par.P;
    energy_pos(counter)= par.P*M;
    for k=1:Nx*Ny   % for all UE positions
        geom.UE=geom.allUE(:,k);            % fix the location        
        d=vecnorm(geom.UE-geom.BS)';        % compute distances
        Es = par.P/(par.Ns*par.Delta_f);    % energy per subcarrier
        rho=par.lambda./(4*pi*d);           % path loss
        SNR=Es*rho.^2/par.N0;               % SNR per subcarrier
        Rate(counter,k)=par.Ns*par.Delta_f*log2(1+max(SNR));    % rate 
        
        % compute the FIM of the distances obtained from TOA 
        J_tau=diag(2*SNR*pi^2*(par.W)^2/(3*(par.c)^2));     % distance estimates from delays               
               
        [PEB_delay(counter,k),J_delay,C_position]=computePEB_delay(par,geom,J_tau,1:M);  % delay only 
        if (abs(imag(PEB_delay(counter,k)))>0)
            PEB_delay(counter,k)=+inf; % if the UE location cannot be estimated
            C_position=eye(2)*(+Inf);
        end  
        PEB_delay_robust(counter,k)=+inf;
        
        sigma_pos(1)=sqrt(C_position(1,1));
        sigma_pos(2)=sqrt(C_position(2,2));
        
        for n=1:2
            PL(n)=findPL(sigma_pos(n),par.TIR/2);
        end
        
                
        PL_XY=sqrt(PL(1)^2+PL(2)^2);

        if abs(imag(PL_XY))>0
            keyboard
        end        
        PEB_delay_robust(counter,k)=PL_XY;        

    end
    
    Rate_cov(counter)=prctile(Rate(counter,:),5);               % 95% coverage (higher is better)
    Rate_avg(counter)=median(Rate(counter,:));                  % median value, since mean is inf in many cases  
    
    PEB_delay_cov(counter)=prctile(PEB_delay(counter,:),95);    % 95% coverage (lower is better)
    PEB_delay_avg(counter)=median(PEB_delay(counter,:));        % median value, since mean is inf in many cases
    PEB_delay_robust_cov(counter)=prctile(PEB_delay_robust(counter,:),95);    % 95% coverage (lower is better)
    PEB_delay_robust_avg(counter)=median(PEB_delay_robust(counter,:));        % median value, since mean is inf in many cases
    

    Precision_cov(counter)=prctile(1./(PEB_delay(counter,:)).^2,95);    % estimate of precision in (m^-2). Not super-correct
    Precision_avg(counter)=median(1./(PEB_delay(counter,:)).^2);        % estimate of precision in (m^-2). Not super-correct
end


% ------------------------------
% 6. Plot the results
% ------------------------------
figure
yyaxis left
ha=semilogy(x_vector,PEB_delay_cov,'b',x_vector,PEB_delay_avg,'b--',x_vector,PEB_delay_robust_avg, 'b:','MarkerIndices',1:10:length(x_vector),'LineWidth',2);                                       
yyaxis left
ylim([0.01 1000]);
ylabel('RMSE [m]','Interpreter','latex','FontSize',14)
yyaxis right
hb=plot(x_vector,Rate_cov,'r',x_vector,Rate_avg,'r--','LineWidth',2);
xlim([1 N_BS_x^2]);
legend('inclusive accuracy','median accuracy','trustworthy accuracy','inclusive rate','median rate','FontSize',12);                                        
%legend('95% PEB','median PEB','95% rate','median rate','Interpreter','latex','FontSize',14);    % <- this didn't work :(                                   
ylabel('rate [Gb/s]','Interpreter','latex','FontSize',14)
switch free_parameter
    case {'power'}       
        xlabel('power in dBm','Interpreter','latex','FontSize',14)   
    case {'bandwidth'}
        xlabel('bandwidth in GHz','Interpreter','latex','FontSize',14)        
    case {'frequency'}        
        xlabel('carrier in GHz','Interpreter','latex','FontSize',14)        
    case {'base stations'}        
        xlabel('number of BSs (CAPEX and OPEX)','Interpreter','latex','FontSize',14)        
end
set(gcf, 'Color', 'w');
