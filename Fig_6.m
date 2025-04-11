%% Code to reproduce Fig. 6 of Bos, Miehl et al., 2025

% Note: This code generates Fig. 6A-C. To generate Fig. 6D-F, you need to
% change these parameters in the code: wSP, rE0, rP0, rS0, and I_mod_S. You
% can do that by simply uncommenting the respective lines.

wEE=0.8; wEP=0.5; wPE=1; wPP=0.6; wSS=0; wES=0; wPS=0.8; wSE=0; 
wSP=0; % Chose these paremeters to get case Fig. 6(A-C), i.e. without PV-to-SST feedback
%wSP=0.2; % Chose these paremeters to get case Fig. 6(D-F), i.e. with PV-to-SST feedback
tauE=10; tauP=10; tauS=10; 
mult_f=1/4;
power=2; 

rE0=3.5; rP0=6; rS0=3; % Chose these paremeters to get case Fig. 6(A-C), i.e. without PV-to-SST feedback
%rE0=5; rP0=4; rS0=3; % Chose these paremeters to get case Fig. 6(D-F), i.e. with PV-to-SST feedback
rE=rE0; rP=rP0; rS=rS0;
rE_stim=rE0; rP_stim=rP0; rS_stim=rS0;

xE0=(rE/mult_f)^(power^(-1))-(wEE*rE-wEP*rP-wES*rS); 
xP0=(rP/mult_f)^(power^(-1))-(wPE*rE-wPP*rP-wPS*rS); 
xS0=(rS/mult_f)^(power^(-1))-(wSE*rE-wSP*rP-wSS*rS);
 
xE=xE0; xP=xP0; xS=xS0;
xE_stim=xE0; xP_stim=xP0; xS_stim=xS0;

dt=0.1; end_sim=1000000;
counter=0; counter2=0;

time_SOM_mod_1=end_sim/2;

I_stim_E=0.3; I_stim_P=0.3; 
I_mod_S=-0.3; % choose I_mod_S=0.3 for Fig. 6A-C and I_mod_S=-0.3 for Fig. 6D-F

rE_save(1,1)=rE0;
rP_save(1,1)=rP0;
rS_save(1,1)=rS0;

% noise parameters
tau_chi_E=50;
tau_chi_P=50;
sigma_E=1;
sigma_P=1;

std_noise_E=6;
std_noise_P=6;
chi_E=xE0; chi_P=xP0; 
chi_E_stim=xE0; chi_P_stim=xP0;

max_step=200;
max_step2=10;
bin_step=0.1;
bin_step_delta=0.01;
bin_step_beta=0.05;
bin_rates=[0:bin_step:max_step];
bin_delta_rates=[-max_step2:bin_step_delta:max_step2];
bin_beta=[0:bin_step_beta:20];

beta_bin_E=zeros(length(bin_beta),2);
count_bin_rates_E=zeros(length(bin_rates),4);
count_bin_rates_P=zeros(length(bin_rates),4);
count_bin_rates_S=zeros(length(bin_rates),4);
count_bin_rates_2D=zeros(length(bin_rates),length(bin_rates),4);
count_bin_delta_E=zeros(length(bin_delta_rates),2);

indx_change=1; indx_change_stim=3;

%%
for tt=dt:dt:end_sim
    counter=counter+1;
    counter2=counter2+1;

    if round(tt/dt)*dt==time_SOM_mod_1
        xS=xS0+I_mod_S;
        xS_stim=xS0+I_mod_S;
        
        indx_change=2;
        indx_change_stim=4;

        counter2=1;
    end

    Erand=randn; Irand=randn;

    chi_E=chi_E+(-chi_E+(xE0+std_noise_E*Erand))/tau_chi_E*dt;
    chi_P=chi_P+(-chi_P+(xP0+std_noise_P*Irand))/tau_chi_P*dt;
    xE=max(sigma_E*chi_E,0);
    xP=max(sigma_P*chi_P,0);
   
    chi_E_stim=chi_E_stim+(-chi_E_stim+(xE0+I_stim_E+std_noise_E*Erand))/tau_chi_E*dt;
    chi_P_stim=chi_P_stim+(-chi_P_stim+(xP0+I_stim_P+std_noise_P*Irand))/tau_chi_P*dt;
    xE_stim=max(sigma_E*chi_E_stim,0);
    xP_stim=max(sigma_P*chi_P_stim,0);

    
    rE_dum=rE; rP_dum=rP; rS_dum=rS;
    rE=rE_dum+dt*(-rE_dum+mult_f*(wEE*rE_dum-wEP*rP_dum-wES*rS_dum+xE)^power)/tauE;
    rP=rP_dum+dt*(-rP_dum+mult_f*(wPE*rE_dum-wPP*rP_dum-wPS*rS_dum+xP)^power)/tauP;
    rS=rS_dum+dt*(-rS_dum+mult_f*(wSE*rE_dum-wSP*rP_dum-wSS*rS_dum+xS)^power)/tauS;
    rE_dum_stim=rE_stim; rP_dum_stim=rP_stim; rS_dum_stim=rS_stim;
    rE_stim=rE_dum_stim+dt*(-rE_dum_stim+mult_f*(wEE*rE_dum_stim-wEP*rP_dum_stim-wES*rS_dum_stim+xE_stim)^power)/tauE;
    rP_stim=rP_dum_stim+dt*(-rP_dum_stim+mult_f*(wPE*rE_dum_stim-wPP*rP_dum_stim-wPS*rS_dum_stim+xP_stim)^power)/tauP;
    rS_stim=rS_dum_stim+dt*(-rS_dum_stim+mult_f*(wSE*rE_dum_stim-wSP*rP_dum_stim-wSS*rS_dum_stim+xS_stim)^power)/tauS;

    if tt<=3000
        rE_save(counter,1)=rE;
        rE_stim_save(counter,1)=rE_stim;
        rP_save(counter,1)=rP;
        rS_save(counter,1)=rS;
    end


    count_bin_delta_E(max_step2/bin_step_delta+1+round((rE_stim-rE)/bin_step_delta),indx_change)=count_bin_delta_E(max_step2/bin_step_delta+1+round((rE_stim-rE)/bin_step_delta),indx_change)+1;

    beta=power*mult_f.*sqrt(rE./mult_f);
    beta_bin_E(1+round(beta/bin_step_beta),indx_change)=beta_bin_E(1+round(beta/bin_step_beta),indx_change)+1;

    count_bin_rates_E(1+round(rE/bin_step),indx_change)=count_bin_rates_E(1+round(rE/bin_step),indx_change)+1;
    count_bin_rates_P(1+round(rP/bin_step),indx_change)=count_bin_rates_P(1+round(rP/bin_step),indx_change)+1;
    count_bin_rates_S(1+round(rS/bin_step),indx_change)=count_bin_rates_S(1+round(rS/bin_step),indx_change)+1;
 
    count_bin_rates_2D(1+round(rE/bin_step),1+round(rP/bin_step),indx_change)=count_bin_rates_2D(1+round(rE/bin_step),1+round(rP/bin_step),indx_change)+1;
    
    count_bin_rates_E(1+round(rE_stim/bin_step),indx_change_stim)=count_bin_rates_E(1+round(rE_stim/bin_step),indx_change_stim)+1;
    count_bin_rates_P(1+round(rP_stim/bin_step),indx_change_stim)=count_bin_rates_P(1+round(rP_stim/bin_step),indx_change_stim)+1;
    count_bin_rates_S(1+round(rS_stim/bin_step),indx_change_stim)=count_bin_rates_S(1+round(rS_stim/bin_step),indx_change_stim)+1;
 
    count_bin_rates_2D(1+round(rE_stim/bin_step),1+round(rP_stim/bin_step),indx_change_stim)=count_bin_rates_2D(1+round(rE_stim/bin_step),1+round(rP_stim/bin_step),indx_change_stim)+1;

end

%%

bound_no_feedback=100;
bound_feedback=120;

mean_noSOM_delta_E=sum(count_bin_delta_E(:,1).*bin_delta_rates')./sum(count_bin_delta_E(:,1));
mean_SOMp_delta_E=sum(count_bin_delta_E(:,2).*bin_delta_rates')./sum(count_bin_delta_E(:,2));

mean_noSOM_nostim_E=sum(count_bin_rates_E(:,1).*bin_rates')./sum(count_bin_rates_E(:,1));
mean_SOMp_nostim_E=sum(count_bin_rates_E(:,2).*bin_rates')./sum(count_bin_rates_E(:,2));
mean_noSOM_stim_E=sum(count_bin_rates_E(:,3).*bin_rates')./sum(count_bin_rates_E(:,3));
mean_SOMp_stim_E=sum(count_bin_rates_E(:,4).*bin_rates')./sum(count_bin_rates_E(:,4));

mean_noSOM_nostim_P=sum(count_bin_rates_P(:,1).*bin_rates')./sum(count_bin_rates_P(:,1));
mean_SOMp_nostim_P=sum(count_bin_rates_P(:,2).*bin_rates')./sum(count_bin_rates_P(:,2));


var_noSOM_nostim_E=sum(count_bin_rates_E(:,1).*(bin_rates'-mean_noSOM_nostim_E).^2)./sum(count_bin_rates_E(:,1));
var_SOMp_nostim_E=sum(count_bin_rates_E(:,2).*(bin_rates'-mean_SOMp_nostim_E).^2)./sum(count_bin_rates_E(:,2));
var_noSOM_stim_E=sum(count_bin_rates_E(:,3).*(bin_rates'-mean_noSOM_stim_E).^2)./sum(count_bin_rates_E(:,3));
var_SOMp_stim_E=sum(count_bin_rates_E(:,4).*(bin_rates'-mean_SOMp_stim_E).^2)./sum(count_bin_rates_E(:,4));

var_noSOM_nostim_P=sum(count_bin_rates_P(:,1).*(bin_rates'-mean_noSOM_nostim_P).^2)./sum(count_bin_rates_P(:,1));
var_SOMp_nostim_P=sum(count_bin_rates_P(:,2).*(bin_rates'-mean_SOMp_nostim_P).^2)./sum(count_bin_rates_P(:,2));

max_colorbar=max(max(max(count_bin_rates_2D(:,:,1)./sum(sum(count_bin_rates_2D(:,:,1))))),max(max(count_bin_rates_2D(:,:,2)./sum(sum(count_bin_rates_2D(:,:,2))))));

% make figures
figure
subplot(3,3,[2,3])
hold on 
plot([1:30000].*dt./1000,rE_save,'r')
plot([1:30000].*dt./1000,rP_save,'b')
plot([1:30000].*dt./1000,rS_save,'g')
hold off
xlabel('Time (s)')
ylabel('FR (Hz)')
ylim([0 max([rE_save;rP_save;rS_save])])

subplot(3,3,4)
imagesc(flipud(count_bin_rates_2D(:,:,1)./sum(sum(count_bin_rates_2D(:,:,1)))))
colorbar
colormap(flipud(gray))
xlim([0 bound_feedback])
ylim([2000-bound_feedback 2000])
title('no stim, no SOM mod')
axis('square')
ylabel('r_E (Hz)')
xlabel('r_P (Hz)')
xticklabels = linspace(0,10,3);
yticklabels = linspace(10,0,3);
xticks = linspace(0, 10./bin_step, numel(xticklabels));
yticks = linspace((max_step-10)./bin_step, max_step./bin_step, numel(yticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
clim([0 max_colorbar])

subplot(3,3,5)
imagesc(flipud(count_bin_rates_2D(:,:,2)./sum(sum(count_bin_rates_2D(:,:,2)))))
colorbar
xlim([0 bound_feedback])
ylim([2000-bound_feedback 2000])
title('no stim, SOM mod')
axis('square')
ylabel('r_E (Hz)')
xlabel('r_P (Hz)')
xticklabels = linspace(0,10,3);
yticklabels = linspace(10,0,3);
xticks = linspace(0, 10./bin_step, numel(xticklabels));
yticks = linspace((max_step-10)./bin_step, max_step./bin_step, numel(yticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
clim([0 max_colorbar])

subplot(3,3,6)
imagesc(flipud(count_bin_rates_2D(:,:,3)./sum(sum(count_bin_rates_2D(:,:,3)))))
colorbar
xlim([0 bound_feedback])
ylim([2000-bound_feedback 2000])
title('stim, no SOM mod')
axis('square')
ylabel('r_E (Hz)')
xlabel('r_P (Hz)')
xticklabels = linspace(0,10,3);
yticklabels = linspace(10,0,3);
xticks = linspace(0, 10./bin_step, numel(xticklabels));
yticks = linspace((max_step-10)./bin_step, max_step./bin_step, numel(yticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
clim([0 max_colorbar])

subplot(3,3,1)
imagesc(flipud(count_bin_rates_2D(:,:,4)./sum(sum(count_bin_rates_2D(:,:,4)))))
colorbar
xlim([0 bound_feedback])
ylim([2000-bound_feedback 2000])
title('stim, SOM mod')
axis('square')
ylabel('r_E (Hz)')
xlabel('r_P (Hz)')
xticklabels = linspace(0,10,3);
yticklabels = linspace(10,0,3);
xticks = linspace(0, 10./bin_step, numel(xticklabels));
yticks = linspace((max_step-10)./bin_step, max_step./bin_step, numel(yticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
clim([0 max_colorbar])

subplot(3,3,7)
hold on
plot(bin_rates,count_bin_rates_E(:,1)./sum(count_bin_rates_E(:,1)),'r--')
plot(bin_rates,count_bin_rates_E(:,2)./sum(count_bin_rates_E(:,2)),'r')
xline(mean_noSOM_nostim_E,'r--')
xline(mean_SOMp_nostim_E,'r')
hold off
xlim([0 bound_feedback/10])
legend(['no SOM mod, Var: ' num2str(round(var_noSOM_nostim_E,1))] , ['SOM mod (+), Var: ' num2str(round(var_SOMp_nostim_E,1))])
xlabel('E FR (Hz)')
ylabel('P(FR)')
title('no stim E')

subplot(3,3,8)
hold on
plot(bin_beta,beta_bin_E(:,1)./sum(beta_bin_E(:,1)),'k--')
plot(bin_beta,beta_bin_E(:,2)./sum(beta_bin_E(:,2)),'k')
hold off
xlim([0 5])
legend(['no SOM mod'] , ['SOM mod (+)'])
xlabel('b')
ylabel('P(b)')
title('no stim')

subplot(3,3,9)
hold on
plot(bin_delta_rates,count_bin_delta_E(:,1)./sum(count_bin_delta_E(:,1)),'k--')
plot(bin_delta_rates,count_bin_delta_E(:,2)./sum(count_bin_delta_E(:,2)),'k')
xline(mean_noSOM_delta_E,'k--')
xline(mean_SOMp_delta_E,'k')
hold off
xlim([0 1])
legend(['no SOM mod'] , ['SOM mod (+)'])
xlabel('Gain (g)')
ylabel('P(g)')
title('stim E')
