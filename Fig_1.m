%% Code to reproduce Fig. 1 of Bos, Miehl et al., 2024


% parameters
power=2; % power of the I/O nonlinearity
mult_f=1/4; % multiplyer of nonlinearity

dt=0.01; % simulation timestep
end_sim=350; % length of stimulation
tauE=10; tauP=10; tauS=10; % time constants for each population
I_mod_S=1; % amount of SST perturbation

%% Case 1 - Inhibition
counter=0;

wEE=0.8; wEP=0.5; wPE=1; wPP=0.6; % E-PV parameters
wSS=0; % SST self-connection
wES=0.2; % inhibitory pathway
wPS=0; % disinhibitory pathway
wSE=0;% FB from E to SST
wSP=0; % FB from PV to SST

rS=2.5; rP=5; rE=5; % initial firing rates for each population

rE_save(1,1)=rE;
rP_save(1,1)=rP;
rS_save(1,1)=rS;

% calculate the inputs for each population from the firing rate steady
% states rX=mult_f*(wXE*rE-wXP*rP-wXS*rS+xX0)^power for X=[E,P,S]
xE0=(rE_save(1,1)/mult_f)^(power^(-1))-(wEE*rE_save(1,1)-wEP*rP_save(1,1)-wES*rS_save(1,1)); 
xP0=(rP_save(1,1)/mult_f)^(power^(-1))-(wPE*rE_save(1,1)-wPP*rP_save(1,1)-wPS*rS_save(1,1)); 
xS0=(rS_save(1,1)/mult_f)^(power^(-1))-(wSE*rE_save(1,1)-wSP*rP_save(1,1)-wSS*rS_save(1,1));
xE=xE0; xP=xP0; xS=xS0;

for tt=dt:dt:end_sim
    counter=counter+1;
    
    if round(tt/dt)*dt==50
        xS=xS0+I_mod_S;
    end
    
    rE_save(counter+1,1)=rE_save(counter,1)+dt*(-rE_save(counter,1)+mult_f*(wEE*rE_save(counter,1)-wEP*rP_save(counter,1)-wES*rS_save(counter,1)+xE)^power)/tauE;
    rP_save(counter+1,1)=rP_save(counter,1)+dt*(-rP_save(counter,1)+mult_f*(wPE*rE_save(counter,1)-wPP*rP_save(counter,1)-wPS*rS_save(counter,1)+xP)^power)/tauP;
    rS_save(counter+1,1)=rS_save(counter,1)+dt*(-rS_save(counter,1)+mult_f*(wSE*rE_save(counter,1)-wSP*rP_save(counter,1)-wSS*rS_save(counter,1)+xS)^power)/tauS;

end

%% Case 1 - Disinhibition
case_dyn=2;
counter=0;

wEE=0.8; wEP=0.5; wPE=1; wPP=0.6; % E-PV parameters
wSS=0; % SST self-connection
wES=0; % inhibitory pathway
wPS=0.2; % disinhibitory pathway
wSE=0;% FB from E to SST
wSP=0; % FB from PV to SST

rS=2.5; rP=5; rE=5;

rE_save(1,case_dyn)=rE;
rP_save(1,case_dyn)=rP;
rS_save(1,case_dyn)=rS;

xE0=(rE_save(1,case_dyn)/mult_f)^(power^(-1))-(wEE*rE_save(1,case_dyn)-wEP*rP_save(1,case_dyn)-wES*rS_save(1,case_dyn)); 
xP0=(rP_save(1,case_dyn)/mult_f)^(power^(-1))-(wPE*rE_save(1,case_dyn)-wPP*rP_save(1,case_dyn)-wPS*rS_save(1,case_dyn)); 
xS0=(rS_save(1,case_dyn)/mult_f)^(power^(-1))-(wSE*rE_save(1,case_dyn)-wSP*rP_save(1,case_dyn)-wSS*rS_save(1,case_dyn));
xE=xE0; xP=xP0; xS=xS0;

for tt=dt:dt:end_sim
    counter=counter+1;
    
    if round(tt/dt)*dt==50
        xS=xS0+I_mod_S;
    end
    
    rE_save(counter+1,case_dyn)=rE_save(counter,case_dyn)+dt*(-rE_save(counter,case_dyn)+mult_f*(wEE*rE_save(counter,case_dyn)-wEP*rP_save(counter,case_dyn)-wES*rS_save(counter,case_dyn)+xE)^power)/tauE;
    rP_save(counter+1,case_dyn)=rP_save(counter,case_dyn)+dt*(-rP_save(counter,case_dyn)+mult_f*(wPE*rE_save(counter,case_dyn)-wPP*rP_save(counter,case_dyn)-wPS*rS_save(counter,case_dyn)+xP)^power)/tauP;
    rS_save(counter+1,case_dyn)=rS_save(counter,case_dyn)+dt*(-rS_save(counter,case_dyn)+mult_f*(wSE*rE_save(counter,case_dyn)-wSP*rP_save(counter,case_dyn)-wSS*rS_save(counter,case_dyn)+xS)^power)/tauS;

end

%% Case 2 - strong PV-to-PV self-connection
case_dyn=3;
counter=0;

wEE=0.8; wEP=1; wPE=1; wPP=1; % E-PV parameters
wSS=0; 
wES=0.5; % inhibitory pathway
wPS=0.6; % disinhibitory pathway
wSE=0;% FB from E to SST
wSP=0; % FB from PV to SST

rS=2.5; rP=5; rE=5;

rE_save(1,case_dyn)=rE;
rP_save(1,case_dyn)=rP;
rS_save(1,case_dyn)=rS;

xE0=(rE_save(1,case_dyn)/mult_f)^(power^(-1))-(wEE*rE_save(1,case_dyn)-wEP*rP_save(1,case_dyn)-wES*rS_save(1,case_dyn)); 
xP0=(rP_save(1,case_dyn)/mult_f)^(power^(-1))-(wPE*rE_save(1,case_dyn)-wPP*rP_save(1,case_dyn)-wPS*rS_save(1,case_dyn)); 
xS0=(rS_save(1,case_dyn)/mult_f)^(power^(-1))-(wSE*rE_save(1,case_dyn)-wSP*rP_save(1,case_dyn)-wSS*rS_save(1,case_dyn));
xE=xE0; xP=xP0; xS=xS0;

for tt=dt:dt:end_sim
    counter=counter+1;
    
    if round(tt/dt)*dt==50
        xS=xS0+I_mod_S;
    end
    
    rE_save(counter+1,case_dyn)=rE_save(counter,case_dyn)+dt*(-rE_save(counter,case_dyn)+mult_f*(wEE*rE_save(counter,case_dyn)-wEP*rP_save(counter,case_dyn)-wES*rS_save(counter,case_dyn)+xE)^power)/tauE;
    rP_save(counter+1,case_dyn)=rP_save(counter,case_dyn)+dt*(-rP_save(counter,case_dyn)+mult_f*(wPE*rE_save(counter,case_dyn)-wPP*rP_save(counter,case_dyn)-wPS*rS_save(counter,case_dyn)+xP)^power)/tauP;
    rS_save(counter+1,case_dyn)=rS_save(counter,case_dyn)+dt*(-rS_save(counter,case_dyn)+mult_f*(wSE*rE_save(counter,case_dyn)-wSP*rP_save(counter,case_dyn)-wSS*rS_save(counter,case_dyn)+xS)^power)/tauS;

end

%% Case 2 - weak PV-to-PV self-connection
case_dyn=4;
counter=0;

wEE=0.8; wEP=1; wPE=1; wPP=0.1; % E-PV parameters
wSS=0; 
wES=0.5; % inhibitory pathway
wPS=0.6; % disinhibitory pathway
wSE=0;% FB from E to SST
wSP=0; % FB from PV to SST

rS=2.5; rP=5; rE=5;

rE_save(1,case_dyn)=rE;
rP_save(1,case_dyn)=rP;
rS_save(1,case_dyn)=rS;

xE0=(rE_save(1,case_dyn)/mult_f)^(power^(-1))-(wEE*rE_save(1,case_dyn)-wEP*rP_save(1,case_dyn)-wES*rS_save(1,case_dyn)); 
xP0=(rP_save(1,case_dyn)/mult_f)^(power^(-1))-(wPE*rE_save(1,case_dyn)-wPP*rP_save(1,case_dyn)-wPS*rS_save(1,case_dyn)); 
xS0=(rS_save(1,case_dyn)/mult_f)^(power^(-1))-(wSE*rE_save(1,case_dyn)-wSP*rP_save(1,case_dyn)-wSS*rS_save(1,case_dyn));
xE=xE0; xP=xP0; xS=xS0;

for tt=dt:dt:end_sim
    counter=counter+1;
    
    if round(tt/dt)*dt==50
        xS=xS0+I_mod_S;
    end
    
    rE_save(counter+1,case_dyn)=rE_save(counter,case_dyn)+dt*(-rE_save(counter,case_dyn)+mult_f*(wEE*rE_save(counter,case_dyn)-wEP*rP_save(counter,case_dyn)-wES*rS_save(counter,case_dyn)+xE)^power)/tauE;
    rP_save(counter+1,case_dyn)=rP_save(counter,case_dyn)+dt*(-rP_save(counter,case_dyn)+mult_f*(wPE*rE_save(counter,case_dyn)-wPP*rP_save(counter,case_dyn)-wPS*rS_save(counter,case_dyn)+xP)^power)/tauP;
    rS_save(counter+1,case_dyn)=rS_save(counter,case_dyn)+dt*(-rS_save(counter,case_dyn)+mult_f*(wSE*rE_save(counter,case_dyn)-wSP*rP_save(counter,case_dyn)-wSS*rS_save(counter,case_dyn)+xS)^power)/tauS;

end

%% Case 3 - low PV firing rate
case_dyn=5;
counter=0;

wEE=0.8; wEP=1; wPE=1; wPP=0.6; % E-PV parameters
wSS=0; 
wES=0.5; % inhibitory pathway
wPS=0.6; % disinhibitory pathway
wSE=0;% FB from E to SST
wSP=0; % FB from PV to SST

rS=2.5; rP=1.8; rE=5;

rE_save(1,case_dyn)=rE;
rP_save(1,case_dyn)=rP;
rS_save(1,case_dyn)=rS;

xE0=(rE_save(1,case_dyn)/mult_f)^(power^(-1))-(wEE*rE_save(1,case_dyn)-wEP*rP_save(1,case_dyn)-wES*rS_save(1,case_dyn)); 
xP0=(rP_save(1,case_dyn)/mult_f)^(power^(-1))-(wPE*rE_save(1,case_dyn)-wPP*rP_save(1,case_dyn)-wPS*rS_save(1,case_dyn)); 
xS0=(rS_save(1,case_dyn)/mult_f)^(power^(-1))-(wSE*rE_save(1,case_dyn)-wSP*rP_save(1,case_dyn)-wSS*rS_save(1,case_dyn));
xE=xE0; xP=xP0; xS=xS0;

for tt=dt:dt:end_sim
    counter=counter+1;
    
    if round(tt/dt)*dt==50
        xS=xS0+I_mod_S;
    end
    
    rE_save(counter+1,case_dyn)=rE_save(counter,case_dyn)+dt*(-rE_save(counter,case_dyn)+mult_f*(wEE*rE_save(counter,case_dyn)-wEP*rP_save(counter,case_dyn)-wES*rS_save(counter,case_dyn)+xE)^power)/tauE;
    rP_save(counter+1,case_dyn)=rP_save(counter,case_dyn)+dt*(-rP_save(counter,case_dyn)+mult_f*(wPE*rE_save(counter,case_dyn)-wPP*rP_save(counter,case_dyn)-wPS*rS_save(counter,case_dyn)+xP)^power)/tauP;
    rS_save(counter+1,case_dyn)=rS_save(counter,case_dyn)+dt*(-rS_save(counter,case_dyn)+mult_f*(wSE*rE_save(counter,case_dyn)-wSP*rP_save(counter,case_dyn)-wSS*rS_save(counter,case_dyn)+xS)^power)/tauS;

end

%% Case 3 - high PV firing rate
case_dyn=6;
counter=0;

wEE=0.8; wEP=1; wPE=1; wPP=0.6; % E-PV parameters
wSS=0; 
wES=0.5; % inhibitory pathway
wPS=0.6; % disinhibitory pathway
wSE=0;% FB from E to SST
wSP=0; % FB from PV to SST

rS=2.5; rP=6.5; rE=5;

rE_save(1,case_dyn)=rE;
rP_save(1,case_dyn)=rP;
rS_save(1,case_dyn)=rS;

xE0=(rE_save(1,case_dyn)/mult_f)^(power^(-1))-(wEE*rE_save(1,case_dyn)-wEP*rP_save(1,case_dyn)-wES*rS_save(1,case_dyn)); 
xP0=(rP_save(1,case_dyn)/mult_f)^(power^(-1))-(wPE*rE_save(1,case_dyn)-wPP*rP_save(1,case_dyn)-wPS*rS_save(1,case_dyn)); 
xS0=(rS_save(1,case_dyn)/mult_f)^(power^(-1))-(wSE*rE_save(1,case_dyn)-wSP*rP_save(1,case_dyn)-wSS*rS_save(1,case_dyn));
xE=xE0; xP=xP0; xS=xS0;

for tt=dt:dt:end_sim
    counter=counter+1;
    
    if round(tt/dt)*dt==50
        xS=xS0+I_mod_S;
    end
    
    rE_save(counter+1,case_dyn)=rE_save(counter,case_dyn)+dt*(-rE_save(counter,case_dyn)+mult_f*(wEE*rE_save(counter,case_dyn)-wEP*rP_save(counter,case_dyn)-wES*rS_save(counter,case_dyn)+xE)^power)/tauE;
    rP_save(counter+1,case_dyn)=rP_save(counter,case_dyn)+dt*(-rP_save(counter,case_dyn)+mult_f*(wPE*rE_save(counter,case_dyn)-wPP*rP_save(counter,case_dyn)-wPS*rS_save(counter,case_dyn)+xP)^power)/tauP;
    rS_save(counter+1,case_dyn)=rS_save(counter,case_dyn)+dt*(-rS_save(counter,case_dyn)+mult_f*(wSE*rE_save(counter,case_dyn)-wSP*rP_save(counter,case_dyn)-wSS*rS_save(counter,case_dyn)+xS)^power)/tauS;

end

%% make figures
fig=figure
subplot(3,2,1)
hold on
plot(0:dt:end_sim,rE_save(:,1),'r')
plot(0:dt:end_sim,rP_save(:,1),'b')
plot(0:dt:end_sim,rS_save(:,1),'g')
hold off
legend('E','P','S')
xlabel('time')
ylabel('rX')
xlim([0 end_sim])
ylim([0 7])

subplot(3,2,2)
hold on
plot(0:dt:end_sim,rE_save(:,2),'r')
plot(0:dt:end_sim,rP_save(:,2),'b')
plot(0:dt:end_sim,rS_save(:,2),'g')
hold off
xlabel('time')
ylabel('rX')
xlim([0 end_sim])
ylim([0 7])

subplot(3,2,3)
hold on
plot(0:dt:end_sim,rE_save(:,3),'r')
plot(0:dt:end_sim,rP_save(:,3),'b')
plot(0:dt:end_sim,rS_save(:,3),'g')
hold off
xlabel('time')
ylabel('rX')
xlim([0 end_sim])
ylim([0 7])

subplot(3,2,4)
hold on
plot(0:dt:end_sim,rE_save(:,4),'r')
plot(0:dt:end_sim,rP_save(:,4),'b')
plot(0:dt:end_sim,rS_save(:,4),'g')
hold off
xlabel('time')
ylabel('rX')
xlim([0 end_sim])
ylim([0 7])

subplot(3,2,5)
hold on
plot(0:dt:end_sim,rE_save(:,5),'r')
plot(0:dt:end_sim,rP_save(:,5),'b')
plot(0:dt:end_sim,rS_save(:,5),'g')
hold off
xlabel('time')
ylabel('rX')
xlim([0 end_sim])
ylim([0 7])

subplot(3,2,6)
hold on
plot(0:dt:end_sim,rE_save(:,6),'r')
plot(0:dt:end_sim,rP_save(:,6),'b')
plot(0:dt:end_sim,rS_save(:,6),'g')
hold off
xlabel('time')
ylabel('rX')
xlim([0 end_sim])
ylim([0 7])

fig.Renderer='Painters';
%saveas(gcf,'Fig_1_cases.pdf')