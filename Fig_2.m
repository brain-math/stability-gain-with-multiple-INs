%% Code to reproduce Fig. 2 of Bos, Miehl et al., 2025

% calculate analytically the response matrix L and the eigenvalues (see Methods)
syms bEx bPx bSx wEEx wEPx wESx wPEx wPPx wPSx wSEx wSPx wSSx dEx dPx dSx

vec_syms=[bEx, bPx, bSx, wEEx, wEPx, wESx, wPEx, wPPx, wPSx, wSEx, wSPx, wSSx]; % vector of parameters (population gains and weights for each neuron/synapse)
W_sym=[[wEEx, -wEPx, -wESx];[wPEx, -wPPx,-wPSx];[wSEx,-wSPx,-wSSx]]; % weight matrix
D_sym=[[1/bEx,0,0];[0,1/bPx,0];[0,0,1/bSx]]; % D=B^(-1)

L_sym=inv(D_sym-W_sym); % response matrix L
det_sym=det(D_sym-W_sym);

J_sym=-D_sym+W_sym; % Jacobian
EVs_sym=eig(J_sym); % Eigenvalues



% Parameters for neuron timescales and I/O function
tauE=10; tauP=10; tauS=10; 
mult_f=1/4;
power=2; 

dt=0.01; % time step
end_sim=750; % length of simulation

w_vec=[0.8,0.5,0,1,0.6,0.8,0,0,0;0.8,0.5,0,1,0.6,0.8,0,0.2,0];
r_vec=[3.5,6,2;5,2,3];
I_mod_vec=[0.3;-0.3];

for kk=1:2 % generate the two cases of Fig. 2

wEE=w_vec(kk,1); wEP=w_vec(kk,2); wES=w_vec(kk,3); wPE=w_vec(kk,4); wPP=w_vec(kk,5); wPS=w_vec(kk,6); wSE=w_vec(kk,7); wSP=w_vec(kk,8); wSS=w_vec(kk,9);

rE0=r_vec(kk,1); rP0=r_vec(kk,2); rS0=r_vec(kk,3);
rE=rE0; rP=rP0; rS=rS0;

% calculate the inputs for each population from the firing rate steady
% states rX=mult_f*(wXE*rE-wXP*rP-wXS*rS+xX0)^power for X=[E,P,S]
xE0=(rE/mult_f)^(power^(-1))-(wEE*rE-wEP*rP-wES*rS); 
xP0=(rP/mult_f)^(power^(-1))-(wPE*rE-wPP*rP-wPS*rS); 
xS0=(rS/mult_f)^(power^(-1))-(wSE*rE-wSP*rP-wSS*rS);
xE=xE0; xP=xP0; xS=xS0;

I_stim_E=0.3; I_stim_P=0.3; I_mod_S=I_mod_vec(kk); % stimulus and modulation perturbations

counter=0;

%% calculate analytical solutions (symbols in Fig. 2)

%% calculate initial condition

% calculate steady state gains
bE=power*mult_f*(wEE*rE-wEP*rP-wES*rS+xE)^(power-1); 
bP=power*mult_f*(wPE*rE-wPP*rP-wPS*rS+xP)^(power-1); 
bS=power*mult_f*(wSE*rE-wSP*rP-wSS*rS+xS)^(power-1);

vec_nums=[bE,bP,bS,wEE,wEP,wES,wPE,wPP,wPS,wSE,wSP,wSS]; %vector to solve equation based on parameters

EVs_num(:,1+4*(kk-1))=double(subs(EVs_sym,vec_syms,vec_nums));

%% Calculate SST perturbation after initial condition

modS_num(1,1)=double(subs(L_sym(1,3),vec_syms,vec_nums)*I_mod_S); % E
modS_num(2,1)=double(subs(L_sym(2,3),vec_syms,vec_nums)*I_mod_S); % P
modS_num(3,1)=double(subs(L_sym(3,3),vec_syms,vec_nums)*I_mod_S); % S

rE_mod(kk,1)=rE0+modS_num(1,1); rP_mod(kk,1)=rP0+modS_num(2,1); rS_mod(kk,1)=rS0+modS_num(3,1);

bE_mod=power*mult_f*(wEE*rE_mod(kk,1)-wEP*rP_mod(kk,1)-wES*rS_mod(kk,1)+xE)^(power-1); 
bP_mod=power*mult_f*(wPE*rE_mod(kk,1)-wPP*rP_mod(kk,1)-wPS*rS_mod(kk,1)+xP)^(power-1); 
bS_mod=power*mult_f*(wSE*rE_mod(kk,1)-wSP*rP_mod(kk,1)-wSS*rS_mod(kk,1)+xS)^(power-1);

vec_nums_mod=[bE_mod,bP_mod,bS_mod,wEE,wEP,wES,wPE,wPP,wPS,wSE,wSP,wSS];

EVs_num(:,2+4*(kk-1))=double(subs(EVs_sym,vec_syms,vec_nums_mod));

%% Calculate E & PV stimulation

stim_init(1,1)=double(subs(L_sym(1,1),vec_syms,vec_nums))*I_stim_E+double(subs(L_sym(1,2),vec_syms,vec_nums))*I_stim_P; % effect of stim on E
stim_init(2,1)=double(subs(L_sym(2,1),vec_syms,vec_nums))*I_stim_E+double(subs(L_sym(2,2),vec_syms,vec_nums))*I_stim_P; % effect of stim on P
stim_init(3,1)=double(subs(L_sym(3,1),vec_syms,vec_nums))*I_stim_E+double(subs(L_sym(3,2),vec_syms,vec_nums))*I_stim_P; % effect of stim on S

rE_stim_init(kk,1)=rE0+stim_init(1,1); rP_stim_init(kk,1)=rP0+stim_init(2,1); rS_stim_init(kk,1)=rS0+stim_init(3,1);

bE_stim_init=power*mult_f*(wEE*rE_stim_init(kk,1)-wEP*rP_stim_init(kk,1)-wES*rS_stim_init(kk,1)+xE)^(power-1); 
bP_stim_init=power*mult_f* (wPE*rE_stim_init(kk,1)-wPP*rP_stim_init(kk,1)-wPS*rS_stim_init(kk,1)+xP)^(power-1); 
bS_stim_init=power*mult_f*(wSE*rE_stim_init(kk,1)-wSP*rP_stim_init(kk,1)-wSS*rS_stim_init(kk,1)+xS)^(power-1);

vec_nums_stim_init=[bE_stim_init,bP_stim_init,bS_stim_init,wEE,wEP,wES,wPE,wPP,wPS,wSE,wSP,wSS];

EVs_num(:,3+4*(kk-1))=double(subs(EVs_sym,vec_syms,vec_nums_stim_init));

%% Calculate E & PV stimulation after SST modulation

stim_mod(1,1)=double(subs(L_sym(1,1),vec_syms,vec_nums_mod))*I_stim_E+double(subs(L_sym(1,2),vec_syms,vec_nums_mod))*I_stim_P; % effect of stim on E
stim_mod(2,1)=double(subs(L_sym(2,1),vec_syms,vec_nums_mod))*I_stim_E+double(subs(L_sym(2,2),vec_syms,vec_nums_mod))*I_stim_P; % effect of stim on P
stim_mod(3,1)=double(subs(L_sym(3,1),vec_syms,vec_nums_mod))*I_stim_E+double(subs(L_sym(3,2),vec_syms,vec_nums_mod))*I_stim_P; % effect of stim on S

rE_stim_mod(kk,1)=rE_mod(kk,1)+stim_mod(1,1); rP_stim_mod(kk,1)=rP_mod(kk,1)+stim_mod(2,1); rS_stim_mod(kk,1)=rS_mod(kk,1)+stim_mod(3,1);

bE_stim_mod=power*mult_f*(wEE*rE_stim_mod(kk,1)-wEP*rP_stim_mod(kk,1)-wES*rS_stim_mod(kk,1)+xE)^(power-1); 
bP_stim_mod=power*mult_f*(wPE*rE_stim_mod(kk,1)-wPP*rP_stim_mod(kk,1)-wPS*rS_stim_mod(kk,1)+xP)^(power-1); 
bS_stim_mod=power*mult_f*(wSE*rE_stim_mod(kk,1)-wSP*rP_stim_mod(kk,1)-wSS*rS_stim_mod(kk,1)+xS)^(power-1);

vec_nums_stim_mod=[bE_stim_mod,bP_stim_mod,bS_stim_mod,wEE,wEP,wES,wPE,wPP,wPS,wSE,wSP,wSS];

EVs_num(:,4+4*(kk-1))=double(subs(EVs_sym,vec_syms,vec_nums_stim_mod));


%% numerical simulation

rE_save(1,1+2*(kk-1))=rE0; rE_save(1,2+2*(kk-1))=rE0;
rP_save(1,1+2*(kk-1))=rP0; rP_save(1,2+2*(kk-1))=rP0;
rS_save(1,1+2*(kk-1))=rS0; rS_save(1,2+2*(kk-1))=rS0;

xS0_2=xS0;

for tt=dt:dt:end_sim
    counter=counter+1;
    
    if round(tt/dt)*dt==50
        xS=xS0+I_mod_S;

    elseif round(tt/dt)*dt==350
        xE=xE0+I_stim_E;
        xP=xP0+I_stim_P;
    end

    rE_save(counter+1,1+2*(kk-1))=rE_save(counter,1+2*(kk-1))+dt*(-rE_save(counter,1+2*(kk-1))+mult_f*(wEE*rE_save(counter,1+2*(kk-1))-wEP*rP_save(counter,1+2*(kk-1))-wES*rS_save(counter,1+2*(kk-1))+xE)^power)/tauE;
    rP_save(counter+1,1+2*(kk-1))=rP_save(counter,1+2*(kk-1))+dt*(-rP_save(counter,1+2*(kk-1))+mult_f*(wPE*rE_save(counter,1+2*(kk-1))-wPP*rP_save(counter,1+2*(kk-1))-wPS*rS_save(counter,1+2*(kk-1))+xP)^power)/tauP;
    rS_save(counter+1,1+2*(kk-1))=rS_save(counter,1+2*(kk-1))+dt*(-rS_save(counter,1+2*(kk-1))+mult_f*(wSE*rE_save(counter,1+2*(kk-1))-wSP*rP_save(counter,1+2*(kk-1))-wSS*rS_save(counter,1+2*(kk-1))+xS)^power)/tauS;

    rE_save(counter+1,2+2*(kk-1))=rE_save(counter,2+2*(kk-1))+dt*(-rE_save(counter,2+2*(kk-1))+mult_f*(wEE*rE_save(counter,2+2*(kk-1))-wEP*rP_save(counter,2+2*(kk-1))-wES*rS_save(counter,2+2*(kk-1))+xE)^power)/tauE;
    rP_save(counter+1,2+2*(kk-1))=rP_save(counter,2+2*(kk-1))+dt*(-rP_save(counter,2+2*(kk-1))+mult_f*(wPE*rE_save(counter,2+2*(kk-1))-wPP*rP_save(counter,2+2*(kk-1))-wPS*rS_save(counter,2+2*(kk-1))+xP)^power)/tauP;
    rS_save(counter+1,2+2*(kk-1))=rS_save(counter,2+2*(kk-1))+dt*(-rS_save(counter,2+2*(kk-1))+mult_f*(wSE*rE_save(counter,2+2*(kk-1))-wSP*rP_save(counter,2+2*(kk-1))-wSS*rS_save(counter,2+2*(kk-1))+xS0_2)^power)/tauS;

end

%%
figure(kk)
subplot(3,3,[2 3])
hold on
plot(0:dt:end_sim,rE_save(:,2+2*(kk-1)),'r--')
plot(0:dt:end_sim,rP_save(:,2+2*(kk-1)),'b--')
plot(0:dt:end_sim,rS_save(:,2+2*(kk-1)),'g--')
plot(0:dt:end_sim,rE_save(:,1+2*(kk-1)),'r')
plot(0:dt:end_sim,rP_save(:,1+2*(kk-1)),'b')
plot(0:dt:end_sim,rS_save(:,1+2*(kk-1)),'g')
scatter(0,rE0,'r')
scatter(0,rP0,'b')
scatter(0,rS0,'g')
scatter(300,rE_mod(kk,1),'r','filled')
scatter(300,rP_mod(kk,1),'b','filled')
scatter(300,rS_mod(kk,1),'g','filled')
scatter(end_sim,rE_stim_mod(kk,1),'r','filled',"square")
scatter(end_sim,rP_stim_mod(kk,1),'b','filled',"square")
scatter(end_sim,rS_stim_mod(kk,1),'g','filled',"square")
scatter(end_sim,rE_stim_init(kk,1),'r',"*")
scatter(end_sim,rP_stim_init(kk,1),'b',"*")
scatter(end_sim,rS_stim_init(kk,1),'g',"*")
xline(50,'k--')
xline(350,'k--')
hold off
xlabel('Time (ms)')
ylabel('rX')
ylim([0 8])
xlim([0 end_sim])

subplot(3,3,4)
hold on
scatter(rP0,rE0,'k')
scatter(rP_mod(kk,1),rE_mod(kk,1),'k','filled')
plot([rP0 rP_mod(kk,1)],[rE0 rE0],'b')
plot([rP0 rP0],[rE0 rE_mod(kk,1)],'r')
plot([rP0 rP_mod(kk,1)],[rE0 rE_mod(kk,1)],'k--')
hold off
xlabel('rP (Hz)')
ylabel('rE (Hz)')

if kk==1
    xlim([5.5 7.75])
    ylim([3 5.25])
elseif kk==2
    xlim([0 4])
    ylim([3 6.5])
end


subplot(3,3,5)
hold on
scatter(rP0,rE0,'k')
scatter(rP_mod(kk,1),rE_mod(kk,1),'k','filled')
scatter(rP_stim_mod(kk,1),rE_stim_mod(kk,1),'k','filled',"square")
scatter(rP_stim_init(kk,1),rE_stim_init(kk,1),'k',"*")
plot([rP_mod(kk,1) rP_stim_mod(kk,1)],[rE_mod(kk,1) rE_mod(kk,1)],'b')
plot([rP_mod(kk,1) rP_mod(kk,1)],[rE_mod(kk,1) rE_stim_mod(kk,1)],'r')
plot([rP_mod(kk,1) rP_stim_mod(kk,1)],[rE_mod(kk,1) rE_stim_mod(kk,1)],'k--')
plot([rP0 rP_stim_init(kk,1)],[rE0 rE_stim_init(kk,1)],'k--')
plot([rP0 rP_stim_init(kk,1)],[rE0 rE0],'b')
plot([rP0 rP0],[rE0 rE_stim_init(kk,1)],'r')
hold off
xlabel('rP (Hz)')
ylabel('rE (Hz)')
if kk==1
    xlim([5.5 7.75])
    ylim([3 5.25])
elseif kk==2
    xlim([0 4])
    ylim([3 6.5])
end
%title(['g no mod= ' num2str(round(rE_stim_init(kk,1)-rE0,2)) ', g mod= ' num2str(round(rE_stim_mod(kk,1)-rE_mod(kk,1),2))])

subplot(3,3,6)
hold on
scatter(real(EVs_num(:,1+4*(kk-1))),imag(EVs_num(:,1+4*(kk-1))),'k')
scatter(real(EVs_num(:,2+4*(kk-1))),imag(EVs_num(:,2+4*(kk-1))),'k','filled')
xline(0,'k--')
yline(0,'k--')
hold off
xlabel('Real')
ylabel('Imag')
%title(['sta no mod= ' num2str(round(max(real(EVs_num(:,1+4*(kk-1)))),2)) ', sta mod= ' num2str(round(max(real(EVs_num(:,2))),2))])
if kk==1
    xlim([-0.8 0.05])
    ylim([-0.5 0.5])
elseif kk==2
    xlim([-0.3 0.05])
    ylim([-0.5 0.5])
end
end
