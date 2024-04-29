%% Code to reproduce Fig. 7 of Bos, Miehl et al., 2024
clear all

%% calculate analytically the response matrix L and the eigenvalues (see Methods)
syms bEx bPx bSx wEEx wEPx wESx wPEx wPPx wPSx wSEx wSPx wSSx dEx dPx dSx

vec_syms=[bEx, bPx, bSx, wEEx, wEPx, wESx, wPEx, wPPx, wPSx, wSEx, wSPx, wSSx]; % vector of parameters (population gains and weights for each neuron/synapse)
W_sym=[[wEEx, -wEPx, -wESx];[wPEx, -wPPx,-wPSx];[wSEx,-wSPx,-wSSx]]; % weight matrix
D_sym=[[1/bEx,0,0];[0,1/bPx,0];[0,0,1/bSx]]; % D=B^(-1)
D_sym2=[[dEx,0,0];[0,dPx,0];[0,0,dSx]]; 

L_sym=inv(D_sym-W_sym); % response matrix L
L_sym2=inv(D_sym2-W_sym);
det_sym=det(D_sym-W_sym);
det_sym2=det(D_sym2-W_sym);

J_sym=-D_sym+W_sym; % Jacobian
EVs_sym=eig(J_sym); % Eigenvalues

%% define parameters 
% tuning parameters
theta_vec=[0:2:180];
sigma_th=20;
theta_p=90;
wFF_E=2;
wFF_P=2;

xE_base=1; % baseline input
xP_base=1;
xS=1.5;

xE_tuned=xE_base+wFF_E.*exp(-(theta_vec-theta_p).^2./sigma_th^2); % tuned input to E
xP_tuned=xP_base+wFF_P.*exp(-(theta_vec-theta_p).^2./sigma_th^2); % tuned input to PV

I_mod_S=-0.3;

% choose I/O function parameters f(x)=mult_f*x^power
power=2;
mult_f=1/4;

tauE=10; tauP=10; tauS=10; 
dt=0.01;
end_sim=100;

% cases: wES&wSE,wPS&wSP,wES&wSP,wPS&wSE, wPS&wSP&wSE,wES&wPS&wSE (see Fig.
% 7)
wES_vec=[0.8,0,0.8,0,0,0.5];
wPS_vec=[0,0.8,0,0.8,0.8,0.8];
wSE_vec=[0.5,0,0,0.5,0.5,0.5];
wSP_vec=[0,0.5,0.5,0,0.5,0];

wEE=0.8; wEP=0.5; wPE=1; wPP=0.6; % E-PV parameters
wSS=0; % SOM-to-SOM connection

for jj=1:length(wES_vec)

wES=wES_vec(jj); % inhibitory pathway
wPS=wPS_vec(jj); % disinhibitory pathway
wSE=wSE_vec(jj);% FB from E to SST
wSP=wSP_vec(jj); % FB from PV to SST

counter=0;

rE_save(1,1:length(theta_vec))=zeros(1,length(theta_vec));
rP_save(1,1:length(theta_vec))=zeros(1,length(theta_vec));
rS_save(1,1:length(theta_vec))=zeros(1,length(theta_vec));

for tt=dt:dt:end_sim
    counter=counter+1;

    rE_save(counter+1,:)=rE_save(counter,:)+dt.*(-rE_save(counter,:)+mult_f.*(wEE.*rE_save(counter,:)-wEP.*rP_save(counter,:)-wES.*rS_save(counter,:)+xE_tuned).^power)/tauE;
    rP_save(counter+1,:)=rP_save(counter,:)+dt.*(-rP_save(counter,:)+mult_f.*(wPE.*rE_save(counter,:)-wPP.*rP_save(counter,:)-wPS.*rS_save(counter,:)+xP_tuned).^power)/tauP;
    rS_save(counter+1,:)=rS_save(counter,:)+dt.*(-rS_save(counter,:)+mult_f.*(wSE.*rE_save(counter,:)-wSP.*rP_save(counter,:)-wSS.*rS_save(counter,:)+xS).^power)/tauS;

end

%% perform the SOM perturbation once the rates stabilized

bE_tuned=power*mult_f.*(wEE.*rE_save(end,:)-wEP.*rP_save(end,:)-wES.*rS_save(end,:)+xE_tuned).^(power-1); 
bP_tuned=power*mult_f.*(wPE.*rE_save(end,:)-wPP.*rP_save(end,:)-wPS.*rS_save(end,:)+xP_tuned).^(power-1); 
bS_tuned=power*mult_f.*(wSE.*rE_save(end,:)-wSP.*rP_save(end,:)-wSS.*rS_save(end,:)+xS).^(power-1);

det_tuned=1./bE_tuned*1./bP_tuned*1./bS_tuned - 1./bP_tuned*1./bS_tuned*wEE + 1./bE_tuned*1./bS_tuned*wPP + 1./bE_tuned*1./bP_tuned*wSS - 1./bS_tuned*wEE*wPP + 1./bS_tuned*wEP*wPE - 1./bP_tuned*wEE*wSS + 1./bP_tuned*wES*wSE + 1./bE_tuned*wPP*wSS - 1./bE_tuned*wPS*wSP - wEE*wPP*wSS + wEE*wPS*wSP + wEP*wPE*wSS - wEP*wPS*wSE - wES*wPE*wSP + wES*wPP*wSE;

L_ES_tuned=-(1./bP_tuned*wES - wEP*wPS + wES*wPP)./det_tuned;
L_PS_tuned=-(1./bE_tuned*wPS - wEE*wPS + wES*wPE)./det_tuned;
L_SS_tuned=(1./bE_tuned*1./bP_tuned - 1./bP_tuned*wEE + 1./bE_tuned*wPP - wEE*wPP + wEP*wPE)./det_tuned;

modS_rE_tuned=L_ES_tuned.*I_mod_S;
modS_rP_tuned=L_PS_tuned.*I_mod_S;
modS_rS_tuned=L_SS_tuned.*I_mod_S;

rE_fit=polyfit(rE_save(end,:),rE_save(end,:)+modS_rE_tuned,1);
rP_fit=polyfit(rP_save(end,:),rP_save(end,:)+modS_rP_tuned,1);
rS_fit=polyfit(rS_save(end,:),rS_save(end,:)+modS_rS_tuned,1);

rE_slope(jj)=rE_fit(1);
rP_slope(jj)=rP_fit(1);
rS_slope(jj)=rS_fit(1);

rE_int_norm(jj)=rE_fit(2)/max(rE_save(end,:)+modS_rE_tuned);
rP_int_norm(jj)=rP_fit(2)/max(rP_save(end,:)+modS_rP_tuned);
rS_int_norm(jj)=rS_fit(2)/max(rS_save(end,:)+modS_rS_tuned);

end


%%
fig1=figure;
subplot(3,3,1)
hold on
plot(theta_vec,rE_save(end,:),'r')
plot(theta_vec,rE_save(end,:)+modS_rE_tuned,'r--')
hold off
xlim([0 180])
ylim([0,3])
xlabel('Theta')
ylabel('rE')

subplot(3,3,2)
hold on
plot(theta_vec,rP_save(end,:),'b')
plot(theta_vec,rP_save(end,:)+modS_rP_tuned,'b--')
hold off
xlim([0 180])
ylim([0,2.5])
xlabel('Theta')
ylabel('rP')


subplot(3,3,3)
hold on
plot(theta_vec,rS_save(end,:),'g')
plot(theta_vec,rS_save(end,:)+modS_rS_tuned,'g--')
hold off
xlim([0 180])
ylim([0,2])
xlabel('Theta')
ylabel('rS')

subplot(3,3,4)
hold on
scatter(rE_save(end,:),rE_save(end,:)+modS_rE_tuned,20,'k','filled')
plot([0,3],[0,3],'k--')
plot([0,3],polyval(rE_fit,[0,3]),'k')
hold off
xlabel('rE before mod')
ylabel('rE after mod')
xlim([0,3])
ylim([0,3])

subplot(3,3,5)
hold on
scatter(rP_save(end,:),rP_save(end,:)+modS_rP_tuned,20,'k','filled')
plot([0,2.5],[0,2.5],'k--')
plot([0,2.5],polyval(rP_fit,[0,2.5]),'k')
hold off
xlabel('rP before mod')
ylabel('rP after mod')
xlim([0,2.5])
ylim([0,2.5])

subplot(3,3,6)
hold on
scatter(rS_save(end,:),rS_save(end,:)+modS_rS_tuned,20,'k','filled')
plot([0,2],[0,2],'k--')
plot([0,2],polyval(rS_fit,[0,2]),'k')
hold off
xlabel('rS before mod')
ylabel('rS after mod')
xlim([0,2])
ylim([0,2])

fig1.Renderer='Painters';

%saveas(gcf,'tuning_curves.pdf')

%%
fig2=figure;
hold on
scatter(rE_slope(1),rE_int_norm(1),100,'r',"square")
scatter(rE_slope(2),rE_int_norm(2),100,'r',"+")
scatter(rE_slope(3),rE_int_norm(3),100,'r')
scatter(rE_slope(5),rE_int_norm(5),100,'r','filled')
scatter(rE_slope(6),rE_int_norm(6),100,'r','diamond')

scatter(rP_slope(1),rP_int_norm(1),100,'b',"square")
scatter(rS_slope(1),rS_int_norm(1),100,'g',"square")


scatter(rP_slope(2),rP_int_norm(2),100,'b',"+")
scatter(rS_slope(2),rS_int_norm(2),100,'g',"+")

scatter(rP_slope(3),rP_int_norm(3),100,'b')
scatter(rS_slope(3),rS_int_norm(3),100,'g')


scatter(rP_slope(5),rP_int_norm(5),100,'b','filled')
scatter(rS_slope(5),rS_int_norm(5),100,'g','filled')

scatter(rP_slope(6),rP_int_norm(6),100,'b','diamond')
scatter(rS_slope(6),rS_int_norm(6),100,'g','diamond')
yline(0,'k--')
xline(1,'k--')
hold off
xlabel('Mult/Div Component')
ylabel('Add/Sub Component (norm)')
legend('wES&wSE','wPS&wSP','wES&wSP','wSP&wPS&wSE','wES&wPS&wSE')

fig2.Renderer='Painters';

%saveas(gcf,'tuning_measures.pdf')