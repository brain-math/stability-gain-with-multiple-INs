%% Code to reproduce Fig. 3 of Bos, Miehl et al., 2025

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
wEE=0.8; wEP=0.5; wPE=1; wPP=0.6; % E-PV parameters
wSS=0; 
wES=0; % inh pathway
wPS=0; % disinh pathway
wSE=0;% FB from E to SST
wSP=0; % FB from PV to SST

rS0=2; % SST firing rate

step_rX=0.01;
rE_vec=[10:-step_rX:0.5]; rP_vec=[0.5:step_rX:10]; % rP-rE space

I_stim_E=0.3; I_stim_P=0.3; I_mod_S=0.3; % perturbation strengths

% choose I/O function parameters f(x)=mult_f*x^power
power=2;
mult_f=1/4;

gain_num_cell=cell(3,1);
gain_num_mod_p_cell=cell(3,1);
gain_num_mod_m_cell=cell(3,1);
modS_rE_num_cell=cell(3,1);
modS_rP_num_cell=cell(3,1);
modS_rS_num_cell=cell(3,1);
maxEVs_num_cell=cell(3,1);
maxEVs_num_mod_p_cell=cell(3,1);
maxEVs_num_mod_m_cell=cell(3,1);


% initialize arrays
gain_num=zeros(length(rE_vec),length(rP_vec));
gain_num_mod_p=zeros(length(rE_vec),length(rP_vec));
gain_num_mod_m=zeros(length(rE_vec),length(rP_vec));
modS_rE_num=zeros(length(rE_vec),length(rP_vec));
modS_rP_num=zeros(length(rE_vec),length(rP_vec));
modS_rS_num=zeros(length(rE_vec),length(rP_vec));
maxEVs_num=zeros(length(rE_vec),length(rP_vec));
maxEVs_num_mod_p=zeros(length(rE_vec),length(rP_vec));
maxEVs_num_mod_m=zeros(length(rE_vec),length(rP_vec));


%% loop through E&PV rates
for jj=1:length(rE_vec)
    for jj2=1:length(rP_vec)

    rE=rE_vec(jj); rP=rP_vec(jj2); rS=rS0;

    % calculate the inputs for each population from the firing rate steady
    % states rX=mult_f*(wXE*rE-wXP*rP-wXS*rS+xX0)^power for X=[E,P,S]
    xE=(rE/mult_f)^(power^(-1))-(wEE*rE-wEP*rP-wES*rS); 
    xP=(rP/mult_f)^(power^(-1))-(wPE*rE-wPP*rP-wPS*rS); 
    xS=(rS/mult_f)^(power^(-1))-(wSE*rE-wSP*rP-wSS*rS);

    % calculate derivatives (bX)
    % use f'(x)=power*mult_f*(...)^(power-1)
    bE=power*mult_f*(wEE*rE-wEP*rP-wES*rS+xE)^(power-1); 
    bP=power*mult_f*(wPE*rE-wPP*rP-wPS*rS+xP)^(power-1); 
    bS=power*mult_f*(wSE*rE-wSP*rP-wSS*rS+xS)^(power-1);

    % response matrix terms (calculated in det_sym2, L_sym2 and J_sym)
    det_num=1/bE*1/bP*1/bS - 1/bP*1/bS*wEE + 1/bE*1/bS*wPP + 1/bE*1/bP*wSS - 1/bS*wEE*wPP + 1/bS*wEP*wPE - 1/bP*wEE*wSS + 1/bP*wES*wSE + 1/bE*wPP*wSS - 1/bE*wPS*wSP - wEE*wPP*wSS + wEE*wPS*wSP + wEP*wPE*wSS - wEP*wPS*wSE - wES*wPE*wSP + wES*wPP*wSE;
    L_EE=(1/bP*1/bS + 1/bS*wPP + 1/bP*wSS + wPP*wSS - wPS*wSP)/det_num;
    L_EP=-(1/bS*wEP + wEP*wSS - wES*wSP)/det_num;
    L_ES=-(1/bP*wES - wEP*wPS + wES*wPP)/det_num;
    L_PS=-(1/bE*wPS - wEE*wPS + wES*wPE)/det_num;
    L_SS=(1/bE*1/bP - 1/bP*wEE + 1/bE*wPP - wEE*wPP + wEP*wPE)/det_num;
    J_num=[[wEE - 1/bE,-wEP,-wES];[wPE,- wPP - 1/bP,-wPS];[wSE,-wSP,- wSS - 1/bS]];

    % calculate gain and max EV
    gain_num(jj,jj2)=L_EE*I_stim_E+L_EP*I_stim_P;
    maxEVs_num(jj,jj2)=max(real(eig(J_num)));

    modS_rE_num(jj,jj2)=L_ES*I_mod_S;
    modS_rP_num(jj,jj2)=L_PS*I_mod_S;
    modS_rS_num(jj,jj2)=L_SS*I_mod_S;

    % calculate stability and gain after positive SST modulation
    bE_mod_p=power*mult_f*(wEE*(rE+modS_rE_num(jj,jj2))-wEP*(rP+modS_rP_num(jj,jj2))-wES*(rS+modS_rS_num(jj,jj2))+xE)^(power-1); 
    bP_mod_p=power*mult_f*(wPE*(rE+modS_rE_num(jj,jj2))-wPP*(rP+modS_rP_num(jj,jj2))-wPS*(rS+modS_rS_num(jj,jj2))+xP)^(power-1); 
    bS_mod_p=power*mult_f*(wSE*(rE+modS_rE_num(jj,jj2))-wSP*(rP+modS_rP_num(jj,jj2))-wSS*(rS+modS_rS_num(jj,jj2))+xS)^(power-1);

    det_num_mod_p=1/bE_mod_p*1/bP_mod_p*1/bS_mod_p - 1/bP_mod_p*1/bS_mod_p*wEE + 1/bE_mod_p*1/bS_mod_p*wPP + 1/bE_mod_p*1/bP_mod_p*wSS - 1/bS_mod_p*wEE*wPP + 1/bS_mod_p*wEP*wPE - 1/bP_mod_p*wEE*wSS + 1/bP_mod_p*wES*wSE + 1/bE_mod_p*wPP*wSS - 1/bE_mod_p*wPS*wSP - wEE*wPP*wSS + wEE*wPS*wSP + wEP*wPE*wSS - wEP*wPS*wSE - wES*wPE*wSP + wES*wPP*wSE;
    L_EE_mod_p=(1/bP_mod_p*1/bS_mod_p + 1/bS_mod_p*wPP + 1/bP_mod_p*wSS + wPP*wSS - wPS*wSP)/det_num_mod_p;
    L_EP_mod_p=-(1/bS_mod_p*wEP + wEP*wSS - wES*wSP)/det_num_mod_p;

    J_num_mod_p=[[wEE - 1/bE_mod_p,-wEP,-wES];[wPE,- wPP - 1/bP_mod_p,-wPS];[wSE,-wSP,- wSS - 1/bS_mod_p]];

    gain_num_mod_p(jj,jj2)=L_EE_mod_p*I_stim_E+L_EP_mod_p*I_stim_P;
    
    if sum(sum(isnan(J_num_mod_p)))>0
         maxEVs_num_mod_p(jj,jj2)=NaN;
    else
        maxEVs_num_mod_p(jj,jj2)=max(real(eig(J_num_mod_p)));
    end

    % calculate stability and gain after negative SST modulation
    bE_mod_m=power*mult_f*(wEE*(rE-modS_rE_num(jj,jj2))-wEP*(rP-modS_rP_num(jj,jj2))-wES*(rS-modS_rS_num(jj,jj2))+xE)^(power-1); 
    bP_mod_m=power*mult_f*(wPE*(rE-modS_rE_num(jj,jj2))-wPP*(rP-modS_rP_num(jj,jj2))-wPS*(rS-modS_rS_num(jj,jj2))+xP)^(power-1); 
    bS_mod_m=power*mult_f*(wSE*(rE-modS_rE_num(jj,jj2))-wSP*(rP-modS_rP_num(jj,jj2))-wSS*(rS-modS_rS_num(jj,jj2))+xS)^(power-1);

    det_num_mod_m=1/bE_mod_m*1/bP_mod_m*1/bS_mod_m - 1/bP_mod_m*1/bS_mod_m*wEE + 1/bE_mod_m*1/bS_mod_m*wPP + 1/bE_mod_m*1/bP_mod_m*wSS - 1/bS_mod_m*wEE*wPP + 1/bS_mod_m*wEP*wPE - 1/bP_mod_m*wEE*wSS + 1/bP_mod_m*wES*wSE + 1/bE_mod_m*wPP*wSS - 1/bE_mod_m*wPS*wSP - wEE*wPP*wSS + wEE*wPS*wSP + wEP*wPE*wSS - wEP*wPS*wSE - wES*wPE*wSP + wES*wPP*wSE;
    L_EE_mod_m=(1/bP_mod_m*1/bS_mod_m + 1/bS_mod_m*wPP + 1/bP_mod_m*wSS + wPP*wSS - wPS*wSP)/det_num_mod_m;
    L_EP_mod_m=-(1/bS_mod_m*wEP + wEP*wSS - wES*wSP)/det_num_mod_m;

    J_num_mod_m=[[wEE - 1/bE_mod_m,-wEP,-wES];[wPE,- wPP - 1/bP_mod_m,-wPS];[wSE,-wSP,- wSS - 1/bS_mod_m]];

    gain_num_mod_m(jj,jj2)=L_EE_mod_m*I_stim_E+L_EP_mod_m*I_stim_P;
    
   if sum(sum(isnan(J_num_mod_m)))>0
         maxEVs_num_mod_m(jj,jj2)=NaN;
   else
        maxEVs_num_mod_m(jj,jj2)=max(real(eig(J_num_mod_m)));
   end

    end
end

gain_num_cell{1,1}=gain_num;
gain_num_mod_p_cell{1,1}=gain_num_mod_p;
gain_num_mod_m_cell{1,1}=gain_num_mod_m;
modS_rE_num_cell{1,1}=modS_rE_num;
modS_rP_num_cell{1,1}=modS_rP_num;
modS_rS_num_cell{1,1}=modS_rS_num;
maxEVs_num_cell{1,1}=maxEVs_num;
maxEVs_num_mod_p_cell{1,1}=maxEVs_num_mod_p;
maxEVs_num_mod_m_cell{1,1}=maxEVs_num_mod_m;


%% eigenvalue panel

rS=2; rP_vec2=[0.5:0.05:10]; rE=5;

bE=mult_f*power*(rE/mult_f)^((power-1)/power);
bS=mult_f*power*(rS/mult_f)^((power-1)/power);

for kk=1:length(rP_vec2)
    
    rP=rP_vec2(kk);
    
    bP=mult_f*power*(rP/mult_f)^((power-1)/power);
   
    vec_nums=[bE,bP,bS,wEE,wEP,wES,wPE,wPP,wPS,wSE,wSP,wSS];

    EVs_bP(1:3,kk)=eig(double(subs(J_sym,vec_syms,vec_nums)));


end

%% make figures

vec=120:-20:0; % define colors for heatmap
NNN = 128;
hex=['#ffffff','#fcbba1','#fc9272','#fb6a4a','#ef3b2c','#cb181d','#99000d']';
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
map = interp1(vec,raw,linspace(120,0,NNN),'pchip');

rE_mat=repmat(rE_vec',1,length(rE_vec));
rP_mat=repmat(rP_vec,length(rP_vec),1);

% figure 1
figure
    
subplot(3,3,2) % stability, norm
mask_threshold=-0.05;
maxEVs_num_set0=min(maxEVs_num_cell{1,1},0);
mask01=maxEVs_num_cell{1,1};
mask01(mask01>mask_threshold)=0;
mask01(mask01<mask_threshold)=1;
norm_maxEVs_num_set0=mask01.*maxEVs_num_set0./min(min(maxEVs_num_set0));


imagesc(norm_maxEVs_num_set0)
colormap(map)
title('max EV, norm')
xlabel('rP (Hz)')
ylabel('rE (Hz)')
axis square

xticklabels = linspace(min(rP_vec),max(rP_vec),3);
xticks = linspace(1, size(maxEVs_num_cell{1,1}, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

yticklabels = linspace(max(rE_vec),min(rE_vec),3);
yticks = linspace(1, size(maxEVs_num_cell{1,1}, 2), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)


subplot(3,3,1) % Gain, norm

gain_num_masked=mask01.*gain_num_cell{1,1};
norm_gain_num_masked=gain_num_masked./max(max(gain_num_masked));

imagesc(norm_gain_num_masked)
%colorbar
colormap(map)
title('Gain, norm')
xlabel('rP (Hz)')
ylabel('rE (Hz)')
axis square

xticklabels = linspace(min(rP_vec),max(rP_vec),3);
xticks = linspace(1, size(gain_num, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

yticklabels = linspace(max(rE_vec),min(rE_vec),3);
yticks = linspace(1, size(gain_num, 2), numel(yticklabels));
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)




% isolines

% ISN boundary
bE_mat=power*mult_f*((rE_mat/mult_f).^(1/power)).^(power-1);
isn_mat=1-bE_mat.*wEE; %<0 ISN regime, >0 non-ISN regime
[val_ISN,loc_ISN]=max(isn_mat>=0,[],1); % find me where isn_mat switches to non-ISN regime

subplot(3,3,5) % stability isolines + stability boundary

[val_stab,loc_stab]=max(norm_maxEVs_num_set0>0,[],1); % find where the system becomes stable
loc_stab(val_stab==0)=NaN; % if a column is only unstable
loc_stab(loc_stab==1)=NaN; % if a column is only stable
[a,loc_noNaNs]=find(~isnan(loc_stab));
loc_stab(isnan(loc_stab))=[];

iso_num_EVs=[0.2,0.4,0.6,0.8];

for uu=1:length(iso_num_EVs)
    [val_iso(uu,:),loc_iso(uu,:)]=max(norm_maxEVs_num_set0>iso_num_EVs(uu),[],1);
    loc_iso(uu,val_iso(uu,:)==0)=NaN;
    loc_iso(uu,loc_iso(uu,:)==1)=NaN;
end


hold on
plot(rP_vec,rE_vec(loc_ISN),'k')
for uu2=1:length(iso_num_EVs)
    dummy_loc_iso=loc_iso(uu2,:);
    [a,loc_iso_noNaNs]=find(~isnan(dummy_loc_iso));
    dummy_loc_iso(isnan(dummy_loc_iso))=[];
    plot(rP_vec(loc_iso_noNaNs),rE_vec(dummy_loc_iso))
end
plot(rP_vec(loc_noNaNs),rE_vec(loc_stab),'k')
title('Stab, iso')
hold off
xlim([min(rP_vec) max(rP_vec)])
ylim([min(rE_vec) max(rE_vec)])
xlabel('rP')
ylabel('rE')
axis square


subplot(3,3,4) % gain insolines

iso_num_gain=[0.02,0.05,0.1,0.2];

gain_num_masked=mask01.*gain_num_cell{1,1};
norm_gain_num_masked=gain_num_masked./max(max(gain_num_masked));

for uu=1:length(iso_num_gain)
    norm_gain_num_masked_NaN=norm_gain_num_masked;
    norm_gain_num_masked_NaN(norm_gain_num_masked_NaN==0)=NaN;
    [val_iso_gain(uu,:),loc_iso_gain(uu,:)]=max(norm_gain_num_masked_NaN<iso_num_gain(uu),[],1);
    loc_iso_gain(uu,val_iso_gain(uu,:)==0)=NaN;
    loc_iso_gain(uu,loc_iso_gain(uu,:)==1)=NaN;
end

hold on
plot(rP_vec,rE_vec(loc_ISN),'k')
for uu2=1:length(iso_num_gain)
    dummy_loc_iso_gain=loc_iso_gain(uu2,:);
    [a,loc_iso_noNaNs_gain]=find(~isnan(dummy_loc_iso_gain));
    dummy_loc_iso_gain(isnan(dummy_loc_iso_gain))=[];
    plot(rP_vec(loc_iso_noNaNs_gain),rE_vec(dummy_loc_iso_gain))
end
plot(rP_vec(loc_noNaNs),rE_vec(loc_stab),'k')
hold off
title('Gain, iso')
xlim([min(rP_vec) max(rP_vec)])
ylim([min(rE_vec) max(rE_vec)])
xlabel('rP')
ylabel('rE')
axis square

fig1.Renderer='Painters';

%% eigenvalue panel
fig2=figure
hex_blue=['#eff3ff','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#084594']'; %blue
raw_blue = sscanf(hex_blue','#%2x%2x%2x',[3,size(hex_blue,1)]).' / 255;
map_blue = interp1(vec,raw_blue,linspace(120,0,NNN),'pchip');

subplot(3,3,1)
hold on
yline(0,'k--')
patch([real(EVs_bP(1,:)) nan],[imag(EVs_bP(1,:)) nan],[rP_vec2 nan],[rP_vec2 nan], 'edgecolor', 'interp','Linewidth',2); 
patch([real(EVs_bP(2,:)) nan],[imag(EVs_bP(2,:)) nan],[rP_vec2 nan],[rP_vec2 nan], 'edgecolor', 'interp','Linewidth',2); 
patch([real(EVs_bP(3,:)) nan],[imag(EVs_bP(3,:)) nan],[rP_vec2 nan],[rP_vec2 nan], 'edgecolor', 'interp','Linewidth',2); 
scatter(real(EVs_bP(1,1)),imag(EVs_bP(1,1)),'k','filled')
scatter(real(EVs_bP(2,1)),imag(EVs_bP(2,1)),'k','filled')
scatter(real(EVs_bP(3,1)),imag(EVs_bP(3,1)),'k','filled')
xline(0,'k--')
hold off
c=colorbar;
colormap(map_blue);
ylim([-0.4 0.4])
xlabel('real')
ylabel('imag')
title('EVs of bE example rE=5')
axis square

fig2.Renderer='Painters';