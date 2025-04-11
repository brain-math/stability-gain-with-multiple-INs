%% Code to reproduce Fig. 5 - figure supplement 2 Bos, Miehl et al., 2025

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
% Parameters Fig. S2A
wEE=0.8; wEP=0.5; wPE=1; wPP=0.6; % E-PV parameters
wSS=0; 
wES=0.8; % inhibitory pathway
wPS=0; % disinhibitory pathway
wSE=0.2;% FB from E to SST
wSP=0; % FB from PV to SST

% Parameters Fig. S2B:
%wPS=0.8; wSP=0.2; wES=0; wSE=0;

rS0_vec=[0.5:0.5:10]; % vector of different SST rates

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

for zz=1:length(rS0_vec) % SST FR loop

    rS0=rS0_vec(zz);


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
    
    gain_num_cell{zz,1}=gain_num;
    gain_num_mod_p_cell{zz,1}=gain_num_mod_p;
    gain_num_mod_m_cell{zz,1}=gain_num_mod_m;
    modS_rE_num_cell{zz,1}=modS_rE_num;
    modS_rP_num_cell{zz,1}=modS_rP_num;
    modS_rS_num_cell{zz,1}=modS_rS_num;
    maxEVs_num_cell{zz,1}=maxEVs_num;
    maxEVs_num_mod_p_cell{zz,1}=maxEVs_num_mod_p;
    maxEVs_num_mod_m_cell{zz,1}=maxEVs_num_mod_m;
end



%%
fig2=figure

for zz2=1:length(rS0_vec)

mask_threshold=-0.05;

maxEVs_num_set0_2=min(maxEVs_num_cell{zz2,1},0); % get me unstable points
mask01_2=maxEVs_num_cell{zz2,1};
mask01_2(mask01_2>mask_threshold)=NaN;
mask01_2(mask01_2<mask_threshold)=1; % have a mask which is NaN for unstable and 1 for stable

maxEVs_num_set0_2_mod_p=min(maxEVs_num_mod_p_cell{zz2,1},0); % same as above for positive modulation
mask01_2_mod_p=maxEVs_num_mod_p_cell{zz2,1};
mask01_2_mod_p(mask01_2_mod_p>mask_threshold)=NaN;
mask01_2_mod_p(mask01_2_mod_p<mask_threshold)=1;

maxEVs_num_set0_2_mod_m=min(maxEVs_num_mod_m_cell{zz2,1},0); % same as above for negative modulation
mask01_2_mod_m=maxEVs_num_mod_m_cell{zz2,1};
mask01_2_mod_m(mask01_2_mod_m>mask_threshold)=NaN;
mask01_2_mod_m(mask01_2_mod_m<mask_threshold)=1;

diff_maxEVs_num_set0_2_mod_p{zz2,1}=(mask01_2.*maxEVs_num_set0_2-mask01_2_mod_p.*maxEVs_num_set0_2_mod_p);
diff_maxEVs_num_set0_2_mod_m{zz2,1}=(mask01_2.*maxEVs_num_set0_2-mask01_2_mod_m.*maxEVs_num_set0_2_mod_m);

gain_num_set0_2=mask01_2.*max(gain_num_cell{zz2,1},0);

gain_num_set0_2_mod_p=mask01_2_mod_p.*max(gain_num_mod_p_cell{zz2,1},0);
diff_gain_num_set0_2_mod_p{zz2,1}=gain_num_set0_2_mod_p-gain_num_set0_2;

gain_num_set0_2_mod_m=mask01_2_mod_m.*max(gain_num_mod_m_cell{zz2,1},0);
diff_gain_num_set0_2_mod_m{zz2,1}=gain_num_set0_2_mod_m-gain_num_set0_2;

nr_points=100;
subset_x=randi(length(diff_maxEVs_num_set0_2_mod_p{zz2,1}(:,1)),nr_points,1); % plot random subset of points
subset_y=randi(length(diff_maxEVs_num_set0_2_mod_p{zz2,1}(1,:)),nr_points,1);

% find changes in same direction
dummy_diff_gain_p=diff_gain_num_set0_2_mod_p{zz2,1};
dummy_diff_gain_p(dummy_diff_gain_p<0)=-1;
dummy_diff_gain_p(dummy_diff_gain_p>0)=1;

dummy_diff_EV_p=diff_maxEVs_num_set0_2_mod_p{zz2,1};
dummy_diff_EV_p(dummy_diff_EV_p<0)=-1;
dummy_diff_EV_p(dummy_diff_EV_p>0)=1;

% calculate percentag of datapoints in each quadrant
thresh_gain=0.1; % threshold to qualify as gain change
thresh_stab=0.01; % threshold to qualify as stability change

Q1_p=sum(sum((diff_maxEVs_num_set0_2_mod_p{zz2,1}<-thresh_stab).*(diff_gain_num_set0_2_mod_p{zz2,1}>thresh_gain)));
Q2_p=sum(sum((diff_maxEVs_num_set0_2_mod_p{zz2,1}>thresh_stab).*(diff_gain_num_set0_2_mod_p{zz2,1}>thresh_gain)));
Q3_p=sum(sum((diff_maxEVs_num_set0_2_mod_p{zz2,1}<-thresh_stab).*(diff_gain_num_set0_2_mod_p{zz2,1}<-thresh_gain)));
Q4_p=sum(sum((diff_maxEVs_num_set0_2_mod_p{zz2,1}>thresh_stab).*(diff_gain_num_set0_2_mod_p{zz2,1}<-thresh_gain)));


Q1_m=sum(sum((diff_maxEVs_num_set0_2_mod_m{zz2,1}<-thresh_stab).*(diff_gain_num_set0_2_mod_m{zz2,1}>thresh_gain)));
Q2_m=sum(sum((diff_maxEVs_num_set0_2_mod_m{zz2,1}>thresh_stab).*(diff_gain_num_set0_2_mod_m{zz2,1}>thresh_gain)));
Q3_m=sum(sum((diff_maxEVs_num_set0_2_mod_m{zz2,1}<-thresh_stab).*(diff_gain_num_set0_2_mod_m{zz2,1}<-thresh_gain)));
Q4_m=sum(sum((diff_maxEVs_num_set0_2_mod_m{zz2,1}>thresh_stab).*(diff_gain_num_set0_2_mod_m{zz2,1}<-thresh_gain)));

total_mod=Q1_p+Q2_p+Q3_p+Q4_p+Q1_m+Q2_m+Q3_m+Q4_m;

Q1_perc=(Q1_p+Q1_m)/total_mod*100;
Q2_perc=(Q2_p+Q2_m)/total_mod*100;
Q3_perc=(Q3_p+Q3_m)/total_mod*100;
Q4_perc=(Q4_p+Q4_m)/total_mod*100;

hold on

scatter(rS0_vec(zz2),Q1_perc,10,'k','filled')
scatter(rS0_vec(zz2),Q2_perc,10,'r','filled')
scatter(rS0_vec(zz2),Q3_perc,10,'g','filled')
scatter(rS0_vec(zz2),Q4_perc,10,'b','filled')

hold off
xlabel('r_S')
ylabel('Q1-Q4 (%)')
axis square
end

fig2.Renderer='Painters';
%saveas(gcf,['suppl_fig_rS_1.pdf'])