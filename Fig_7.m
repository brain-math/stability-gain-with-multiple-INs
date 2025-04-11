%% Code to reproduce Fig. 7 & Fig. 7 - figure supplement 1 of Bos, Miehl et al., 2025

%% analytically calculate the respective M terms, gain and det
syms bEx bPx bSx wEEx wEPx wESx wPEx wPPx wPSx wSEx wSPx wSSx dEx dPx dSx

vec_syms=[bEx, bPx, bSx, wEEx, wEPx, wESx, wPEx, wPPx, wPSx, wSEx, wSPx, wSSx];
W_sym=[[wEEx, -wEPx, -wESx];[wPEx, -wPPx,-wPSx];[wSEx,-wSPx,-wSSx]];
D_sym=[[1/bEx,0,0];[0,1/bPx,0];[0,0,1/bSx]];
D_sym2=[[dEx,0,0];[0,dPx,0];[0,0,dSx]];

L_sym=inv(D_sym-W_sym);
det_sym=det(D_sym-W_sym);

L_sym2=inv(D_sym2-W_sym);
det_sym2=det(D_sym2-W_sym);

J_sym=-D_sym+W_sym;


%% choose weights & rates
power=2; 
mult_f=1/4;
I_stim_E=0.3; I_stim_P=0.3;

wSS=0; %SOM-to-SOM weight

rE=3; % E rate
rP=5; % PV rate
rS=0.5; % SOM rate

bE=mult_f*power*(rE/mult_f)^((power-1)/power);
bP=mult_f*power*(rP/mult_f)^((power-1)/power);
bS=mult_f*power*(rS/mult_f)^((power-1)/power);

wEE_vec=[0.5, 0.8]; % non-ISN (1-wEE*bE>0) & ISN regime (1-wEE*bE<0) [Fig. 6 is non-ISN and Fig. S3 is ISN regime]

wES_vec=[0.5,0.1]; % inhibitory (wES>wPS) versus disinhibitory (wES<wPS) regimes
wPS_vec=[0.1,0.5];

w_vec=[0:0.05:1];

for gg=1:length(wEE_vec)
    wEE=wEE_vec(gg);
    for gg3=1:length(wES_vec)
        
        for gg2=1:length(w_vec)
            % gain and stability for different wPP
            wEP=0.5; wPE=0.5; wPP=w_vec(gg2); % EI terms
            wES=wES_vec(gg3); wPS=wPS_vec(gg3); % FF terms
            wSE=0.5; wSP=0.5; % FB terms

            vec_nums=[bE,bP,bS,wEE,wEP,wES,wPE,wPP,wPS,wSE,wSP,wSS];

            gain_wPP(gg,gg2,gg3)=double(subs(L_sym(1,1),vec_syms,vec_nums))*I_stim_E+double(subs(L_sym(1,2),vec_syms,vec_nums))*I_stim_P;
            maxEVs_wPP(gg,gg2,gg3)=max(real(eig(double(subs(J_sym,vec_syms,vec_nums)))));

            % gain and stability for different wPE
            wEP=0.5; wPE=w_vec(gg2); wPP=0.5; % EI terms
            wES=wES_vec(gg3); wPS=wPS_vec(gg3); % FF terms
            wSE=0.5; wSP=0.5; % FB terms

            vec_nums=[bE,bP,bS,wEE,wEP,wES,wPE,wPP,wPS,wSE,wSP,wSS];

            gain_wPE(gg,gg2,gg3)=double(subs(L_sym(1,1),vec_syms,vec_nums))*I_stim_E+double(subs(L_sym(1,2),vec_syms,vec_nums))*I_stim_P;
            maxEVs_wPE(gg,gg2,gg3)=max(real(eig(double(subs(J_sym,vec_syms,vec_nums)))));

            % gain and stability for different wEP
            wEP=w_vec(gg2); wPE=0.5; wPP=0.5; % EI terms
            wES=wES_vec(gg3); wPS=wPS_vec(gg3); % FF terms
            wSE=0.5; wSP=0.5; % FB terms

            vec_nums=[bE,bP,bS,wEE,wEP,wES,wPE,wPP,wPS,wSE,wSP,wSS];

            gain_wEP(gg,gg2,gg3)=double(subs(L_sym(1,1),vec_syms,vec_nums))*I_stim_E+double(subs(L_sym(1,2),vec_syms,vec_nums))*I_stim_P;
            maxEVs_wEP(gg,gg2,gg3)=max(real(eig(double(subs(J_sym,vec_syms,vec_nums)))));

            % gain and stability for different wSE
            wEP=0.5; wPE=0.5; wPP=0.5; % EI terms
            wES=w_vec(gg2); wPS=wPS_vec(gg3); % FF terms
            wSE=0.5; wSP=0.5; % FB terms

            vec_nums=[bE,bP,bS,wEE,wEP,wES,wPE,wPP,wPS,wSE,wSP,wSS];

            gain_wES(gg,gg2,gg3)=double(subs(L_sym(1,1),vec_syms,vec_nums))*I_stim_E+double(subs(L_sym(1,2),vec_syms,vec_nums))*I_stim_P;
            maxEVs_wES(gg,gg2,gg3)=max(real(eig(double(subs(J_sym,vec_syms,vec_nums)))));

            % gain and stability for different wPS
            wEP=0.5; wPE=0.5; wPP=0.5; % EI terms
            wES=wES_vec(gg3); wPS=w_vec(gg2); % FF terms
            wSE=0.5; wSP=0.5; % FB terms

            vec_nums=[bE,bP,bS,wEE,wEP,wES,wPE,wPP,wPS,wSE,wSP,wSS];

            gain_wPS(gg,gg2,gg3)=double(subs(L_sym(1,1),vec_syms,vec_nums))*I_stim_E+double(subs(L_sym(1,2),vec_syms,vec_nums))*I_stim_P;
            maxEVs_wPS(gg,gg2,gg3)=max(real(eig(double(subs(J_sym,vec_syms,vec_nums)))));

            % gain and stability for different wSE
            wEP=0.5; wPE=0.5; wPP=0.5; % EI terms
            wES=wES_vec(gg3); wPS=wPS_vec(gg3); % FF terms
            wSE=w_vec(gg2); wSP=0.5; % FB terms

            vec_nums=[bE,bP,bS,wEE,wEP,wES,wPE,wPP,wPS,wSE,wSP,wSS];

            gain_wSE(gg,gg2,gg3)=double(subs(L_sym(1,1),vec_syms,vec_nums))*I_stim_E+double(subs(L_sym(1,2),vec_syms,vec_nums))*I_stim_P;
            maxEVs_wSE(gg,gg2,gg3)=max(real(eig(double(subs(J_sym,vec_syms,vec_nums)))));

            % gain and stability for different wSP
            wEP=0.5; wPE=0.5; wPP=0.5; % EI terms
            wES=wES_vec(gg3); wPS=wPS_vec(gg3); % FF terms
            wSE=0.5; wSP=w_vec(gg2); % FB terms

            vec_nums=[bE,bP,bS,wEE,wEP,wES,wPE,wPP,wPS,wSE,wSP,wSS];

            gain_wSP(gg,gg2,gg3)=double(subs(L_sym(1,1),vec_syms,vec_nums))*I_stim_E+double(subs(L_sym(1,2),vec_syms,vec_nums))*I_stim_P;
            maxEVs_wSP(gg,gg2,gg3)=max(real(eig(double(subs(J_sym,vec_syms,vec_nums)))));
            
            % gain and stability for different wSS
            wEP=0.5; wPE=0.5; wPP=0.5; % EI terms
            wES=wES_vec(gg3); wPS=wPS_vec(gg3); % FF terms
            wSE=0.5; wSP=0.5; % FB terms
            wSS=w_vec(gg2); % SST FB

            vec_nums=[bE,bP,bS,wEE,wEP,wES,wPE,wPP,wPS,wSE,wSP,wSS];

            gain_wSS(gg,gg2,gg3)=double(subs(L_sym(1,1),vec_syms,vec_nums))*I_stim_E+double(subs(L_sym(1,2),vec_syms,vec_nums))*I_stim_P;
            maxEVs_wSS(gg,gg2,gg3)=max(real(eig(double(subs(J_sym,vec_syms,vec_nums)))));
            

        end
    end
end



%% figure

vec=120:-20:0;
NNN = 128;
hex_red=['#fee5d9','#fcbba1','#fc9272','#fb6a4a','#ef3b2c','#cb181d','#99000d']'; %red
raw_red = sscanf(hex_red','#%2x%2x%2x',[3,size(hex_red,1)]).' / 255;
map_red = interp1(vec,raw_red,linspace(120,0,NNN),'pchip');

hex_blue=['#eff3ff','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#084594']'; %blue
raw_blue = sscanf(hex_blue','#%2x%2x%2x',[3,size(hex_blue,1)]).' / 255;
map_blue = interp1(vec,raw_blue,linspace(120,0,NNN),'pchip');

hex_green=['#edf8e9','#c7e9c0','#a1d99b','#74c476','#41ab5d','#238b45','#005a32']'; %green
raw_green = sscanf(hex_green','#%2x%2x%2x',[3,size(hex_green,1)]).' / 255;
map_green = interp1(vec,raw_green,linspace(120,0,NNN),'pchip');

ylim_nonisn=[-0.5 1];

fig=figure
subplot(2,2,1)
hold on 
patch([(-1)*maxEVs_wPS(1,:,1) nan],[gain_wPS(1,:,1) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_blue);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wPS - Inh & non-ISN')
ylim(ylim_nonisn)
xlim([0 0.9])

subplot(2,2,2)
hold on 
%xline(0,'k--','Linewidth',2)
patch([(-1)*maxEVs_wPS(2,:,1) nan],[gain_wPS(2,:,1) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_blue);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wPS - Inh & ISN')
xlim([0 0.45])
ylim([-5 10])

subplot(2,2,3)
hold on 
%xline(0,'k--','Linewidth',2)
patch([(-1)*maxEVs_wPS(1,:,2) nan],[gain_wPS(1,:,2) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_blue);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wPS - Disinh & non-ISN')
ylim(ylim_nonisn)
xlim([0 0.9])

subplot(2,2,4)
hold on 
%xline(0,'k--','Linewidth',2)
patch([(-1)*maxEVs_wPS(2,:,2) nan],[gain_wPS(2,:,2) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_blue);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wPS - Disinh & ISN')
xlim([0 0.45])
ylim([-5 10])

fig.Renderer='Painters';

%saveas(gcf,'panel_fig_wPS.pdf')

fig2=figure
subplot(2,2,1)
hold on 
patch([(-1)*maxEVs_wES(1,:,1) nan],[gain_wES(1,:,1) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_red);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wES - Inh & non-ISN')
ylim(ylim_nonisn)
xlim([0 0.9])

subplot(2,2,2)
hold on 
patch([(-1)*maxEVs_wES(2,:,1) nan],[gain_wES(2,:,1) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_red);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wES - Inh & ISN')
xlim([0 0.45])
ylim([-5 10])

subplot(2,2,3)
hold on 
patch([(-1)*maxEVs_wES(1,:,2) nan],[gain_wES(1,:,2) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_red);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wES - Disinh & non-ISN')
ylim(ylim_nonisn)
xlim([0 0.9])

subplot(2,2,4)
hold on 
patch([(-1)*maxEVs_wES(2,:,2) nan],[gain_wES(2,:,2) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_red);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wES - Disinh & ISN')
xlim([0 0.45])
ylim([-5 10])

fig2.Renderer='Painters';

%saveas(gcf,'panel_fig_wES.pdf')

fig3=figure
subplot(2,2,1)
hold on 
patch([(-1)*maxEVs_wSE(1,:,1) nan],[gain_wSE(1,:,1) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_green);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wSE - Inh & non-ISN')
ylim(ylim_nonisn)
xlim([0 0.9])

subplot(2,2,2)
hold on 
patch([(-1)*maxEVs_wSE(2,:,1) nan],[gain_wSE(2,:,1) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_green);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wSE - Inh & ISN')
xlim([0 0.45])
ylim([-5 10])

subplot(2,2,3)
hold on 
patch([(-1)*maxEVs_wSE(1,:,2) nan],[gain_wSE(1,:,2) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_green);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wSE - Disinh & non-ISN')
ylim(ylim_nonisn)
xlim([0 0.9])

subplot(2,2,4)
hold on 
patch([(-1)*maxEVs_wSE(2,:,2) nan],[gain_wSE(2,:,2) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_green);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wSE - Disinh & ISN')
xlim([0 0.45])
ylim([-5 10])

fig3.Renderer='Painters';

%saveas(gcf,'panel_fig_wSE.pdf')

fig4=figure
subplot(2,2,1)
hold on 
patch([(-1)*maxEVs_wSP(1,:,1) nan],[gain_wSP(1,:,1) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_green);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wSP - Inh & non-ISN')
ylim(ylim_nonisn)
xlim([0 0.9])

subplot(2,2,2)
hold on 
patch([(-1)*maxEVs_wSP(2,:,1) nan],[gain_wSP(2,:,1) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_green);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wSP - Inh & ISN')
xlim([0 0.45])
ylim([-5 10])

subplot(2,2,3)
hold on
patch([(-1)*maxEVs_wSP(1,:,2) nan],[gain_wSP(1,:,2) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_green);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wSP - Disinh & non-ISN')
ylim(ylim_nonisn)
xlim([0 0.9])

subplot(2,2,4)
hold on 
patch([(-1)*maxEVs_wSP(2,:,2) nan],[gain_wSP(2,:,2) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_green);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wSP - Disinh & ISN')
xlim([0 0.45])
ylim([-5 10])

fig4.Renderer='Painters';

%saveas(gcf,'panel_fig_wSP.pdf')

fig5=figure
subplot(2,2,1)
hold on 
patch([(-1)*maxEVs_wPP(1,:,1) nan],[gain_wPP(1,:,1) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_blue);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wPP - Inh & non-ISN')
ylim(ylim_nonisn)
xlim([0 0.9])

subplot(2,2,2)
hold on 
patch([(-1)*maxEVs_wPP(2,:,1) nan],[gain_wPP(2,:,1) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_blue);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wPP - Inh & ISN')
xlim([0 0.45])
ylim([-5 10])

subplot(2,2,3)
hold on 
patch([(-1)*maxEVs_wPP(1,:,2) nan],[gain_wPP(1,:,2) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_blue);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wPP - Disinh & non-ISN')
ylim(ylim_nonisn)
xlim([0 0.9])

subplot(2,2,4)
hold on 
patch([(-1)*maxEVs_wPP(2,:,2) nan],[gain_wPP(2,:,2) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_blue);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wPP - Disinh & ISN')
xlim([0 0.45])
ylim([-5 10])

fig5.Renderer='Painters';

%saveas(gcf,'panel_fig_wPP.pdf')

fig6=figure
subplot(2,2,1)
hold on 
patch([(-1)*maxEVs_wPE(1,:,1) nan],[gain_wPE(1,:,1) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_blue);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wPE - Inh & non-ISN')
ylim(ylim_nonisn)
xlim([0 0.9])

subplot(2,2,2)
hold on 
patch([(-1)*maxEVs_wPE(2,:,1) nan],[gain_wPE(2,:,1) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_blue);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wPE - Inh & ISN')
xlim([0 0.45])
ylim([-5 10])

subplot(2,2,3)
hold on 
patch([(-1)*maxEVs_wPE(1,:,2) nan],[gain_wPE(1,:,2) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_blue);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wPE - Disinh & non-ISN')
ylim(ylim_nonisn)
xlim([0 0.9])

subplot(2,2,4)
hold on 
patch([(-1)*maxEVs_wPE(2,:,2) nan],[gain_wPE(2,:,2) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_blue);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wPE - Disinh & ISN')
xlim([0 0.45])
ylim([-5 10])

fig6.Renderer='Painters';

%saveas(gcf,'panel_fig_wPE.pdf')

fig7=figure
subplot(2,2,1)
hold on 
patch([(-1)*maxEVs_wEP(1,:,1) nan],[gain_wEP(1,:,1) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_red);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wEP - Inh & non-ISN')
ylim(ylim_nonisn)
xlim([0 0.9])

subplot(2,2,2)
hold on 
patch([(-1)*maxEVs_wEP(2,:,1) nan],[gain_wEP(2,:,1) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_red);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wEP - Inh & ISN')
xlim([0 0.45])
ylim([-5 10])

subplot(2,2,3)
hold on 
patch([(-1)*maxEVs_wEP(1,:,2) nan],[gain_wEP(1,:,2) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_red);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wEP - Disinh & non-ISN')
ylim(ylim_nonisn)
xlim([0 0.9])

subplot(2,2,4)
hold on 
patch([(-1)*maxEVs_wEP(2,:,2) nan],[gain_wEP(2,:,2) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_red);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wEP - Disinh & ISN')
xlim([0 0.45])
ylim([-5 10])

fig7.Renderer='Painters';

%saveas(gcf,'panel_fig_wEP.pdf')

fig3=figure
subplot(2,2,1)
hold on 
patch([(-1)*maxEVs_wSS(1,:,1) nan],[gain_wSS(1,:,1) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_green);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wSS - Inh & non-ISN')
ylim(ylim_nonisn)
xlim([0 0.9])

subplot(2,2,2)
hold on 
patch([(-1)*maxEVs_wSS(2,:,1) nan],[gain_wSS(2,:,1) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_green);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wSS - Inh & ISN')
xlim([0 0.45])
ylim([-5 10])

subplot(2,2,3)
hold on 
patch([(-1)*maxEVs_wSS(1,:,2) nan],[gain_wSS(1,:,2) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_green);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wSS - Disinh & non-ISN')
ylim(ylim_nonisn)
xlim([0 0.9])

subplot(2,2,4)
hold on 
patch([(-1)*maxEVs_wSS(2,:,2) nan],[gain_wSS(2,:,2) nan],[w_vec nan],[w_vec nan], 'edgecolor', 'interp','Linewidth',2); 
c=colorbar;
colormap(map_green);
hold off
xlabel('Stability')
ylabel('Gain')
c.Label.String ='w';
title('wSS - Disinh & ISN')
xlim([0 0.45])
ylim([-5 10])

fig3.Renderer='Painters';

%saveas(gcf,'panel_fig_wSS.pdf')
