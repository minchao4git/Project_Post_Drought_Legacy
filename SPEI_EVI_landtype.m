%% For the whole summer period May-Aug by land cover type
sz=12;
m=13;
m_rng=[5:8];
yscal=0.95; xscal=0.95;


figure('color','white','Position',[330   422   860   552]);
% Normal years
r_ind=repmat((r_spei_evi(:,:,m)>=0.5),[1 1 1 20]);
spei_ex_m_y0=squeeze(nanmean(DATA_SPEI_out(:,:,m_rng,:),3)).*r_ind.*ex_y0;
evi_ex_m_y0=squeeze(nanmean(DATA_VI_05rs_dtds_out(:,:,m_rng,:,dsid),3)).*r_ind.*ex_y0;
a=spei_ex_m_y0(:);
b=evi_ex_m_y0(:);
nnan_idc=~isnan(a) & ~isnan(b) & (a~=0) & (b~=0);
hold on;
scatter(spei_ex_m_y0(:),evi_ex_m_y0(:),sz,'b');
x_rng2=[nanmin(spei_ex_m_y0(nnan_idc)) nanmax(spei_ex_m_y0(nnan_idc))];
md2=fitlm(spei_ex_m_y0(nnan_idc),evi_ex_m_y0(nnan_idc));
plot([-3 3],[0 0],'--k');
plot([0 0],[-3 3],'--k');

% === Extreme years
% figure;
r_ind=repmat((r_spei_evi(:,:,m)>=0.5),[1 1 1 20]);
spei_ex_m_y1=squeeze(nanmean(DATA_SPEI_out(:,:,m_rng,:),3)).*r_ind.*ex_y1;
evi_ex_m_y1=squeeze(nanmean(DATA_VI_05rs_dtds_out(:,:,m_rng,:,dsid),3)).*r_ind.*ex_y1;
a=spei_ex_m_y1(:);
b=evi_ex_m_y1(:);
nnan_idc=~isnan(a) & ~isnan(b) & (a~=0) & (b~=0);
hold on;
scatter(spei_ex_m_y1(nnan_idc),evi_ex_m_y1(nnan_idc),sz,'r');
x_rng1=[nanmin(spei_ex_m_y1(nnan_idc)) nanmax(spei_ex_m_y1(nnan_idc))];
md1=fitlm(spei_ex_m_y1(nnan_idc),evi_ex_m_y1(nnan_idc));
plot([-3 3],[0 0],'--k');
plot([0 0],[-3 3],'--k');

% Post-extreme years
r_ind=repmat((r_spei_evi(:,:,m)>=0.5),[1 1 1 20]);
spei_ex_m_pst=squeeze(nanmean(DATA_SPEI_out(:,:,m_rng,:),3)).*r_ind.*ex_pst;
evi_ex_m_pst=squeeze(nanmean(DATA_VI_05rs_dtds_out(:,:,m_rng,:,dsid),3)).*r_ind.*ex_pst;
a=spei_ex_m_pst(:);
b=evi_ex_m_pst(:);
nnan_idc=~isnan(a) & ~isnan(b) & (a~=0) & (b~=0);
hold on;
scatter(spei_ex_m_pst(:),evi_ex_m_pst(:),sz,'g');
x_rng3=[nanmin(spei_ex_m_pst(nnan_idc)) nanmax(spei_ex_m_pst(nnan_idc))];
md3=fitlm(spei_ex_m_pst(nnan_idc),evi_ex_m_pst(nnan_idc));
plot([-3 3],[0 0],'--k');
plot([0 0],[-3 3],'--k');
xlim([-2 2]);ylim([-0.15 0.1]);
title('SPEI vs. EVI for normal, extreme, and post-extreme years for May-Aug');

% Key parameters
s1=md1.Coefficients.SE(2); n1=md1.NumObservations; a1=md1.Coefficients.Estimate(1); cof1=md1.Coefficients.Estimate(2);
s2=md2.Coefficients.SE(2); n2=md2.NumObservations; a2=md2.Coefficients.Estimate(1); cof2=md2.Coefficients.Estimate(2);
s3=md3.Coefficients.SE(2); n3=md2.NumObservations; a3=md3.Coefficients.Estimate(1); cof3=md3.Coefficients.Estimate(2);

plot(x_rng2,a2+cof2*x_rng2,'color',[0 0 255]/255,'LineWidth',3,'LineStyle','--');
plot(x_rng1,a1+cof1*x_rng1,'color',[170 0 0]/255,'LineWidth',3,'LineStyle','--');
plot(x_rng3,a3+cof3*x_rng3,'color',[0 85 0]/255,'LineWidth',3,'LineStyle','--');

% Significance test
Snrm_ex=sqrt(power(s2,2)/n2+power(s1,2)/n1);
Snrm_pst=sqrt(power(s2,2)/n2+power(s3,2)/n3);

Tnrm_ex=(cof1-cof2)/Snrm_ex;
Tnrm_pst=(cof3-cof2)/Snrm_pst;

sig_nrm_ex=0; sig_nrm_pst=0;
if abs(Tnrm_ex) > 1.96
    sig_nrm_ex=1;
    fw_nrm_ex='bold';
    sig_sym_nrm_ex='*';
else
    fw_nrm_ex='normal';
    sig_sym_nrm_ex='';
end
if abs(Tnrm_pst) > 1.96
    sig_nrm_pst=1;
    fw_nrm_pst='bold';
    sig_sym_nrm_pst='*';
else
    fw_nrm_pst='normal';
    sig_sym_nrm_pst='';
end

text(-1.9,0.08,sprintf('EVI_N_o_r_m_a_l_Y=%3.3fxSPEI+%2.2f',cof2,a2),'color','b','FontSize',9);
text(-1.75,0.06,sprintf('EVI_E_x_Y=%3.3fxSPEI+%2.2f^%s',cof1,a1,sig_sym_nrm_ex),'color','r','FontSize',9, 'FontWeight',fw_nrm_ex);
text(-1.85,0.04,sprintf('EVI_P_s_t_E_x_Y=%3.3fxSPEI+%2.2f^%s',cof3,a3,sig_sym_nrm_pst),'color','g','FontSize',9, 'FontWeight',fw_nrm_pst);
a=get(gca); gca_w = (a.XLim(2)-a.XLim(1)); gca_h = (a.YLim(2)-a.YLim(1));
text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,'a)','FontSize',11,'FontName','Arial');
box on;

xlabel('SPEI anomalies'); ylabel('EVI anomalies');