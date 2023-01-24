% This script is to explain the bell-shape of the recovery curve (drought occuring in
% EGS)
subgs=1;

[ex_m ex_pst_m norm_m eem_map ex_pst_m_gs]=def_extreme_m(squeeze(DATA_05rs_dtds_out_gs(:,:,:,:,3)),1, 'SAT', subgs);
% identify the extremes occured in 1st or 2nd GS
ex_m_gs=find_ex_m_gs(ex_m, ex_pst_m_gs);

a=lgs_map(:,:,2);
a=(a==9);
b=squeeze(nanmax(nanmax(ex_pst_m_gs.*(ex_m_gs==1),[],3),[],4));
b=(b>6);

c=lc_dom_grp;
c=(c==6);

out=a.*b.*c;

[s1 s2]=size(out);

n=0;
for i=1:s1
   for j=1:s2
       if out(i,j)==1
           n=n+1;
           x_list(n)=i;
           y_list(n)=j;
       end
   end
end
% x=173; y=29;

d=reshape(dtmp_dtds(:,:,:,:,2),[ 720   102 12*20]);
e=reshape(DATA_05rs_dtds_out_gs(:,:,:,:,2),[ 720   102 12*19]);

for i=[2, 3, 4, 5, 9, 10, 11, 12, 13, 15, 19]
    x=x_list(i); y=y_list(i);
    lc_dom_grp(x,y)

    ex_pst_ts=reshape(ex_pst_m(x,y,:,:),[1 12*19]);
    ex_pst_gs_ts=reshape(ex_pst_m_gs(x,y,:,:),[1 12*19]);
    ex_ts=reshape(squeeze(ex_m(x,y,:,:)),[1 12*19]);
    ex_m_gs_ts=reshape(squeeze(ex_m_gs(x,y,:,:)),[1 12*19]);

    spei_ts=reshape(squeeze(DATA_05rs_dtds_out_gs(x,y,:,:,3)),[1 12*19]);
    % spei_ts_all=reshape(squeeze(DATA_SPEI_out(x,y,:,:)),[1 12*40]);
    nrm_ts=reshape(norm_m(x,y,:,:),[1 12*19]);

    clm_vi=nanmean((squeeze(DATA_VI_05rs_sm_out(x,y,:,:,2,3))),2);
    clm_vi_ts=reshape(repmat(clm_vi,[1 20]),[1 12*20]);
    ats=reshape(squeeze(DATA_VI_05rs_sm_out(x,y,:,:,2,3)),[1 12*20]);


    figure('Position',[186         299        1410         540]);
    subplot(3,1,1);
    hold on;
    plot(clm_vi_ts(3:230),'--k','Linewidth',2);
    plot(ats(3:230),'-m','Linewidth',2);
    title(sprintf('i: %d, x: %d, y: %d',i, x,y));

    subplot(3,1,2);
    hold on;
    plot([1 12*20],[0 0],'--k','Linewidth',2);
    plot(squeeze(e(x,y,:)),'-b','Linewidth',2);
    plot(squeeze(ex_m_gs_ts)/10,'-m','Linewidth',2);
    plot(ex_pst_ts/50, 'Linewidth',2);

    subplot(3,1,3);
    plot(ex_ts,'-b','Linewidth',1.5);

end

% Example list: 
plot_list=[2, 3, 4, 5, 9, 10, 11, 12, 13, 15, 19];

%%
for i=[4]
    x=x_list(i); y=y_list(i);
    lc_dom_grp(x,y)

    ex_pst_ts=reshape(ex_pst_m(x,y,:,:),[1 12*19]);
    ex_pst_gs_ts=reshape(ex_pst_m_gs(x,y,:,:),[1 12*19]);
    ex_ts=reshape(squeeze(ex_m(x,y,:,:)),[1 12*19]);
    ex_m_gs_ts=reshape(squeeze(ex_m_gs(x,y,:,:)),[1 12*19]);

    spei_ts=reshape(squeeze(DATA_05rs_dtds_out_gs(x,y,:,:,3)),[1 12*19]);
    % spei_ts_all=reshape(squeeze(DATA_SPEI_out(x,y,:,:)),[1 12*40]);
    nrm_ts=reshape(norm_m(x,y,:,:),[1 12*19]);

    clm_vi=nanmean((squeeze(DATA_VI_05rs_sm_out(x,y,:,:,2,3))),2);
    clm_vi_ts=reshape(repmat(clm_vi,[1 20]),[1 12*20]);
    ats=reshape(squeeze(DATA_VI_05rs_sm_out(x,y,:,:,2,3)),[1 12*20]);

    ex_list=find(ex_ts==1);
    
    figure('color','w','Position',[181         416        1246         418]);
    subplot(2,1,1);
    hold on;
    for n=1:length(ex_list)
        plot([ex_list(n) ex_list(n)], [-100 100],'--r','Linewidth',1.5);
    end
    plot(clm_vi_ts(3:230),'--k','Linewidth',2);
    plot(ats(3:230),'-m','Linewidth',2);
    title('Example of a savanna site in NH, seasonal variability of VI, drought occures in EGS');
    ylim([0.25 0.65]);
    xlim([73 167]);
    box on;

    subplot(2,1,2);
    hold on;
    for n=1:length(ex_list)
        plot([ex_list(n) ex_list(n)], [-100 100],'--r','Linewidth',1.5);
    end
    plot([1 12*20],[0 0],'--k','Linewidth',2);
    plot(squeeze(e(x,y,:)),'-b','Linewidth',2);
    gs1=squeeze(ex_m_gs_ts==1)/4; gs1(gs1==0)=nan;
    gs2=squeeze(ex_m_gs_ts==2)/4; gs2(gs2==0)=nan;
    plot(gs1,'-m','Linewidth',3);
    plot(gs2,'-b','Linewidth',3);
%     plot(ex_pst_ts/50, 'Linewidth',2);
    title('VI anomalies');
    ylim([-0.25 0.3]);
    xlim([73 167]);
    xlabel('months');
    box on;
end
% export_fig(('Example_savanna_drought_EGS.png'), '-png', '-r600');
% sos_map(x,y,2)
% eos_map(x,y,2)

%%
% This script is to explain positive anomalies after drought occuring in
% LGS
subgs=2;

[ex_m ex_pst_m norm_m eem_map ex_pst_m_gs]=def_extreme_m(squeeze(DATA_05rs_dtds_out_gs(:,:,:,:,3)),1, 'SAT', subgs);
% identify the extremes occured in 1st or 2nd GS
ex_m_gs=find_ex_m_gs(ex_m, ex_pst_m_gs);

a=lgs_map(:,:,2);
a=(a==9);
b=squeeze(nanmax(nanmax(ex_pst_m_gs.*(ex_m_gs==1),[],3),[],4));
b=(b>2);

c=lc_dom_grp;
c=(c==6);

out=a.*b.*c;

figure;
imagesc((out));


%%
[s1 s2]=size(out);

n=0;
for i=1:s1
   for j=1:s2
       if out(i,j)==1
           n=n+1;
           x_list(n)=i;
           y_list(n)=j;
       end
   end
end

d=reshape(dtmp_dtds(:,:,:,:,2),[ 720   102 12*20]);
e=reshape(DATA_05rs_dtds_out_gs(:,:,:,:,2),[ 720   102 12*19]);

for i=[2, 10, 11, 12, 13, 15, 17, 37]
    x=x_list(i); y=y_list(i);
    lc_dom_grp(x,y)

    ex_pst_ts=reshape(ex_pst_m(x,y,:,:),[1 12*19]);
    ex_pst_gs_ts=reshape(ex_pst_m_gs(x,y,:,:),[1 12*19]);
    ex_ts=reshape(squeeze(ex_m(x,y,:,:)),[1 12*19]);
    ex_m_gs_ts=reshape(squeeze(ex_m_gs(x,y,:,:)),[1 12*19]);

    spei_ts=reshape(squeeze(DATA_05rs_dtds_out_gs(x,y,:,:,3)),[1 12*19]);
    nrm_ts=reshape(norm_m(x,y,:,:),[1 12*19]);

    clm_vi=nanmean((squeeze(DATA_VI_05rs_sm_out(x,y,:,:,2,3))),2);
    clm_vi_ts=reshape(repmat(clm_vi,[1 20]),[1 12*20]);
    ats=reshape(squeeze(DATA_VI_05rs_sm_out(x,y,:,:,2,3)),[1 12*20]);


    figure('Position',[186         299        1410         540]);
    subplot(3,1,1);
    hold on;
    plot(clm_vi_ts(3:230),'--k','Linewidth',2);
    plot(ats(3:230),'-m','Linewidth',2);
    title(sprintf('i: %d, x: %d, y: %d',i, x,y));

    subplot(3,1,2);
    hold on;
    plot([1 12*20],[0 0],'--k','Linewidth',2);
    plot(squeeze(e(x,y,:)),'-b','Linewidth',2);
    plot(squeeze(ex_m_gs_ts)/10,'-m','Linewidth',2);
    plot(ex_pst_ts/50, 'Linewidth',2);

    subplot(3,1,3);
    plot(ex_ts,'-b','Linewidth',1.5);

end

% Example list: 
plot_list=[2, 10, 11, 12, 13, 15, 17, 37];

%%
for i=[12]
    x=x_list(i); y=y_list(i);
    lc_dom_grp(x,y)

    ex_pst_ts=reshape(ex_pst_m(x,y,:,:),[1 12*19]);
    ex_pst_gs_ts=reshape(ex_pst_m_gs(x,y,:,:),[1 12*19]);
    ex_ts=reshape(squeeze(ex_m(x,y,:,:)),[1 12*19]);
    ex_m_gs_ts=reshape(squeeze(ex_m_gs(x,y,:,:)),[1 12*19]);

    spei_ts=reshape(squeeze(DATA_05rs_dtds_out_gs(x,y,:,:,3)),[1 12*19]);
    nrm_ts=reshape(norm_m(x,y,:,:),[1 12*19]);

    clm_vi=nanmean((squeeze(DATA_VI_05rs_sm_out(x,y,:,:,2,3))),2);
    clm_vi_ts=reshape(repmat(clm_vi,[1 20]),[1 12*20]);
    ats=reshape(squeeze(DATA_VI_05rs_sm_out(x,y,:,:,2,3)),[1 12*20]);

    ex_list=find(ex_ts==1);
    
    figure('color','w','Position',[181         416        1246         418]);
    subplot(2,1,1);
    hold on;
    for n=1:length(ex_list)
        plot([ex_list(n) ex_list(n)], [-100 100],'--r','Linewidth',1.5);
    end
    plot(clm_vi_ts(3:230),'--k','Linewidth',2);
    plot(ats(3:230),'-m','Linewidth',2);
    title('Example of a savanna site in NH, seasonal variability of VI, drought occures in LGS');
    ylim([0.15 0.65]);
    xlim([109 205]);
    box on;

    subplot(2,1,2);
    hold on;
    for n=1:length(ex_list)
        plot([ex_list(n) ex_list(n)], [-100 100],'--r','Linewidth',1.5);
    end
    plot([1 12*20],[0 0],'--k','Linewidth',2);
    plot(squeeze(e(x,y,:)),'-b','Linewidth',2);
    gs1=squeeze(ex_m_gs_ts==1)/4; gs1(gs1==0)=nan;
    gs2=squeeze(ex_m_gs_ts==2)/4; gs2(gs2==0)=nan;
    plot(gs1,'-m','Linewidth',3);
    plot(gs2,'-b','Linewidth',3);
%     plot(ex_pst_ts/50, 'Linewidth',2);
    title('VI anomalies');
    ylim([-0.3 0.3]);
    xlim([109 205]);
    xlabel('months');
    box on;
end
