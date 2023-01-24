function VI_SPEI_r_plot(etype_map)
    %% Fig. 2. GS-based VI-SPEI Correlation coefficient
	global ax_show lons lats;
    global r_spei_evi p_spei_evi;
    
    dtmp=r_spei_evi;
    dtmp(isnan(dtmp))=0;
    
    % Significant points
    [s1 s2 s3]=size(p_spei_evi);
    lons_sig=repmat(lons,[1 1 s3]);
    lats_sig=repmat(lats,[1 1 s3]);
    lons_sig(isnan(p_spei_evi))=nan;
    lats_sig(isnan(p_spei_evi))=nan;
    lons_sig(p_spei_evi>=0.05)=nan;
    lats_sig(p_spei_evi>=0.05)=nan;
    lons_sig(p_spei_evi<=0)=nan;
    lats_sig(p_spei_evi<=0)=nan;

    etype_EGS=etype_map{1,1};
    etype_LGS=etype_map{2,1};
    % negative effects from EGS and LGS
    etype_EGS_neg=squeeze(etype_EGS(:,:,1));
    etype_LGS_neg=squeeze(etype_LGS(:,:,1));
    etype_EGS_neg(isnan(etype_EGS_neg))=0;
    etype_LGS_neg(isnan(etype_LGS_neg))=0;
    
    % postive effects from EGS and LGS
    etype_EGS_pstv=squeeze(etype_EGS(:,:,2));
    etype_LGS_pstv=squeeze(etype_LGS(:,:,2));
    etype_EGS_pstv(isnan(etype_EGS_pstv))=0;
    etype_LGS_pstv(isnan(etype_LGS_pstv))=0;
    
    figure('color','w','Position',[38  118   1600   878]);
    gap_h=0.03; gap_w=0.002;
    gap=[gap_h gap_w]; marg_h=[0.08 0.05]; marg_w=[0.1 0.1];
    ha = tight_subplot(3,1,gap,marg_h,marg_w);

    yscal=0.95; xscal=0.02;
    cb_xR=1.0;
    cb_yR=0.5;
    cb_wR=0.4;
    cb_tR=2.5;
    
    color1=[0 0 0]/255; % All the significant points
    color_neg=[100 0 0]/255; % for the existance of negative points
    mrk_neg='s';
    lw=0.6;  % linewidth
    lw_neg=0.4;
    mks=4;
    mks2=4;
    axes(ha(1));
    hold on;
    m=15;
    bg_show=squeeze(dtmp(:,:,m))';
    geoplot(ax_show, bg_show);
    
    % Mark significant points for the entire GS
    a=squeeze(lats_sig(1:2:end,1:2:end,m));
    b=squeeze(lons_sig(1:2:end,1:2:end,m));
    plotm(a(:), b(:), 'Color',color1,'LineStyle','none', 'Marker', '+', 'MarkerSize',mks2, 'LineWidth',lw);
    
    a=squeeze(lats_sig(1:2:end,1:2:end,m)).*(etype_EGS_neg(1:2:end,1:2:end)|etype_LGS_neg(1:2:end,1:2:end));
    b=squeeze(lons_sig(1:2:end,1:2:end,m)).*(etype_EGS_neg(1:2:end,1:2:end)|etype_LGS_neg(1:2:end,1:2:end));
    plotm(a(:), b(:), 'Color',color_neg,'LineStyle','none', 'Marker', mrk_neg, 'MarkerSize',mks, 'LineWidth',lw_neg);
    
    title('VI~SPEI (r), Entire GS');
    cb1=colorbar('Eastoutside');
    crng=[-0.8 0.8];
    colormap(gca, flipud(cbrewer('div','RdBu',33)));
    caxis([crng(1) crng(2)]);
    a=get(gca); gca_w = (a.XLim(2)-a.XLim(1)); gca_h = (a.YLim(2)-a.YLim(1));
    text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,'a)','FontSize',11,'FontName','Arial');
    resizeCB(cb1, cb_xR, cb_yR, cb_wR, cb_tR, '',nan,nan,9);

    axes(ha(2));
    hold on;
    m=13;
    bg_show=squeeze(dtmp(:,:,m))';
    geoplot(ax_show, bg_show);
    
    % Mark significant points for EGS
    a=squeeze(lats_sig(1:2:end,1:2:end,m));
    b=squeeze(lons_sig(1:2:end,1:2:end,m));
    plotm(a(:), b(:), 'Color',color1,'LineStyle','none', 'Marker', '+', 'MarkerSize',mks2, 'LineWidth',lw);
    
    a=squeeze(lats_sig(1:2:end,1:2:end,m)).*(etype_EGS_neg(1:2:end,1:2:end));
    b=squeeze(lons_sig(1:2:end,1:2:end,m)).*(etype_EGS_neg(1:2:end,1:2:end));
    plotm(a(:), b(:), 'Color',color_neg,'LineStyle','none', 'Marker', mrk_neg, 'MarkerSize',mks, 'LineWidth',lw_neg);
    
    title('VI~SPEI (r), EGS');
    cb1=colorbar('Eastoutside');
    crng=[-0.8 0.8];
    colormap(gca, flipud(cbrewer('div','RdBu',33)));
    caxis([crng(1) crng(2)]);
    a=get(gca); gca_w = (a.XLim(2)-a.XLim(1)); gca_h = (a.YLim(2)-a.YLim(1));
    text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,'b)','FontSize',11,'FontName','Arial');
    resizeCB(cb1, cb_xR, cb_yR, cb_wR, cb_tR, '',nan,nan,9);

    axes(ha(3));
    hold on;
    m=14;
    bg_show=squeeze(dtmp(:,:,m))';
    geoplot(ax_show, bg_show);
    
    % Mark significant points for LGS
    a=squeeze(lats_sig(1:2:end,1:2:end,m));
    b=squeeze(lons_sig(1:2:end,1:2:end,m));
    plotm(a(:), b(:), 'Color',color1,'LineStyle','none', 'Marker', '+', 'MarkerSize',mks2, 'LineWidth',lw);
    a=squeeze(lats_sig(1:2:end,1:2:end,m)).*(etype_LGS_neg(1:2:end,1:2:end));
    b=squeeze(lons_sig(1:2:end,1:2:end,m)).*(etype_LGS_neg(1:2:end,1:2:end));
    plotm(a(:), b(:), 'Color',color_neg,'LineStyle','none', 'Marker', mrk_neg, 'MarkerSize',mks, 'LineWidth',lw_neg);
    
    title('VI~SPEI (r), LGS');
    cb1=colorbar('Eastoutside');
    crng=[-0.8 0.8];
    colormap(gca, flipud(cbrewer('div','RdBu',33)));
    caxis([crng(1) crng(2)]);
    a=get(gca); gca_w = (a.XLim(2)-a.XLim(1)); gca_h = (a.YLim(2)-a.YLim(1));
    text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,'c)','FontSize',11,'FontName','Arial');
    resizeCB(cb1, cb_xR, cb_yR, cb_wR, cb_tR, '',nan,nan,9);
    
end
