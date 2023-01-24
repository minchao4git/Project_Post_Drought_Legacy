function Recovery_summary(beta_lc, ptype, gs)
    lc_nam={'ENF','EBF','DF','MF','WS', 'SAV','GRS','WL', 'CRP', 'CRO', 'SIB', 'BAR', 'WAT'};
    eco_color=[
               [0 120 0];     % ENF
               [85 170 0];    % EBF
               [170 255 0];   % DF
               [170 85 255];  % MF 
               [170 0 127];   % WS
               [170 170 0];   % SAV
               [255 255 0];   % GRS
               [0 85 255];    % WL
               [255 170 0];   % CRP
               [255 170 255];  % CRO
               [0 170 255];   % SIB
               [120 120 120]; % BAR
               [250   0   0]; % WAT
               ]./255;

    ec_color={[0 120 0]/255, [255 170 0]/255};
    figure('color','w', 'Position',[385   293   607   436]);
    hold on;
    nec=0;
    line_h=nan(13,1);
    recover_tim=nan(13,1);
    valid_ec=[];
    mon_stack=zeros(13,1);
    
    switch ptype
        case 'SAT'
            yrng=[-18 8];
            ylim(yrng);
            ylabel('Response of VI to neutral climate after extremes (%)');
            dot_gap=0.6;
        case 'ET'
            yrng=[-18 8];
            ylim(yrng);
            ylabel('Response of ET to neutral climate after extremes (%)');
            dot_gap=0.6;
        case 'DGVM'
            yrng=[-40 8];
            ylim(yrng);
            ylabel('Response of GPP to neutral climate after extremes (%)');
            dot_gap=1.2;
        otherwise
            error('Wrong processing type!');
    end
    
    switch gs
        case 1
            pstm_rng=1:13;
            xshift=0;
            gs_str=sprintf('1st\_%s',ptype);
            gs_str_ttl=sprintf('1st');
        case 2
            pstm_rng=1:13;
            xshift=1;
            beta_lc(1,:)=nan;
            gs_str=sprintf('2nd\_%s',ptype);
            gs_str_ttl=sprintf('2nd');
        otherwise
            pstm_rng=1:13;
            xshift=0;
            gs_str='All';
    end
    
    for lc=1:13
        % Skip the ecosystem with all nan months
        if sum(isnan(beta_lc(1:13,lc)))==13
            continue;
        end
        
        % Skip the ecosystem without post-extreme data
        if isnan(beta_lc(2,lc))
            continue;
        end
        
        valid_ec=[valid_ec lc];
    end
    
    % grass and tree
    grs_ec_tmp=[7 9:10 12];
    tre_ec_tmp=[1:6];
    
    % Only include the valid land class
    grs_ec=[];tre_ec=[];
    for v=1:length(grs_ec_tmp)
        if any(valid_ec==grs_ec_tmp(v))
            grs_ec=[grs_ec grs_ec_tmp(v)];
        end
    end
    for v=1:length(tre_ec_tmp)
        if any(valid_ec==tre_ec_tmp(v))
            tre_ec=[tre_ec tre_ec_tmp(v)];
        end
    end
    
    ec_grp={tre_ec;grs_ec};
    
    for ec=1:2
        pt_pos_low=nan(1,2);
        pt_pos_up=nan(1,2);
        
        ec_mn=nanmin(beta_lc(pstm_rng,ec_grp{ec}),[],2); % 1:13, not including normal month
        ec_mean=nanmean(beta_lc(pstm_rng,ec_grp{ec}),2);
        ec_mx=nanmax(beta_lc(pstm_rng,ec_grp{ec}),[],2);

        npstm=length(ec_mn);
        % align the length of the growing season
        for m=1:npstm
            if (isnan(ec_mx(m)) || isnan(ec_mean(m)) || isnan(ec_mn(m)))
                ec_mx(m)=nan;
                ec_mean(m)=nan;
                ec_mn(m)=nan;
            end
        end

        % get lower and upper bound of the polygon
        n=0;
        for m=1:npstm
            if ~isnan(ec_mn(m))
                n=n+1;
                pt_pos_low(n,1)=m;
                pt_pos_low(n,2)=ec_mn(m);
            end
        end

        n=0;
        for m=1:npstm
            if ~isnan(ec_mx(m))
                n=n+1;
                pt_pos_up(n,1)=m;
                pt_pos_up(n,2)=ec_mx(m);
            end
        end

        v1 = [pt_pos_low; flipud(pt_pos_up)];
        f1 = 1:size(pt_pos_low,1)*2;
        patch('Faces',f1,'Vertices',v1,'FaceColor',ec_color{ec},'FaceAlpha',.3, 'EdgeColor','none');
        plot(((1:size(pt_pos_low)+xshift)), ec_mean(~isnan(ec_mean)), 'color',ec_color{ec},'LineWidth',2);
    end
    
    for lc=valid_ec
        
        plot(beta_lc(pstm_rng,lc),'color',eco_color(lc,:));
        line_h(lc)=scatter(pstm_rng,beta_lc(pstm_rng,lc),40,'MarkerFaceColor', eco_color(lc,:),'MarkerEdgeColor', [0 0 0]/255);

        % show the recover time
        if gs~=2
            for m=pstm_rng
                if beta_lc(m,lc) >= beta_lc(14,lc)
                    recover_tim(lc,1)=m;

                    % push to month stack, this is to avoid overlap of the scatter
                    % points
                    mon_stack(m)=mon_stack(m)+1;
                    break;
                end
            end

            if mon_stack(m)==1
                plot([recover_tim(lc,1) recover_tim(lc,1)],[yrng(1)+1+dot_gap*(mon_stack(m)-1) beta_lc(m,lc)],'color',[0.5 0.5 0.5],'LineStyle','--');
            end
            scatter(recover_tim(lc,1),yrng(1)+1+dot_gap*(mon_stack(m)-1),40,'MarkerFaceColor', eco_color(lc,:),'MarkerEdgeColor', [0 0 0]/255);
        end
        
    end
    plot([-100 100],[0 0],'k--');
    norm_min=nanmin(squeeze(beta_lc(14,valid_ec)));
    norm_mean=squeeze(nanmean(beta_lc(14,valid_ec)));
    norm_max=nanmax(squeeze(beta_lc(14,valid_ec)));
    % figure;
    plot([-100 100],[norm_max norm_max],'b--');
    plot([-100 100],[norm_mean norm_mean],'b-','LineWidth',2);
    plot([-100 100],[norm_min norm_min],'b--');
    xlim([0 14]);

    legend(line_h(valid_ec), lc_nam(valid_ec), 'Location','Southeast');
    legend boxoff;

    set(gca,'XTick',1:13);
    if gs==2
        set(gca,'XTickLabel',{'2nd GS','1','2','3','4','5','6','7','8','9','10','11','12'});
    else
        set(gca,'XTickLabel',{'Ex. mon.','1','2','3','4','5','6','7','8','9','10','11','12'});
    end
    % set(gca,'XTickLabelRotation',45);

    box on;
    title(sprintf('Response to extreme and recover time (%s GS)',gs_str_ttl));
    xlabel('Post-extreme growing season months');

end