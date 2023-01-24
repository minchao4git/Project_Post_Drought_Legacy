function r_by_lc_clmzone()

    % ===== R between Veg. and JLI grouped by landclass and climate zone ===== 
    global lc_dom_grp clmzone_grp r_spei_evi;
    global clm_cl_x;

    [s1 s2]=size(clmzone_grp);

    % Calculation
    data_calc=nan(s1,s2,3);
    clm_cl_x=nan(3,8);

    % number of gridpoint
    clm_cl_cnt=nan(3,8,5);

    % group in bar chart
    clm_cl_rpt=nan(3,8,5);
    clm_cl_stderr=nan(8,3,5);
    clm_lc_p=nan(8,3,5);

    nds=1;
    nv=3;
    data_calc(:,:,1,1:nds)=squeeze(r_spei_evi(:,:,14)); % May and June
    data_calc(:,:,2,1:nds)=squeeze(r_spei_evi(:,:,15)); % July and Aug
    data_calc(:,:,3,1:nds)=squeeze(r_spei_evi(:,:,13)); % May~Aug

    for clm=1:7
        for lc=1:8

            % spatial points
            amask=clmzone_grp; amask(amask~=clm)=nan;amask(amask==clm)=1;
            bmask=lc_dom_grp; bmask(bmask~=lc)=nan;bmask(bmask==lc)=1;
            dtmp=data_calc.*repmat(amask,[1 1 nv nds]).*repmat(bmask,[1 1 nv nds]);

            for v=1:3  % variable presented in the calculation
                % dataset level
                for dsi=1:nds
                    atmp=squeeze(dtmp(:,:,v,dsi));
                    atmp=atmp(~isnan(atmp));

                    % The mean across space
                    clm_cl_rpt(clm,lc,v,dsi)=nanmean(atmp);
                end
                
                % mean across dataset
                btmp=squeeze(nanmean(dtmp(:,:,v,:),4));
                btmp=btmp(~isnan(btmp));
                
                % number of gridpoints
                clm_cl_cnt(clm,lc,v)=length(btmp);
                % Standard error across space
                clm_cl_stderr(clm,lc,v)=std(btmp,0,1,'omitnan');

                % Confidence interval across space
                clm_lc_p(clm,lc,v)= clm_cl_stderr(clm,lc,v); % Here I use 1 SD as error bar instead
            end
            a=nanmean(dtmp(:,:,1,:),4);a=a(~isnan(a));
            b=nanmean(dtmp(:,:,2,:),4);b=b(~isnan(b));
            clm_cl_sig(clm,lc)=ttest2(a(:),b(:));
            
        end
    end


    % ===== plotting ===== 
    color_lc_grp=[[0 129 27];
                      [124 255 100]; 
                      [142 204 51];
                      [214 236 163]; 
                      [244 181 120]; 
                      [0 104 150]; 
                      [255 236 131]; 
                      [170 255 255]
        ]/255;

    lc_names={'EF','DF','MF','WS', 'GRS','WL', 'CRP', 'SIB'};
    title_nam={'Boreal','Boreal','Temperate','Temperate','Medit.','Medit.'};
    clmz_name={'Warm-summer','Cold-summer','Warm-summer','Cold-summer', 'Warm-summer', 'Cold-summer'};
    labels={'a)','b)','c)','d)','e)','f)','g)'};
    nv_show=1:4;
    figure('Position',[2475 247 467 499],'Color','w');

    t = tiledlayout(3,1,'TileSpacing','Compact');
    % --- land class ---
    for clm_grp=3:4

        % Sorted by number of grid points in 1st season map
        % [B sort_i(:,clm_grp)]=sort(clm_cl_cnt(clm_grp,:,1),'descend','MissingPlacement','last');
        % sorted by annual values of land class
        [B sort_i(:,clm_grp)]=sort(nanmean(clm_cl_rpt(clm_grp,:,3,:),4),'descend','MissingPlacement','last');
        xtmp=sort_i(:,clm_grp);

        % --- climatezone A ---
        nexttile;

        n=0;
        x=nan;
        for i=1:length(xtmp)
          if ~isnan(nanmean(clm_cl_rpt(clm_grp,xtmp(i),3,:),4)) && clm_cl_cnt(clm_grp,xtmp(i),3) >= 10 % skim based on the CV values
             n=n+1;
             x(n)=xtmp(i);
             clm_cl_x(clm_grp,n)=xtmp(i); % output for decomp_VegClm_relation_plot
          end
        end

        hold on
        width=0.8;
        xoff=0;

        % Show relevant varaible in the barchat column
        b1=bar((1:length(x))-xoff, squeeze(nanmean(clm_cl_rpt(clm_grp,x,[3 1 2],:),4)),width);
        b1(2).FaceColor = [22 170 116]/255;
        b1(3).FaceColor = [255 232 99]/255;
        
        % First error bar in the group
        x_bar=[b1(1,1).XEndPoints]';
        y_bar=[b1(1,1).YEndPoints; b1(1,2).YEndPoints; b1(1,3).YEndPoints]';
        er = errorbar(x_bar, squeeze(nanmean(clm_cl_rpt(clm_grp,x,[3],:),4)), squeeze(clm_lc_p(clm_grp,x,[3])));
        for v=1:1
            er(1,v).Color = [0 0 0];
            er(1,v).LineStyle = 'none'; % remove the plot line
        end
        
        % The second and the third bar in the group
        x_bar_2_3=[b1(1,2).XEndPoints; b1(1,3).XEndPoints]';

        for i=1:length(x)
            er_2_3 = errorbar(x_bar_2_3(i,:), squeeze(nanmean(clm_cl_rpt(clm_grp,x(i),[1 2],:),4)), squeeze(clm_lc_p(clm_grp,x(i),[1 2])));
            er_2_3.Color = [0 0 0];
            er_2_3.LineStyle = 'none'; % remove the plot line
            
            % solid line for significant values
            if clm_cl_sig(clm_grp,x(i))==0
                er_2_3.Bar.LineStyle = 'dashed';
            elseif clm_cl_sig(clm_grp,x(i))==1
                er_2_3.Bar.LineStyle = 'solid';
            end
        end

        set(gca, 'XTick', 1:length(x));
        set(gca, 'XTickLabel',{lc_names{x}});
        xtickangle(0);

        plot([-100 100],[0 0],'-k');

        hold off
        box on;
        ylabel(sprintf('%s (r)',clmz_name{clm_grp}));

        if mod(clm_grp,2)==1
           title(sprintf('SPEI~EVI(r), %s',title_nam{clm_grp}),'Fontsize',12);
        end

        if mod(clm_grp,2)==0
            xlabel('Landclass');
        end
        xlim([0.2 6.8]);
        if clm_grp==1
            ylim([-0.4 0.6]);
        elseif clm_grp==2
            ylim([-0.4 0.6]);
        else
            ylim([-0.4 0.6]);
        end

        if clm_grp==4
            legend(b1,{'May-Aug', 'May, Jun', 'Jul, Aug'});
            legend boxoff;
        end
        % labels for sub-plots
        yscal=0.92; xscal=0.03;
        a=get(gca);
        gca_w = (a.XLim(2)-a.XLim(1));
        gca_h = (a.YLim(2)-a.YLim(1));
        text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,labels{clm_grp},'FontSize',11,'FontName','Arial');
        
    end
end