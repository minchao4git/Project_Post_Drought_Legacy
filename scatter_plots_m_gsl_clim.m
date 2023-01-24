function [alpha_lc, beta_lc, mdinfo_lc, p_slope_diff, b_diff_m_all_wg]=scatter_plots_m_gsl_clim(ptype, gs, clim)
    global lc_dom_grp clmzone_grp;
	global ex_m ex_pst_m_gs norm_m ex_pst_m;
    global nds lgs_map;
    global DATA_ESA_CCI_LC_out;

    % Plot by land class for warm summer boreal and temperate only
    lc_nam={'ENF','EBF','DF','MF','WS', 'SAV','GRS','WL', 'CRP', 'CRO', 'SIB', 'BAR', 'WAT'};
    label_str={'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)','p)','r)','s)','t)','u)'};
    ttl_str={'SPEI', 'Temp.','Rad.','SM'};
    %        gsl, gsm, lc, climvars
    alpha_lc=nan(3,14,13, 4);
    beta_lc=nan(3,14,13, 4);
    mdinfo_lc=cell(3,26,13, 4); % 14 months + 12 normal months
    %        gsl, lc, climvars
    p_slope_diff=nan(3,13,4); % land class level
    
    % group level
    px_m_all_w_x=[];
    px_m_all_w_y=[];
    px_m_all_g_x=[];
    px_m_all_g_y=[];
    
    nor_m_all_w_x=[];
    nor_m_all_w_y=[];
    nor_m_all_g_x=[];
    nor_m_all_g_y=[];
    
    % gsl, difference (p) between post-extreme and normal months, woody(1) and grass(2),
    % climate, and parameters(b1,b2, standard error, p-value)
    b_diff_m_all_wg=nan(3,2,4,4); % difference
    
    clim_n=0;
    
    for clim_id=3:6 % spei, t2m, ssrd, soil moisture
    
        clim_n=clim_n+1;
        
        fh=figure('Position',[134  207  1512  649],'Color','w');
        fh.Visible=false;
        gap_h=0.07; gap_w=0.02;
        gap=[gap_h gap_w]; marg_h=[0.08 0.05]; marg_w=[0.04 0.04];
        ha=tight_subplot(3,7,gap,marg_h,marg_w);
        set(ha,'Visible','off');
        fi=0;
        yscal=0.92; xscal=0.04;

        xlabel_str='';
        ylabel_str='';
        show_ledg=0;

        gsl_rng={[1 4],[5 8],[9 12]};

        for gsl=1:3
            lgs1=gsl_rng{gsl}(1);
            lgs2=gsl_rng{gsl}(2);

            switch ptype
                case 'DGVM'
                    ndatasrc=1;
                case 'SAT'
                    ndatasrc=nds;
                case 'ET'
                    ndatasrc=nds;
                otherwise
                    error('Wrong processing type!');
            end

            lgs_mask=lgs_map(:,:,ndatasrc);
            lgs_mask(isnan(lgs_mask))=0;

            % landclass level
            ntile=1;
            
            for lc=[1:7 9]
                ha_tmp=ha(ntile+(gsl-1)*7);
                ha_tmp.Visible=false;
                
                axes(ha_tmp);

                set(gcf, 'Visible', 'off');
                set(fh, 'Visible', 'off');

                if ntile==1
                    switch gsl
                        case 1
                            row_str='GS length 1-4 months';
                        case 2
                            row_str='GS length 5-8 months';
                        case 3
                            row_str='GS length 9-12 months';
                        otherwise
                            error('error!');
                    end
                else
                    row_str='';
                end

                if gsl==3
                    xlabel_str=sprintf('%s anom. (z-score)',ttl_str{clim_n});
                else
                    xlabel_str='';
                end

                if ntile==5 && gsl==3
                    show_ledg=1;
                else
                    show_ledg=0;
                end

                if gsl==3 && lc==2 % skip EBF, only on post-extreme months
                    continue;
                end

                land_mask=((lc_dom_grp==lc).*((lgs_mask>=lgs1).*(lgs_mask<=lgs2))); % Skip the southern hemisphere cold regions
                
                if lc==9 % croplands
                    land_mask=land_mask.*(DATA_ESA_CCI_LC_out<=0.5); % TEST-MW: Consider irrigated croplands
                end
                
                ifnopst=false;
                [alpha_lc(gsl,:,lc, clim_n), beta_lc(gsl,:,lc, clim_n), ifnopst, mdinfo, p_slope_diff(gsl,lc, clim_n), px_m_1, nor_m_1]=...
                    plot_extreme_m_clim(ex_m, ex_pst_m_gs, norm_m,land_mask,sprintf('%s',lc_nam{lc}),row_str,xlabel_str, show_ledg,ptype, gs, clim_id);
                
                % woody and grass
                if length(px_m_1) > 0
                    switch lc
                        case {1,2,3,4,5,6}
                            % woody
                            px_m_all_w_x=[px_m_all_w_x px_m_1(1,:)];
                            px_m_all_w_y=[px_m_all_w_y px_m_1(2,:)];

                            nor_m_all_w_x=[nor_m_all_w_x nor_m_1(1,:)];
                            nor_m_all_w_y=[nor_m_all_w_y nor_m_1(2,:)];
                        case {7,9}
                            % grass
                            px_m_all_g_x=[px_m_all_g_x px_m_1(1,:)];
                            px_m_all_g_y=[px_m_all_g_y px_m_1(2,:)];

                            nor_m_all_g_x=[nor_m_all_g_x nor_m_1(1,:)];
                            nor_m_all_g_y=[nor_m_all_g_y nor_m_1(2,:)];
                        otherwise
                            error('Wrong class!');
                    end
                end
                
                % can not transfer the whole cell array from a function,
                % needs to do it one by one
                for i=1:26
                    mdinfo_lc{gsl,i,lc, clim_n}=mdinfo{i};
                end
                if ~ifnopst
                    fi=fi+1;
                    a=get(gca); gca_w = (a.XLim(2)-a.XLim(1)); gca_h = (a.YLim(2)-a.YLim(1));
                    text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,label_str{fi},'FontSize',11,'FontName','Arial');
                    ntile=ntile+1;
                end

            end % landclass level
            
            
            % Compute the slope and their differences
            % woody
            if length(px_m_all_w_x)>0
                md1=fitlm(px_m_all_w_x, px_m_all_w_y,'RobustOpts','on');
                md2=fitlm(nor_m_all_w_x, nor_m_all_w_y,'RobustOpts','on');
                b1=md1.Coefficients.Estimate(2);
                b2=md2.Coefficients.Estimate(2);
                p1=md1.Coefficients.pValue(2);
                p2=md2.Coefficients.pValue(2);
                SE_b1=md1.Coefficients.SE(2);
                SE_b2=md2.Coefficients.SE(2);
                SE_b_diff=sqrt(SE_b1^2 + SE_b2^2);
                T=abs((b1-b2))/SE_b_diff;
                DF=(md1.NumObservations+md2.NumObservations-4);

                b_diff_m_all_wg(gsl,1,clim_n, 1)=b1;
                b_diff_m_all_wg(gsl,1,clim_n, 2)=b2;
                b_diff_m_all_wg(gsl,1,clim_n, 3)=SE_b_diff;
                b_diff_m_all_wg(gsl,1,clim_n, 4)=2*(1-tcdf(T,DF));
                b_diff_m_all_wg(gsl,1,clim_n, 5)=p1;
                b_diff_m_all_wg(gsl,1,clim_n, 6)=p2;
            end

            % grass
            if length(px_m_all_g_x)>0
                md1=fitlm(px_m_all_g_x, px_m_all_g_y,'RobustOpts','on');
                md2=fitlm(nor_m_all_g_x, nor_m_all_g_y,'RobustOpts','on');
                b1=md1.Coefficients.Estimate(2);
                b2=md2.Coefficients.Estimate(2);
                p1=md1.Coefficients.pValue(2);
                p2=md2.Coefficients.pValue(2);
                SE_b1=md1.Coefficients.SE(2);
                SE_b2=md2.Coefficients.SE(2);
                SE_b_diff=sqrt(SE_b1^2 + SE_b2^2);
                T=abs((b1-b2))/SE_b_diff;
                DF=(md1.NumObservations+md2.NumObservations-4);

                b_diff_m_all_wg(gsl,2,clim_n,1)=b1;
                b_diff_m_all_wg(gsl,2,clim_n,2)=b2;
                b_diff_m_all_wg(gsl,2,clim_n,3)=SE_b_diff;
                b_diff_m_all_wg(gsl,2,clim_n,4)=2*(1-tcdf(T,DF));
                b_diff_m_all_wg(gsl,2,clim_n,5)=p1;
                b_diff_m_all_wg(gsl,2,clim_n,6)=p2;
            end
        end % gsl
        
        % Save the *.fig first, exported as image file after necessay
        % layout adjustification
        fh.Visible=true;
        savefig(fh,sprintf('./figures/scatter_plots_gs%d_clim%d.fig', gs, clim_id));
        close(fh);
        
    end % climate variable
end