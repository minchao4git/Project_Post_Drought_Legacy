function [alpha_lc]=scatter_plots_m_gsl(ptype, gs, clim)
    global lc_dom_grp clmzone_grp;
	global ex_m ex_pst_m_gs norm_m ex_pst_m;
    global nds lgs_map;

    % Plot by land class for warm summer boreal and temperate only
    lc_nam={'ENF','EBF','DF','MF','WS', 'SAV','GRS','WL', 'CRP', 'CRO', 'SIB', 'BAR', 'WAT'};
    label_str={'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)','p)','r)','s)','t)','u)'};
    figure('Position',[134  207  1512  649],'Color','w');

    gap_h=0.07; gap_w=0.02;
    gap=[gap_h gap_w]; marg_h=[0.08 0.05]; marg_w=[0.04 0.04];
    ha = tight_subplot(3,7,gap,marg_h,marg_w);
    fi=0;
	yscal=0.92; xscal=0.04;

	xlabel_str='';
	ylabel_str='';
	show_ledg=0;
    alpha_lc=nan(3,14,13);
    
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

        ntile=1;
        for lc=[1:7 9]
            axes(ha(ntile+(gsl-1)*7));
            
            if ntile==1
                switch gsl
                    case 1
                        ylabel_str='GS length 1-4 months';
                    case 2
                        ylabel_str='GS length 5-8 months';
                    case 3
                        ylabel_str='GS length 9-12 months';
                    otherwise
                        error('error!');
                end
            else
                ylabel_str='';
            end

            if gsl==3
                xlabel_str='SPEI anom. (z-score)';
            else
                xlabel_str='';
            end

            if ntile==6 && gsl==2
                show_ledg=1;
            else
                show_ledg=0;
            end

            if gsl==3 && lc==2 % very small amount of EBF, skip this
                continue;
            end
            
            land_mask=((lc_dom_grp==lc).*((lgs_mask>=lgs1).*(lgs_mask<=lgs2))); % Skip the southern hemisphere cold regions
            ifnopst=false;
            [alpha_lc(gsl,:,lc), ifnopst]=plot_extreme_m(ex_m, ex_pst_m_gs, norm_m,land_mask,sprintf('%s',lc_nam{lc}),ylabel_str,xlabel_str, show_ledg,ptype, gs);
            if ~ifnopst
                fi=fi+1;
                a=get(gca); gca_w = (a.XLim(2)-a.XLim(1)); gca_h = (a.YLim(2)-a.YLim(1));
                text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,label_str{fi},'FontSize',11,'FontName','Arial');
                ntile=ntile+1;
            end

        end
    end
end