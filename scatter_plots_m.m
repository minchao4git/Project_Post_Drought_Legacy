function [beta_lc]=scatter_plots_m(ptype, gs, gsl, clim)
    global lc_dom_grp clmzone_grp;
	global ex_m ex_pst_m_gs norm_m ex_pst_m;
    global nds lgs_map;

    % Plot by land class for warm summer boreal and temperate only
    lc_nam={'ENF','EBF','DF','MF','WS', 'SAV','GRS','WL', 'CRP', 'CRO', 'SIB', 'BAR', 'WAT'};
    label_str={'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)','p)'};
    figure('Position',[129  354 1562 540],'Color','w');

    gap_h=0.07; gap_w=0.04;
    gap=[gap_h gap_w]; marg_h=[0.08 0.05]; marg_w=[0.1 0.05];
    ha = tight_subplot(2,5,gap,marg_h,marg_w);
    fi=0;

    ntile=0;
	xlabel_str='';
	ylabel_str='';
	show_ledg=0;
    beta_lc=nan(14,13);
    
    gsl_rng={[1 4],[5 8],[9 12]};
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

    for lc=[1:7 9]
        fi=fi+1;

        ntile=ntile+1;
        axes(ha(ntile));

        if fi<2
            switch ptype
                case 'DGVM'
                    ylabel_str='Standardized GPP anomalies (%)';
                case 'SAT'
                    ylabel_str='Standardized EVI anomalies (%)';
                case 'ET'
                    ylabel_str='Standardized ET anomalies (%)';
                otherwise
                    error('Wrong processing type!');
            end
        else
            ylabel_str='';
        end

        if fi>=4
            xlabel_str='SPEI anomalies';
        else
            xlabel_str='';
        end

        if fi==8
            show_ledg=1;
        end

        land_mask=((lc_dom_grp==lc).*((lgs_mask>=lgs1).*(lgs_mask<=lgs2))); % Skip the southern hemisphere cold regions
        [beta_lc(:,lc)]=plot_extreme_m(ex_m, ex_pst_m_gs, norm_m,land_mask,sprintf('%s',lc_nam{lc}),label_str{fi},ylabel_str,xlabel_str, show_ledg,ptype, gs);

    end

end