function Recovery_summary_big(alpha_lc,ptype,subgs)

    % We only plot the ecosystem available in both 1st and 2nd GS
    [s1 s2 s3]=size(alpha_lc);
    for i=1:s2
        for j=1:s3
            a1=alpha_lc{1,i,j};
            a2=alpha_lc{2,i,j};

            if isempty(a1)
                continue;
            end

            for gsl=1:3
                for l=1:13
                    if isnan(a1(gsl,2,l)) || isnan(a2(gsl,2,l)) % one of the 1st post-extreme month value is nan, then don't show either of the GS
                        a1(gsl,:,l)=nan;
                        a2(gsl,:,l)=nan;
                    end
                end
            end

            alpha_lc{1,i,j}=a1;
            alpha_lc{2,i,j}=a2;
        end
    end
    
    clim=1;
    use_fig=0;
    
    label_str={'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)','p)'};
    
    figure('Position',[284   259   896   703],'Color','w');
    gap_h=0.05; gap_w=0.04;
    gap=[gap_h gap_w]; marg_h=[0.08 0.08]; marg_w=[0.1 0.05];
    ha = tight_subplot(3,2,gap,marg_h,marg_w);
    
    fig_pos=[1 3 5 2 4 6];
	yscal=0.9; xscal=0.02;
    if subgs==1
        gsl_rng=1:3;
        gs_str='EGS';
    else
        gsl_rng=2:3;
        label_str={'','','a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)','p)'};
        gs_str='LGS';
    end
    
    for gs=1:2
        for gsl=gsl_rng

            fid=fig_pos((gs-1)*3+gsl);
            axes(ha(fid));
            plot_ttl=0;
            plot_xlabel=0;
            plot_ylabel=0;
            if gsl==1
                plot_ttl=1;
            end
            
            if gs==1
               plot_ylabel=1;
            end
            
            if gsl==3
               plot_xlabel=1;
            end
            
            switch gsl
                case 1
                    row_str='1-4 months';
                case 2
                    row_str='5-8 months';
                case 3
                    row_str='9-12 months';
                otherwise
                    error('error!');
            end
            
            dtmp=alpha_lc{gs, clim,subgs};
            
            if subgs==2 && gsl==2
                dtmp(2,6,9)=nan;
            end
            % Only use the one for SPEI (index 4 is 1 in dtmp), using other
            % have similar results
            Recovery_summary(squeeze(dtmp(gsl,:,:,1)),ptype, gs,use_fig, plot_ttl, plot_xlabel, plot_ylabel,'\\alpha_l_e_g (%%)', row_str);
            a=get(gca); gca_w = (a.XLim(2)-a.XLim(1)); gca_h = (a.YLim(2)-a.YLim(1));
            text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,label_str{fid},'FontSize',11,'FontName','Arial');
        end
    end
end