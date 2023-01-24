function scatter_plots(m)
    global lc_dom_grp clmzone_grp;
	global ex_y1 ex_y0 ex_pst;
    %% Scatter plot by land cover
    land_mask=(~isnan(lc_dom_grp)); % All land class
    figure('color','white','Position',[330   422   860   552]);
    plot_extreme(m, ex_y1, ex_y0, ex_pst,land_mask,'SPEI vs. EVI (GS), All','a)');

    %% Plot by land class for warm summer boreal and temperate only
    lc_nam={'ENF','EBF','DF','MF','WS', 'SAV','GRS','WL', 'CRP', 'CRO', 'SIB', 'BAR', 'WAT'};
    label_str={'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)'};
    figure('Position',[134  199  1552  700],'Color','w');
    t = tiledlayout(3,5,'TileSpacing','Compact');
    fi=0;
%     for lc=[5 7 2 3 1 4]
    for lc=1:13
        nexttile;
        fi=fi+1;
        land_mask=((lc_dom_grp==lc));
        try
            plot_extreme(m, ex_y1, ex_y0, ex_pst,land_mask,sprintf('%s',lc_nam{lc}),label_str{fi});
        catch
            continue;
        end
    end
    %%
end