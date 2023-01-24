function LandClass_ClimteZone()
    global DATA_LC_out DATA_KG_CLMZ_05rs_out;
    global lc_dom_grp clmzone_grp;
    global ax_show;
    %% Climate zone and landclass grouping
    % ====> IGBP Land cover
    % > The calculation is based on the regridded dataset from MODIS IGBP land
    % cover type
    % > Two appraoches are used to achieve the dominant landcover type, and
    % both yield similar results
    % > The landcover type based on 0.5 grid and 0.05 grid are quite similar,
    % and also consistent with landclass from IGBP: 
    % https://e4ftl01.cr.usgs.gov/MOTA/MCD12C1.006/2008.01.01/BROWSE.MCD12C1.A2008001.006.2018053184623.1.jpg

    % ----------------------------
    %  Eight groups of land class 
    % ----------------------------
    % Skip igbp 1: water boday

    % Method I: Freqency approach: get the stable land cover
    [s1 s2 s3 s4]=size(DATA_LC_out);
    nyr_lc=s3;
    for i=1:s1
        for j=1:s2
            for y=1:nyr_lc
                [dump lc_dom_ts(i,j,y)]=max(DATA_LC_out(i,j,y,:)); % get the maximun fraction
            end
        end
    end

    lc_dom=mode(lc_dom_ts,3); % get the most freqent value (not the function mod() :-) )
    lc_dom_stable=nan(s1,s2);
    for i=1:s1
        for j=1:s2

            if sum((find(lc_dom_ts(i,j,:)==lc_dom(i,j))>0))==nyr_lc
                lc_dom_stable(i,j)=lc_dom(i,j);
            end
        end
    end
    % 
    % % Method II: Mean approach
    % for i=1:s1
    %     for j=1:s2
    %         [dump lc_dom(i,j)]=max(nanmean(DATA_LC_out(i,j,:,:),3));
    %     end
    % end

    % grouping
    lc_dom_grp=nan(s1,s2);
    for i=1:s1
        for j=1:s2
            switch lc_dom(i,j) % !!! Note: not the IGBP ID, is the variable index from preproc_MODIS_Landcover_Global
                % Evergreen needleleaf forests
                case {2}
                    lc_grp=1;
                % Evergreen broadleaf forests
                case {3}
                    lc_grp=2;
                % Decidous
                case {4 5}
                    lc_grp=3;
                % Mixed forest
                case {6}
                    lc_grp=4;
                % Closed shrubland, Open shrubland and Woody savannas
                case {7,8,9}
                    lc_grp=5;
                % Savannas
                case {10}
                    lc_grp=6;
                % Grasslands
                case {11}
                    lc_grp=7;
                % Permanent wetland
                case {12}
                    lc_grp=8;
                % Croplands (Skip the Urban and builtup area)
                case {13}
                    lc_grp=9;
                % Croplands and natural vegetation mosaic (Skip the Urban and builtup area)
                case {15}
                    lc_grp=10;
                % Snow and ice
                case {16}
                    lc_grp=11;
                % Barren sparsely veg.
                case {17}
                    lc_grp=12;
                % Water bodies
                case {1}
                    lc_grp=13;
                otherwise
                    lc_grp=nan;
            end
            lc_dom_grp(i,j)=lc_grp;
        end
    end

    color_lc_grp=[[0 129 27]; 
                  [0 170 0]; 
                  [124 255 100]; 
                  [158 121 237]; 
                  [180 216 19];
                  [244 181 120]; 
                  [214 236 163]; 
                  [0 104 150];
                  [255 236 131];
                  [255 170 255];
                  [170 255 255];
                  [184 184 184];
                  [0 85 255];
                  [255 255 255];
    ]/255;
    %%
    % -----------------------------------
    %     Three groups of climate zones 
    % -----------------------------------
     % grouping (see the function get_cls_id() for the explaination of climate ID)
    clmzone_grp=nan(s1,s2);
    for i=1:s1
        for j=1:s2
            switch DATA_KG_CLMZ_05rs_out(i,j)
                % Af (11): Tropical rainforest climate
                % Am (12): Tropical momsoon climate
                case {11, 12}
                    cl_grp=1;
                % Precip.-limited GS
                case {13, 14, 22, 27, 34, 35, 36}
                    cl_grp=2;
                % Temp.-limited GS
                case {21, 26, 31, 32, 33, 37, 38, 39, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 61, 62}
                    cl_grp=3;
                otherwise
                    cl_grp=nan;
            end

            clmzone_grp(i,j)=cl_grp;
        end
    end
    % You can skip the group with number of grid point < 10
    color_clm_grp=[[0 170 127];
                   [255 197 81];
                   [0 85 255];
                   [255 255 255];
    ]/255;

    %%
    figure('color','white','Position',[259 273 1058  602]);
    gap_h=0.05; gap_w=0.005;
    gap=[gap_h gap_w]; marg_h=[0.08 0.05]; marg_w=[0.1 0.15];
    ha = tight_subplot(2,1,gap,marg_h,marg_w);
    yscal=0.92; xscal=0.03;
    cb_xR=1.0;
    cb_yR=1.0;
    cb_wR=0.6;
    cb_hR=1.0;
    ax_show.frame=0;

    % --- Land class ---
    axes(ha(1));
    geoplot(ax_show, lc_dom_grp');

    colormap(gca,color_lc_grp);
    caxis([0.5 14.5]);
    cbh1=colorbar('eastoutside');
    cbh1.Ticks = [1:14];
    cbh1.TickLabels={'ENF','EBF','DF','MF','WS', 'SAV','GRS','WL', 'CRP', 'CRO', 'SIB', 'BAR', 'WAT', 'No data'};
    a=get(gca);
    gca_w = (a.XLim(2)-a.XLim(1));
    gca_h = (a.YLim(2)-a.YLim(1));
    text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,'a)','FontSize',11,'FontName','Arial');
    title('Landclass','FontSize',11);
    resizeCB(cbh1, cb_xR, cb_yR, cb_wR, cb_hR, '',0.32,30,9);

end
