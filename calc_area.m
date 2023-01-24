
lat_v=23:0.5:73.5;
lon_v=-179.75:0.5:179.75;
[LON_A, LAT_A]=meshgrid(lon_v,lat_v);
A=cdtarea(LAT_A, LON_A,'km2');

[s1,s2]=size(DATA_KG_CLMZ_05rs_out);
 % grouping
clmzone_grp=nan(s1,s2);
for i=1:s1
    for j=1:s2
        switch DATA_KG_CLMZ_05rs_out(i,j)
            % Dfc, ET (snow, fully humid, warm and cool summer, polar tundra)
            case {43, 62}
                cl_grp=1;
            % Dfb, Cfb, Cfc (warm temperate, fully humid, hot-warm-cool summer)
            % 35: Csb
            case {42, 32, 33,       26, 41, 42, 35, 51, 50, 49, 21, 46, 27}
                cl_grp=2;
            % Cfa, Csa, Csc (warm temperate, summer dry, hot-warm-cool summer)
            % BSk (arid steppe cold arid)
            case {22, 31, 34,   35,  36, 26, 37, 38}
                cl_grp=3;
            case {13,14}
                cl_grp=4;
            otherwise
                cl_grp=nan;
        end

        clmzone_grp(i,j)=cl_grp;
    end
end
figure;
imagesc(flipud(clmzone_grp'));
caxis([0 4]);

figure;
subplot(4,1,1)
lat=squeeze(lats_sig(:,:,m)).*(etype_EGS_neg(:,:)|etype_LGS_neg(:,:));
lon=squeeze(lons_sig(:,:,m)).*(etype_EGS_neg(:,:)|etype_LGS_neg(:,:));
b=lat';
imagesc(~isnan(b))
subplot(4,1,2)
imagesc((abs(bg_show)>0))
subplot(4,1,3)
imagesc(A)
subplot(4,1,4)
imagesc(flipud(rot90(clmzone_grp)));caxis([0 4]);



clmzone=flipud(rot90(clmzone_grp));
sum(sum((~isnan(b).*A),1),2)/sum(sum((abs(bg_show)>0).*A,1),2)




