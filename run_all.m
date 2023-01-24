% load('./data/Workspace_2021-10-18_SAT_DataInput.mat', '*');
%% ---- MATLAB NC Lib ---
clear global *;
clearvars;
sys_id='Tiberino'; % Linux_Dell, ra_T420, tiberino

global data_src;
global lib_path;
global wrk_dir lib_local_dir;

switch sys_id
    case 'Linux_Dell' 
        lib_path='/home/minchao/work/matlab_lib/';
        data_src='/media/minchao/DS1/DATA/obs/';
        wrk_dir='.';
    case 'Linux_Dell_dlocal'
        lib_path='/home/minchao/work/matlab_lib/';
        data_src='/Data/obs';
        wrk_dir='/home/minchao/work/GoogleDriver/Research/a1_Studies_and_projects/Legacy_effects/';
        lib_local_dir=sprintf('%s/lib/',wrk_dir);
    case 'ra_T420'
        lib_path='/home/minchao/NewWork/matlab_lib/';
        data_src='/media/minchao/DS1/DATA/obs/';
        wrk_dir='.';
    case 'Tiberino'
        lib_path='/home/minchao/work/matlab_lib/';
        data_src='/data/fyris/minchao/';
        wrk_dir='.';
        lib_local_dir=sprintf('%s/lib/',wrk_dir);
    case 'Bi'
        wrk_dir='.';
        lib_local_dir=sprintf('%s/lib/',wrk_dir);
        lib_path='/home/sm_minwu/work/matlab_lib';
        data_src='/nobackup/rossby20/sm_minwu/DATA/OBS/';
    otherwise
        error('Not valid system!!');
end

addpath(lib_path);
addpath(lib_local_dir);
addpath(sprintf('%s/sm_grini/nctoolbox/nctoolbox-1.1.1/',lib_path));
setup_nctoolbox;
addpath(sprintf('%s/CDT-Chad_Greene/v1.01/cdt/',lib_path));
addpath(sprintf('%s/CDT-Chad_Greene/v1.01/cdt/cdt_data/',lib_path));
addpath(sprintf('%s/cbrewer/',lib_path));
addpath(sprintf('%s/export_fig/',lib_path));
addpath(sprintf('%s/ATSA_Meko/',lib_path));
addpath(sprintf('%s/freezeColors/',lib_path));

global domain_def;
global chosen_prd chosen_prd_dgvm;
global nyr nyr_gs;
global nds ds_s ds_e nmd mdi;
global nyr_dgvm nyr_dgvm_gs;
global stime;
global TT;

global DATA_CRU_05rs_out;
global sos_map pos_map eos_map lgs_map;

global DATA_VI_05rs_out DATA_VI_05rs_sm_out;
global DATA_VI_05rs_bw_out DATA_VI_05rs_bw_sm_out;
global DATA_VI_05rs_bw_out;

% VI datasets
global DATA_05rs_out_gs;
global DATA_05rs_dtds_out DATA_05rs_dtds_out_gs DATA_05rs_dtds_out_sgs;
global DATA_05rs_dtdspw_out;

% DGVMs dataset
global DATA_DGVM_05rs_out;
global DATA_DGVM_05rs_dtds_out DATA_DGVM_05rs_dtds_out_gs;

% EC site dataset
global DATA_EC_out_gs;
global DATA_EC_dtds_out DATA_EC_dtds_out_gs;
global DATA_EC_dtdspw_out;

global DATA_SPEI_out DATA_SPEI_out_s;
global DATA_LC_out DATA_KG_CLMZ_05rs_out;
global lc_dom_grp clmzone_grp all_subgrp;
global ax_show;
global r_spei_evi p_spei_evi;
global r_spei_evi_asyn p_spei_evi_asyn;
global ex_y1 ex_y0 ex_pst;
global ex_y1_s ex_y0_s ex_pst_s;
global DATA_ESA_CCI_LC_out;

global DATA_ERA5_05rs_out;
global DATA_Evap_05rs_out DATA_SMRoot_out;
global lons lats;

global r_spei_gpp p_spei_gpp;
global DATA_05rs_dtds_out_gs_inc;
global ex_m ex_pst_m norm_m eem_map;
global ex_m_s ex_pst_m_s norm_m_s eem_map_s;

% Domain definition, generate data mask
%            Lon.W., Lon.E., Lat.S, Lat.N
% domain_def = [-13.0 33.0 34.0 72.5 ];
% domain_def = [-13.0 35.0 34.0 72.5 ];
% domain_def = [-180.0 180.0 -60.0 80 ];
domain_def = [-180.0 180.0 23.0 74.0 ]; % Tropics of cancer
%            Lon.W., Lon.E., Lat.S, Lat.N

% For Satellite dataset
chosen_prd=2001:2020;
nyr=chosen_prd(end)-chosen_prd(1)+1;
nyr_gs=nyr-1;

% For DGVMs dataset
chosen_prd_dgvm=1981:2019;
nyr_dgvm=chosen_prd_dgvm(end)-chosen_prd_dgvm(1)+1;
nyr_dgvm_gs=nyr_dgvm-1;

% For FLUXNET dataset
global chosen_prd_s nyr_s nyr_s_gs;
chosen_prd_s=1991:2019;
nyr_s=chosen_prd_s(end)-chosen_prd_s(1)+1;
nyr_s_gs=nyr_s-1;

% then the period for SPEI:
global chosen_prd_spei;
chosen_prd_spei=chosen_prd_dgvm(1):chosen_prd(end);

% Monthly calendar
stime = datetime(chosen_prd(1),1,15);
tstep = calmonths(1);
TT = timetable('Size',[nyr*12 1],'VariableTypes',{'double'},...
               'TimeStep',tstep,'StartTime',stime);
ds_s=1; ds_e=2; nds=ds_e-ds_s+1; 
nmd=3;
mdi=1; % 4: Double logistic from TIMESAT (NOTE: We don't deal with this at this moment)

% Parameters for calculating trend
global window_wd nw;
window_wd=20;
nw=(nyr-window_wd)+1;
nw_gs=nw-1;

% Setting for plotting
ax_show=struct('south_s',domain_def(3),'north_s',domain_def(4),'west_s',domain_def(1),'east_s',domain_def(2),'resolution_s',0.5, ...
               'south_d',domain_def(3),'north_d',domain_def(4),'west_d',domain_def(1),'east_d',domain_def(2),'resolution_d',0.5, ...
               'countries',0,'contourf', 0, 'lake',1,'coast',1,'lonlat',0,'frame',1, 'landmask',0);

%% SPEI
[lons lats DATA_SPEI_out]=preproc_SPEI_Global('RS_05', chosen_prd_spei, 'spei06'); % <-- SPEI
% save('/home/minchao/work/GoogleDriver/Research/a1_Studies_and_projects/Legacy_effects/data/Data_SPEI.mat','DATA_SPEI_out');
% load('/home/minchao/work/GoogleDriver/Research/a1_Studies_and_projects/Legacy_effects/data/Data_SPEI.mat');
% save('/home/minchao/work/GoogleDriver/Research/a1_Studies_and_projects/Legacy_effects/data/Data_SPEI_global.mat','DATA_SPEI_out');
% load('./data/Data_SPEI_global.mat');

%% ======= Satellite-based Vegetation Index ======
% preproc_GIMMS_NDVI3g_Global('RS_05'); % <-- NDVI
% preproc_VIP_NDVI_EVI2_Global('RS_05'); % <-- EVI2
preproc_MODIS_NDVI_EVI_Global('RS_05'); % <-- MODIS NDVI and EVI

 % NDVI, EVI: positive values generally indicate the presence of vegetation 
 % (with greater values indicating healthier vegetation). 
 % Negative values generally indicate a lack of vegetation (water, rock, soil, etc.).
 % LAI should > 0
%   DATA_VI_05rs_sm_out=DATA_VI_05rs_out; % No smoothing for monthly data
 
 DATA_VI_05rs_out(DATA_VI_05rs_out<=0)=nan;
 DATA_VI_05rs_sm_out(DATA_VI_05rs_sm_out<=0)=nan;
 DATA_VI_05rs_sm_out(isnan(DATA_VI_05rs_out))=nan; % keep these two array align
 
% use mean to make grid point nan if any of the grid point 
% (spatially and temporally) with nan occurs
global mask_gs;
%                                                              mdi
% mask_gs=squeeze(mean(mean(mean(DATA_VI_05rs_sm_out(:,:,:,:,:,:),3),4),5));
% mask_gs(~isnan(mask_gs))=1;
% save('./data/mask_gs.mat','mask_gs');
% load('./data/mask_gs.mat');
[s1 s2 s3 s4 s5]=size(DATA_VI_05rs_sm_out);
mask_gs=ones(s1,s2); % Skip the masking

%% Climate variable
% Actual evapotranspiration
DATA_ERA5_05rs_out=preproc_ERA5(domain_def, 'RS_05'); % <-- Europe domain
preproc_GLEAM(chosen_prd);
% preproc_GLDAS21(chosen_prd);
preproc_CRU(chosen_prd);
%% Growing season
SOS_POS_EOS_LOS('SAT');
% save('/home/minchao/work/GoogleDriver/Research/a1_Studies_and_projects/Legacy_effects/data/Data_GS.mat', 'sos_map', 'pos_map', 'eos_map', 'lgs_map');
% load('./data/Data_GS.mat');
%% ---------- calculate -------------
% chosen_prd should be within chosen_prd_spei
preproc_Calc('monthly_dsclim',0,1, chosen_prd, chosen_prd_spei,'SAT'); % <-- Paper I:  % 1st param: 1: prewhitened, 0: orig <-- Paper I
                                          % 2nd param: 1: use growing season definition, 0: do not use

%% Calculation
% Identify the extreme year during the GS
% For Satellite dataset
% [ex_y1 ex_y0 ex_pst]=def_extreme(squeeze(DATA_05rs_dtds_out_gs(:,:,:,:,3)), 'SAT');

% Calculate correlation
clearvars DATA_SPEI_out_s; % remove the non-used variable to get more memory
calc_r_spei_evi('SAT');
% save('/home/minchao/work/GoogleDriver/Research/a1_Studies_and_projects/Legacy_effects/data/Data_r_spei_evi.mat', 'r_spei_evi', 'p_spei_evi');
% load('/home/minchao/work/GoogleDriver/Research/a1_Studies_and_projects/Legacy_effects/data/Data_r_spei_evi.mat');
%% Land cover type
preproc_climatezone; % <--
preproc_MODIS_Landcover_Global; % <--
LandClass_ClimteZone;

preproc_ESA_CCI_LC;

figure;
imagesc(rot90(DATA_ESA_CCI_LC_out));
%% Analysis based on annual data
% Scatter plots by land cover and climate zone
% scatter_plots(15);

%% Analysis based on monthly data
% clearvars DATA_SPEI_out;

%% Identify the extreme month during the GS (extreme occurs in EGS)
%  Fig. Sx Sensitivity during the post-extreme
% global ex_pst_m_gs;
% global ex_m_gs;
% nclm=1;
% nsubgs=2;
% alpha_lc=cell(2,nclm,nsubgs);
% beta_lc=cell(2,nclm,nsubgs);
% mdinfo_lc=cell(2,nclm,nsubgs);
% p_slope_diff_lc=cell(2,nclm,nsubgs);
% etype_map=cell(2,1);
% b_diff_m_all_wg=cell(2,nclm,nsubgs);

for subgs=1:1  % the sub-gs in which extremes occur

    [ex_m ex_pst_m norm_m eem_map ex_pst_m_gs etype_map{subgs,1}]=def_extreme_m(squeeze(DATA_05rs_dtds_out_gs(:,:,:,:,3)),1, 'SAT', subgs);

    % identify the extremes occured in 1st or 2nd GS
    ex_m_gs=find_ex_m_gs(ex_m, ex_pst_m_gs);

    for clim=1:nclm
        % 1st GS
        [alpha_lc{1,clim,subgs}, beta_lc{1,clim,subgs}, mdinfo_lc{1,clim,subgs}, p_slope_diff_lc{1,clim,subgs}, b_diff_m_all_wg{1,clim,subgs}]=...
            scatter_plots_m_gsl_clim('SAT',1, clim);
% 
%         % 2nd GS
%         [alpha_lc{2,clim,subgs}, beta_lc{2,clim,subgs}, mdinfo_lc{2,clim,subgs}, p_slope_diff_lc{2,clim,subgs}, b_diff_m_all_wg{2,clim,subgs}]=...
%             scatter_plots_m_gsl_clim('SAT',2, clim);
    end
end

% save('./data/Workspace_2022-01-22_SAT_bf_plotsummary_SPEI03.mat', '*', '-V7.3');
% load('./data/Workspace_2022-01-21_SAT_bf_plotsummary.mat');
% load('./data/Workspace_2022-01-22_SAT_bf_plotsummary_SPEI03.mat');
% load('./data/Workspace_2022-01-18_SAT_bf_plotsummary_SPEI06.mat');

%     export_fig 'scatter_plots_m_gsl_1stGS_SM_EGS.png' -png -r600;
%     export_fig 'scatter_plots_m_gsl_1stGS_RAD_EGS.png' -png -r600;
%     export_fig 'scatter_plots_m_gsl_1stGS_TEMP_EGS.png' -png -r600;
%     export_fig 'scatter_plots_m_gsl_1stGS_SPEI_EGS.png' -png -r600;
%     export_fig 'scatter_plots_m_gsl_2ndGS_SPEI_EGS.png' -png -r600;

%% Fig. 1?
% VI_SPEI_r_plot(etype_map);
VI_SPEI_r_plot_v1(etype_map,1);
% VI_SPEI_r_plot_v1(etype_map,1);

%% Fig. 2: One big summary plot
% Recovery_summary_big(alpha_lc,'SAT',1); % Drought occurs in EGS
Recovery_summary_big_SE(alpha_lc,mdinfo_lc, 'SAT',1); % Drought occurs in EGS, plot SE of alpha also

% Fig. Sx
% Recovery_summary_big(alpha_lc,'SAT',2); % Drought occurs in LGS
Recovery_summary_big_SE(alpha_lc,mdinfo_lc,'SAT',2); % Drought occurs in LGS, plot SE of alpha also

%% Fig. 3 Sensitivity during the post-extreme
% EGS_pst_1stGS=beta_lc{1, 1, 1}; % 1st GS, EGS post-extreme
% EGS_pst_2ndGS=beta_lc{2, 1, 1}; % 2nd GS, EGS post-extreme
% 
% plot_sensitivity(EGS_pst_1stGS);
% plot_sensitivity(EGS_pst_2ndGS);

plot_sensitivity_big(beta_lc, mdinfo_lc,1, b_diff_m_all_wg); % EGS
plot_sensitivity_big(beta_lc, mdinfo_lc,2, b_diff_m_all_wg); % LGS

% plot_sensitivity_big_times(beta_lc, mdinfo_lc);
% plot_sensitivity_big_wg(b_diff_m_all_wg);




%% --> Evaportranspiration index
beta_lc_et=cell(2);
[beta_lc_et{1}]=scatter_plots_m('ET',1);
Recovery_summary(beta_lc_et{1},'ET', 1);

[beta_lc_et{2}]=scatter_plots_m('ET',2);
Recovery_summary(beta_lc_et{2},'ET', 2);

%% Identify the extreme month during the GS (extreme occurs in LGS)
global ex_pst_m_gs;
[ex_m ex_pst_m norm_m eem_map ex_pst_m_gs]=def_extreme_m(squeeze(DATA_05rs_dtds_out_gs(:,:,:,:,3)),1, 'SAT', 2);

global ex_m_gs;
ex_m_gs=test_func(ex_m, ex_pst_m_gs);

[beta_lc]=scatter_plots_m('SAT',1);
Recovery_summary(beta_lc,'SAT', 1);

[beta_lc]=scatter_plots_m('SAT',2);
Recovery_summary(beta_lc,'SAT', 2);

%% ====== FLUXNET data ======
% save('./data/Workspace_2021-07-30.mat', '*');
% load('./data/Workspace_2021-07-30.mat', '*');
load('/home/minchao/work/GoogleDriver/Research/a1_Studies_and_projects/SPAC_in_ESMs/Analysis/data/Fluxnet_output_MM_230sites.mat');

%% Plot the fluxnet site locations
% only get the sites for the northern hemisphere
[s1 s2 s3 s4]=size(Data_EC_MM);
global stn_names_chosen Data_EC_MM_chosen;
global Data_EC_MM_p;

stn_names_chosen=cell(4,1);
Data_EC_MM_chosen=nan(s1,s2,1,s4);
i=0;
for n=1:length(stn_names)
    if str2num(stn_names{3,n}) >= ax_show.south_s
        i=i+1;
        stn_names_chosen(:,i)=stn_names(:,n);
        Data_EC_MM_chosen(:,:,i,:)=Data_EC_MM(:,:,n,:);
    end
end
Data_EC_MM_p=permute(Data_EC_MM_chosen,[3 1 2 4]);
%%
plot_sites(stn_names_chosen);

%% Get the corresponding SPEI
global Data_EC_SPEI_MM;
[Data_EC_SPEI_MM]=get_site_info(stn_names_chosen, DATA_SPEI_out_s);

% Calculate site GS GPP anomalies
global sos_map_s pos_map_s eos_map_s lgs_map_s;
preproc_Calc_site('monthly_dtdsclim',0,1);

% Identify extreme years for sites
[ex_y1_s ex_y0_s ex_pst_s]=def_extreme(DATA_EC_dtds_out_gs(:,:,:,:,18));

[r_spei_gpp p_spei_gpp]=calc_r_spei_gpp(DATA_EC_dtds_out_gs);

%% Analysis based on annual data for sites
global x_p y_p;
global ex_y1_s_p ex_y0_s_p ex_pst_s_p;
select_sites(15);
% Make the scatter plot
plot_extreme_site();

%% Analysis based on monthly data for sites
% Identify extreme month during the GS for sites
[ex_m_s ex_pst_m_s norm_m_s eem_map_s]=def_extreme_m(DATA_EC_dtds_out_gs(:,1,:,:,18),0);

[beta_lc_s]=scatter_plots_m_site();
% Recovery_summary(beta_lc);

%%
% save('/home/minchao/work/GoogleDriver/Research/a1_Studies_and_projects/Legacy_effects/data/Workspace_2021-07-30.mat','*');
% load('/home/minchao/work/GoogleDriver/Research/a1_Studies_and_projects/Legacy_effects/data/Workspace_2021-07-30.mat');

% All site GS-based monthly calculation: 
% save('/home/minchao/work/GoogleDriver/Research/a1_Studies_and_projects/Legacy_effects/data/Workspace_2021-08-19.mat','*');
% load('/home/minchao/work/GoogleDriver/Research/a1_Studies_and_projects/Legacy_effects/data/Workspace_2021-08-19.mat');

% All site GS-based: 
%     extrm_thr=-1.5;
%     norm_thr=-1.0;
% save('/home/minchao/work/GoogleDriver/Research/a1_Studies_and_projects/Legacy_effects/data/Workspace_2021-08-11.mat','*');
% load('/home/minchao/work/GoogleDriver/Research/a1_Studies_and_projects/Legacy_effects/data/Workspace_2021-08-11.mat');

% All site: 
%     extrm_thr=-1.5;
%     norm_thr=-1.5;
% save('/home/minchao/work/GoogleDriver/Research/a1_Studies_and_projects/Legacy_effects/data/Workspace_2021-08-09.mat','*');

% All site: 
%     extrm_thr=-1.5;
%     norm_thr=-1.0;
% save('/home/minchao/work/GoogleDriver/Research/a1_Studies_and_projects/Legacy_effects/data/Workspace_2021-07-30.mat','*');

%%
%% ======= DGVMs dataset ======
% preproc_DGVMs_ISIMIP('RS_05'); % 
preproc_DGVMs_TRENDY('RS_05'); % 

global mask_gs;
[s1 s2 s3 s4 s5]=size(DATA_DGVM_05rs_out);
mask_gs=ones(s1,s2); % Skip the masking

%% Climate variable
% Actual evapotranspiration
% preproc_GLEAM();

preproc_CRU(chosen_prd_dgvm);
%% Growing season
SOS_POS_EOS_LOS('DGVM');

%% ---------- calculate ------------- 
 preproc_Calc('monthly_dsclim',0,1, chosen_prd_dgvm, chosen_prd_spei,'DGVM'); % <-- Paper I:  % 1st param: 1: prewhitened, 0: orig <-- Paper I
%                                                     % 2nd param: 1: use growing season definition, 0: do not use
%% Calculation
% Identify the extreme year during the GS
% For Satellite dataset
% [ex_y1 ex_y0 ex_pst]=def_extreme(squeeze(DATA_05rs_dtds_out_gs(:,:,:,:,3)), 'SAT');
% For DGVMs dataset
% [ex_y1 ex_y0 ex_pst]=def_extreme(squeeze(DATA_05rs_dtds_out_gs(:,:,:,:,2)), 'DGVM');

% Calculate correlation
calc_r_spei_evi('SAT');
% calc_r_spei_evi('DGVM');
%% Analysis based on monthly data
clearvars DATA_SPEI_out;
% save('./data/Workspace_2021-09-07_ORCHIDEE.mat', '*');
% save('./data/Workspace_2021-09-07_SAT.mat', '*');
% load('./data/Workspace_2021-09-07_ORCHIDEE.mat', '*');
% load('./data/Workspace_2021-09-07_LPJ-GUESS.mat', '*');
% load('./data/Workspace_2021-09-07_SAT.mat', '*');

% Identify the extreme month during the GS
global ex_pst_m_gs;
[ex_m ex_pst_m norm_m eem_map]=def_extreme_m(squeeze(DATA_05rs_dtds_out_gs(:,:,:,:,2)),1, 'DGVM');

global ex_m_gs;
ex_m_gs=test_func(ex_m, ex_pst_m_gs);

% [beta_lc]=scatter_plots_m('DGVM');
% Recovery_summary(beta_lc,'DGVM');