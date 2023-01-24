function []=SOS_POS_EOS_LOS(ptype)
    % Input
    global nyr nds mdi;
    global nyr_dgvm;
    global DATA_VI_05rs_sm_out;
    global DATA_DGVM_05rs_out;
    
    global DATA_CRU_05rs_out;
    global ax_show;
    global mask_gs;
    
    % Output
    global sos_map pos_map eos_map lgs_map;
    
    if strcmp(ptype,'SAT')
        nyr_arry=nyr;
        datain=DATA_VI_05rs_sm_out(:,:,:,:,:,mdi); %
        ndatasrc=nds;
    elseif strcmp(ptype,'DGVM')
        nyr_arry=nyr_dgvm;
        datain=DATA_DGVM_05rs_out; %
        ndatasrc=1;
    end

    [s1 s2 s3 s4 s5]=size(datain);
    dtmp_map_sm=nan(s1,s2,12,nyr_arry,ndatasrc);
	sos_map=nan(s1,s2,ndatasrc);
    pos_map=nan(s1,s2,ndatasrc);
    eos_map=nan(s1,s2,ndatasrc);
    lgs_map=nan(s1,s2,ndatasrc);
    
    mskID=1;
    
    dtmp_map_sm=datain;
    
    for dsid=1:ndatasrc
        
        fprintf(sprintf('==> DSID: %d, calculating SOS ...\n',dsid));
       
        clm_ssn=squeeze(nanmean(dtmp_map_sm(:,:,:,:,dsid),4)); % all years average
        
        if strcmp(ptype,'SAT')
            dtmp_mean_ds=squeeze(nanmean(DATA_VI_05rs_sm_out(:,:,:,:,dsid,mdi),4));
        elseif strcmp(ptype,'DGVM')
            dtmp_mean_ds=squeeze(nanmean(DATA_DGVM_05rs_out,4));
        end
       
        switch nds
            case 1 % NDVI
                clm_ssn_max=1.0;
            case 2 % EVI2
                clm_ssn_max=0.9;
            otherwise
                % Max. values of the entire domain
                clm_ssn_max=squeeze(nanmax(nanmax(nanmax(dtmp_mean_ds))));
                fprintf(sprintf('--- Max. value: %f, \n',clm_ssn_max));
        end
        
        % Max. values of the growing season for each gridpoint
        clm_ssn_max_grid=squeeze(nanmax(dtmp_mean_ds,[],3));
        
        for i=1:s1
            for j=1:s2
                
                assn=squeeze(clm_ssn(i,j,:));
                
                if nansum(isnan(assn))==12 % skip the ocean
                    continue;
                end
                
                if clm_ssn_max_grid(i,j)<=clm_ssn_max*0.2 % skip bare land (% 2021-10-15 revision 0.10->0.15)
                    continue;
                end
                
                atemp=squeeze(nanmean(DATA_CRU_05rs_out(i,j,:,:),4));
                
                [a_mx a_mxi]=max(assn);
                [a_mn a_mni]=min(assn);
                
                assn(isnan(assn))=a_mn; % this is to ensure nan is the smallest, nan is the largest in direct comparision

                if (a_mx-a_mn)<clm_ssn_max_grid(i,j)*0.3
                    lgs_map(i,j,dsid)=12; % growing season all year around

                    sos_map(i,j,dsid)=nan;
                    pos_map(i,j,dsid)=nan;
                    eos_map(i,j,dsid)=nan;
                    continue;
                end
                
                thr_sos=(a_mx-a_mn)*0.1+a_mn;
                thr_eos=(a_mx-a_mn)*0.1+a_mn;

                three_ssn=[assn; assn; assn];
                three_clm=[atemp; atemp; atemp];
                
                mark_3ssn=nan(36,1);
                frs_days=15;

                for m=1:35

                    if three_ssn(m)>=thr_sos && three_ssn(m)<three_ssn(m+1) && three_ssn(m)<=a_mx ... 
                       && three_clm(m)<=frs_days
                   
                       mark_3ssn(m)=1;

                    elseif three_ssn(m)<=a_mx && three_ssn(m)>three_ssn(m+1) && three_ssn(m)>=thr_eos ... 
                           && three_clm(m)<=frs_days
                       
                       mark_3ssn(m)=2;

                    end
                end
                
                % detect SOS
                found_sos=0;m_sos=nan;
                for m=2:36
                    if mark_3ssn(m-1)~=1 && mark_3ssn(m)==1
                        m_sos=m;
                        found_sos=1;
                        break;
                    end
                end
                if ~found_sos
                    fprintf(sprintf('SOS not found!! (i:%d j:%d)\n', i,j));
                    continue;
                end
                
                % detect EOS based on the found SOS month
                % start from m_sos to ensure m_eos >= m_sos
                found_eos=0;m_eos=nan;
                for m=m_sos:35
                    if mark_3ssn(m)==2 && mark_3ssn(m+1)~=2
                        m_eos=m;
                        found_eos=1;
                        break;
                    end
                end
                if ~found_eos
                    fprintf(sprintf('EOS not found!! (i:%d j:%d)\n', i,j));
                    continue;
                end
                
                % to constraint m_eos, data between m_sos and m_eos should
                % not be nan and 
                if sum(isnan(mark_3ssn(m_sos:m_eos)))>0
                    fprintf(sprintf('Error: invalid value between m_sos and m_eos !! (i:%d j:%d)\n', i,j));
                    continue;
                end
                
                % detect POS within the detected growing season
                [dummy mxi]=max(mark_3ssn(m_sos:m_eos));% retruned index starts from 1
                m_pos=m_sos+mxi-1;
                if isnan(m_pos)
                    fprintf(sprintf('POS not found!! (i:%d j:%d)\n', i,j));
                    continue;
                end

                lgs_map(i,j,dsid)=m_eos-m_sos+1;
                
                sos_map(i,j,dsid)=mod(m_sos,12);
                pos_map(i,j,dsid)=mod(m_pos,12);
                eos_map(i,j,dsid)=mod(m_eos,12);
                
                sos_map(sos_map==0)=12;
                pos_map(pos_map==0)=12;
                eos_map(eos_map==0)=12;

            end
        end
    end % dsid
    
    % to be consistent, we always use mdi=4 as the data mask
    
    mask_tmp=squeeze(repmat(squeeze(mask_gs(:,:,mskID)),[1 1 ndatasrc]));
    
    lgs_map(lgs_map==0)=nan;
    
    size(mask_tmp)
    size(sos_map)
    sos_map=sos_map.*mask_tmp;
    pos_map=pos_map.*mask_tmp;
    eos_map=eos_map.*mask_tmp;
    lgs_map=lgs_map.*mask_tmp;
    if strcmp(ptype,'SAT')
        ds_rng={[2]}; % EVI
    elseif strcmp(ptype,'DGVM')
        ds_rng={[1]}; % 
    end
    ax_show.frame=0;
    for dsid=1:length(ds_rng)
        yscal=0.95; xscal=0.02;
        cb_xR=1.0;
        cb_yR=1.0;
        cb_wR=0.4;
        cb_tR=1.0;

        figure('color','w','Position',[371 197  1033 802]);
        gap_h=0.03; gap_w=0.002;
        gap=[gap_h gap_w]; marg_h=[0.08 0.05]; marg_w=[0.1 0.1];
        ha = tight_subplot(4,1,gap,marg_h,marg_w);

        axes(ha(1));
        hold on;
        bg_show=nanmean(sos_map(:,:,ds_rng{dsid}),3)';
        geoplot(ax_show, bg_show);

        title('Start of growing season');
        ctmp=colormap(gca,jet(24));
        ctmp(3:26,:)=ctmp(1:24,:);ctmp(1:2,:)=[1 1 1;1 1 1];
        colormap(gca,ctmp);
        caxis([-0.5 12.5]);

        cb1=colorbar('Eastoutside');
        cb1.Ticks=[0:12];
        cb1.TickLabels={sprintf('No data/seasonality'),'Jan','','Mar','','May','','Jul','','Sep','','Nov',''};
        a=get(gca); gca_w = (a.XLim(2)-a.XLim(1)); gca_h = (a.YLim(2)-a.YLim(1));
        text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,'a)','FontSize',11,'FontName','Arial');
        resizeCB(cb1, cb_xR, cb_yR, cb_wR, cb_tR, '',nan,nan,9);

        axes(ha(2));
        hold on;
        bg_show=nanmean(pos_map(:,:,ds_rng{dsid}),3)';
        geoplot(ax_show, bg_show);

        title('Peak of growing season');
        ctmp=colormap(gca,jet(24));
        ctmp(3:26,:)=ctmp(1:24,:);ctmp(1:2,:)=[1 1 1;1 1 1];
        colormap(gca,ctmp);
        caxis([-0.5 12.5]);

        cb2=colorbar('Eastoutside');
        cb2.Ticks=[0:12];
        cb2.TickLabels={sprintf('No data/seasonality'),'Jan','','Mar','','May','','Jul','','Sep','','Nov',''};
        a=get(gca); gca_w = (a.XLim(2)-a.XLim(1)); gca_h = (a.YLim(2)-a.YLim(1));
        text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,'b)','FontSize',11,'FontName','Arial');
        resizeCB(cb2, cb_xR, cb_yR*0.6, cb_wR, cb_tR*3, '',nan,nan,9);
        
        axes(ha(3));
        hold on;
        bg_show=nanmean(eos_map(:,:,ds_rng{dsid}),3)';
        geoplot(ax_show, bg_show);

        title('End of growing season');
        ctmp=colormap(gca,jet(24));
        ctmp(3:26,:)=ctmp(1:24,:);ctmp(1:2,:)=[1 1 1;1 1 1];
        colormap(gca,ctmp);
        caxis([-0.5 12.5]);

        cb2=colorbar('Eastoutside');
        cb2.Ticks=[0:12];
        cb2.TickLabels={sprintf('No data/seasonality'),'Jan','','Mar','','May','','Jul','','Sep','','Nov',''};
        a=get(gca); gca_w = (a.XLim(2)-a.XLim(1)); gca_h = (a.YLim(2)-a.YLim(1));
        text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,'c)','FontSize',11,'FontName','Arial');
        resizeCB(cb2, cb_xR, cb_yR, cb_wR, cb_tR*2, '',nan,nan,9);

        axes(ha(4));
        hold on;
        bg_show=nanmean(lgs_map(:,:,ds_rng{dsid}),3)';
        geoplot(ax_show, bg_show);

        title('Length of growing season');
        a=jet(30);
        ctmp=a(1:24,:);
        ctmp(3:26,:)=ctmp(1:24,:);ctmp(1:2,:)=[1 1 1;1 1 1];
        colormap(gca,ctmp);
        caxis([-0.5 12.5]);

        cb3=colorbar('Eastoutside');
        cb3.Ticks=[0:12];
        cb3.TickLabels={'No data','1','','3','','5','','7','','9','','11',''};
        a=get(gca); gca_w = (a.XLim(2)-a.XLim(1)); gca_h = (a.YLim(2)-a.YLim(1));
        text(a.XLim(1)+gca_w*xscal, a.YLim(1)+gca_h*yscal,'d)','FontSize',11,'FontName','Arial');
        resizeCB(cb3, cb_xR, cb_yR, cb_wR, cb_tR, '',5,2,8);
    end
end
