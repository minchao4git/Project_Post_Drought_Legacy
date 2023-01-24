function calc_r_spei_evi(ptype)
    % Check the point-based correlation coefficent bewtween SPEI and EVI anomalies

    % Input
    global DATA_05rs_dtds_out_gs DATA_05rs_dtds_out_sgs;
    global sos_map pos_map lgs_map;
    
    % Output
	global r_spei_evi p_spei_evi;
    global r_spei_evi_asyn p_spei_evi_asyn;
        
    switch ptype
        case 'SAT'
            dsid=2; % using EVI
        case 'DGVM'
            dsid=1; %
        otherwise
            error('Wrong dataset type!!');
    end
    spei_id=3;
    
    [s1 s2 nm ny nv]=size(DATA_05rs_dtds_out_gs);
    DATA_05rs_dtds_out_sgs=nan(s1,s2,17,ny,nv);
    
    r_spei_evi=nan(s1,s2,12);
    p_spei_evi=nan(s1,s2,12);
    r_spei_evi_asyn=nan(s1,s2);
    p_spei_evi_asyn=nan(s1,s2);
    
    n=0;
    for i=1:s1
        for j=1:s2
            n=n+1;
            if mod(n,1000)==0
                fprintf('calc_r_spei_evi: %d points processed ...\n',n);
            end

            m_s=sos_map(i,j,dsid);
            m_p=pos_map(i,j,dsid);
            m_l=lgs_map(i,j,dsid);

            if isnan(m_l)
                continue;
            end

            [m_rng]=get_sub_ssn(m_s, m_p, m_l); % also get the end of growing season

            % Synchronous correlation
            for m=1:length(m_rng)
                
                if isnan(m_rng{m})
                    continue;
                end
                
                spei_tmp=squeeze(nanmean(DATA_05rs_dtds_out_gs(i,j,m_rng{m},:,spei_id),3));
                evi_tmp=squeeze(nanmean(DATA_05rs_dtds_out_gs(i,j,m_rng{m},:,dsid),3));
                
                [r p]=corrcoef(spei_tmp,evi_tmp);
                r_spei_evi(i,j,m)=r(1,2);
                p_spei_evi(i,j,m)=p(1,2);
                
                % Sub-growing season mean
                DATA_05rs_dtds_out_sgs(i,j,m,:,:)=squeeze(nanmean(DATA_05rs_dtds_out_gs(i,j,m_rng{m},:,:),3));
                
            end
            
            % Asynchronous correlation
            if ~isnan(m_rng{13})
                spei_asyn=squeeze(nanmean(DATA_05rs_dtds_out_sgs(i,j,13,:,spei_id),3));   % EGS
                evi_asyn=squeeze(nanmean(DATA_05rs_dtds_out_sgs(i,j,14,:,dsid),3)); % LGS

                [r p]=corrcoef(spei_asyn,evi_asyn);
                r_spei_evi_asyn(i,j)=r(1,2);
                p_spei_evi_asyn(i,j)=p(1,2);
            end
        end
    end
end