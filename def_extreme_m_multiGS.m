function [ex_m ex_pst_m norm_m eem_map]=def_extreme_m(data_ref, use_effm, ptype)
    
    global pos_map_s lgs_map_s;
    global pos_map lgs_map;
    global r_spei_evi p_spei_evi;
    global r_spei_gpp p_spei_gpp;
        
    switch ptype
        case 'SAT'
            lgs_map_tmp=lgs_map;
            pos_map_tmp=pos_map;
            r_tmp=r_spei_evi;
            p_tmp=p_spei_evi;
            dsid=2; % EVI
        case 'DGVM'
            lgs_map_tmp=lgs_map;
            pos_map_tmp=pos_map;
            r_tmp=r_spei_evi;
            p_tmp=p_spei_evi;
            dsid=1; %
        case 'FLUXNET'
            fprintf('This is an EC dataset\n');

            lgs_map_tmp=lgs_map_s;
            pos_map_tmp=pos_map_s;
            r_tmp=permute(r_spei_gpp,[1 3 2]);
            p_tmp=permute(p_spei_gpp,[1 3 2]);
            dsid=1;
        otherwise
            error('Must be a 3D or 4D dataset!!');
    end
    
	% Parameters definition
    [s1 s2 nm ny]=size(data_ref);
    extrm_thr=-1.5;
    norm_thr=-1.0;
    
    ex_m=nan(s1,s2,nm*ny);
    ex_pst_m=nan(s1,s2,nm*ny);
	norm_m=nan(s1,s2,nm*ny);
    
    ex_m_ind=((data_ref)<=extrm_thr); % extreme months indicator
    nrm_m_ind=((data_ref)>=norm_thr); % normal month indicator

    eem_map=nan(s1,s2,100);
    
    n=0;
    array_tmp=nan(s1,s2,nm*ny);
    for i=1:s1
        for j=1:s2
            n=n+1;
            if mod(n,1000)==0
                fprintf('def_extreme: %d points processed ...\n',n);
            end
            
            nmon_chk=pos_map_tmp(i,j,dsid);
            if isnan(nmon_chk) || nmon_chk<=0
                continue;
            end
            
            % only check on the focused GS months and set others nan
            atmp=squeeze(ex_m_ind(i,j,:,:));
            atmp=int16(atmp);
            atmp((nmon_chk+1):end,:)=-1;
            ex_ts=reshape(atmp, [1 nm*ny]);
            
            % Time series to mark all the indicators
            ind_ts=nan(1,nm*ny);
            
            spei_ts=reshape(squeeze(data_ref(i,j,:,:)),[1 nm*ny]);
            r_ts=reshape(repmat(squeeze(r_tmp(i,j,1:12)),[1 ny]),[1 nm*ny]);
            p_ts=reshape(repmat(squeeze(p_tmp(i,j,1:12)),[1 ny]),[1 nm*ny]);
            
            % Scan the months
            ex_id=0;
            eemi=0;
            for m=1:nm*ny
                
                % find the Effective Extreme Month 
                % e.g., with effective vegetation response (high SPEI-VI correlation)
                if use_effm==1
                    eff_cond=(r_ts(m) > 0.5) && (p_ts(m) < 0.05); % Used by EVI
                else
                    eff_cond=1;
                end
                
                if (ex_ts(m)==1) && (eff_cond) 

                    ex_id=ex_id-1; % negative values for extreme indicator
                    ind_ts(1,m)=ex_id;
                    
                    % mark down the growing season month ID
                    eemi=eemi+1;
                    gsm=mod(m,12);
                    if gsm==0
                        gsm=12;
                    end
                    eem_map(i,j,eemi)=gsm;
                    
                    % starting from the one month after the extreme months, 
                    % and mark the post-extreme months
                    pst_id=0;
                    for m_pst_ex=(m+1):nm*ny
                        
                        % Only check the GS period
                        if ~isnan(spei_ts(m_pst_ex))
                            % not extreme, mark the post-extreme month ID
                            % (postive values)
                            if (ex_ts(m_pst_ex)~=1)
                                pst_id=pst_id+1;
                                ind_ts(1,m_pst_ex)=pst_id;
                            else % extreme, advance the m, and go back to the outer loop
                                m=m_pst_ex;
                                break;
                            end
                        end
                    end
                end
            end
            
            array_tmp(i,j,:)=ind_ts;
        end
    end
    
    nmon_pstex=12;
    ex_m=(array_tmp<0).*~isnan(array_tmp);
    ex_m=reshape(ex_m,[s1 s2 nm ny]);
    
    % output post-extreme month array 
    ex_pst_m=array_tmp;
    ex_pst_m(ex_pst_m<=0)=nan;
    ex_pst_m(array_tmp>nmon_pstex)=nan;
	ex_pst_m=reshape(ex_pst_m,[s1 s2 nm ny]);
    
    % output normal month array indicators
    norm_m=(array_tmp>nmon_pstex).*~isnan(array_tmp);
	norm_m=reshape(norm_m,[s1 s2 nm ny]);
    
end