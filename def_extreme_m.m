function [ex_m ex_pst_m norm_m eem_map ex_pst_m_gs etype_map]=def_extreme_m(data_ref, use_effm, ptype, subgs)
%     data_ref=DATA_05rs_dtds_out_gs(:,:,:,:,3);
    % Input variables
    % data_ref: extreme index dataset (gs version)
    % must be 4D array: 
    % d1 and d2, spatial or site/varaibles dimension
    % d3: month, d4: years
    
    global sos_map_s pos_map_s lgs_map_s;
    global sos_map pos_map lgs_map;
    global r_spei_evi p_spei_evi;
    global r_spei_gpp p_spei_gpp;
    global DATA_05rs_dtds_out_gs;
    
        
    switch ptype
        case 'SAT'
            lgs_map_tmp=lgs_map;
            sos_map_tmp=sos_map;
            pos_map_tmp=pos_map;
            r_tmp=r_spei_evi;
            p_tmp=p_spei_evi;
            dsid=2; % use EVI gs definition
        case 'DGVM'
            lgs_map_tmp=lgs_map;
            sos_map_tmp=sos_map;
            pos_map_tmp=pos_map;
            r_tmp=r_spei_evi;
            p_tmp=p_spei_evi;
            dsid=1; %
        case 'FLUXNET'
            fprintf('This is an EC dataset\n');

            % Add one more dimension to align with the VI dataset dimension
            lgs_map_tmp=lgs_map_s;
            sos_map_tmp=sos_map_s;
            pos_map_tmp=pos_map_s;
            r_tmp=permute(r_spei_gpp,[1 3 2]);
            p_tmp=permute(p_spei_gpp,[1 3 2]);
            dsid=1;
        otherwise
            error('Must be a 3D or 4D dataset!!');
    end
    
	% Parameters definition
    [s1 s2 nm ny]=size(data_ref);
    extrm_thr=-1.2;
    norm_thr=-1.0;
    
    ex_m_tmp=nan(s1,s2,nm*ny);
    ex_pst_m=nan(s1,s2,nm*ny);
	norm_m=nan(s1,s2,nm*ny);
    
    ex_m_ind=((data_ref)<=extrm_thr); % extreme months indicator
    nrm_m_ind=((data_ref)>=norm_thr); % normal month indicator

    % VI anomalies
    vi_m_amom=squeeze(DATA_05rs_dtds_out_gs(:,:,:,:,2));
    
    % map for extreme timing in growing season months
    eem_map=nan(s1,s2,100);
    
    % map for the types of droughts
    etype_map=nan(s1,s2,2); % dimension 3: index 1 for negative impacts, index 2 for no or postive impacts
    
    n=0;
    array_tmp=nan(s1,s2,nm*ny);
    array_tmp1=nan(s1,s2,nm*ny);
    
    % Markdown the effective month
    eff_m_tmp=nan(s1,s2,nm*ny);
    for i=1:s1
        for j=1:s2
            n=n+1;
            if mod(n,1000)==0
                fprintf('def_extreme: %d points processed ...\n',n);
            end
            
            nmon_chk=pos_map_tmp(i,j,dsid)-sos_map_tmp(i,j,dsid)+1;
            if isnan(nmon_chk) || nmon_chk<=0
                continue;
            end
            atmp=squeeze(ex_m_ind(i,j,:,:));
            vitmp=squeeze(vi_m_amom(i,j,:,:));

            % only check on the focused GS months and set others nan
            if subgs==1
                % for EGS
                atmp((nmon_chk+2):end,:)=false;
            elseif subgs==2
                % for LGS
                atmp(1:(nmon_chk),:)=false;
            end
            
            ex_ts=reshape(atmp, [1 nm*ny]);
            vi_ts=reshape(vitmp, [1 nm*ny]);
            
            % Time series to mark all the indicators
            ind_ts=nan(1,nm*ny);
            ind_ts_gs=nan(1,nm*ny);
            
            spei_ts=reshape(squeeze(data_ref(i,j,:,:)),[1 nm*ny]);
            r_ts=reshape(repmat(squeeze(r_tmp(i,j,1:12)),[1 ny]),[1 nm*ny]);
            p_ts=reshape(repmat(squeeze(p_tmp(i,j,1:12)),[1 ny]),[1 nm*ny]);
            
            % Scan the months
            ex_id=0;
            eemi=0;
            
            m=1;
            while m<=nm*ny % NOTE: loop vaiable can not be changed in a for loop, but while loop does not have such a problem!
                
                if use_effm==1
                    eff_cond=(r_ts(m) > 0.5) && (p_ts(m) < 0.05); % Used by EVI
                    eff_m_tmp(i,j,m)=1;
                % otherwise, we regard all months as effective months
                else
                    eff_cond=1;
                    eff_m_tmp(i,j,m)=1;
                end
                
                break_flag=0;
                if (ex_ts(m)) && (eff_cond)

                    if vi_ts(m) < 0
                        etype_map(i,j,1)=1;
                    end
                    
                    if vi_ts(m) > 0
                        etype_map(i,j,2)=1;
                    end
                    
                    if vi_ts(m) < 0
                        
                    ex_id=ex_id-1; % negative values for extreme indicator
                    ind_ts(1,m)=ex_id;
                    
                    % mark down the growing season month ID
                    eemi=eemi+1;
                    gsm=mod(m,12);
                    if gsm==0
                        gsm=12;
                    end
                    eem_map(i,j,eemi)=gsm;
                    
                    pst_id=0;    % extreme event level post-extreme ID
                    pst_gs_id=0; % GS level post-extreme ID
                    for m_pst_ex=(m+1):nm*ny
                        
                        % Only check the GS period
                        if ~isnan(spei_ts(m_pst_ex))
                            % not extreme, mark the post-extreme month ID
                            % (postive values)
                            if (ex_ts(m_pst_ex)==false)
                                pst_id=pst_id+1;
                                pst_gs_id=pst_gs_id+1;
                                
                                ind_ts(1,m_pst_ex)=pst_id;
                                ind_ts_gs(1,m_pst_ex)=pst_gs_id;
                                
                            else % extreme, advance the m, and go back to the outer loop
                                m=m_pst_ex;
                                break_flag=1;
                                break;
                            end
                        else
                            pst_gs_id=0;
                        end
                    end
                    
                    end % skip the non-effective extreme
                end
                
                if break_flag==0
                    m=m+1;
                end
            end
            
            array_tmp(i,j,:)=ind_ts;
            array_tmp1(i,j,:)=ind_ts_gs;
        end
    end
    
    % This threshold would affect the definition of normal month and
    % post-extreme months
    nmon_pstex=12;

    % output extreme month array indicator
    ex_m_tmp=(array_tmp<0).*eff_m_tmp; % we need this eff_m_tmp to do spatial and temporal filtering (get the signficant r point)
    ex_m=reshape(ex_m_tmp,[s1 s2 nm ny]);
    
    % output post-extreme month array
    % (NOT indicator, need to keep the post-extreme ID for later analysis)
    ex_pst_m=array_tmp; % for post-extreme months, we don't need to check effective month
    ex_pst_m(ex_pst_m<=0)=nan;
    ex_pst_m(array_tmp>nmon_pstex)=nan;
	ex_pst_m=reshape(ex_pst_m,[s1 s2 nm ny]);
    
    % Dont know what the below used for
    % GS level post-extreme ID
    ex_pst_m_gs=array_tmp1;
    ex_pst_m_gs(ex_pst_m_gs<=0)=nan;
    ex_pst_m_gs(array_tmp1>nmon_pstex)=nan;
	ex_pst_m_gs=reshape(ex_pst_m_gs,[s1 s2 nm ny]);
    
    % output normal month array indicators
    % i.e. : not extreme and post-extreme months
    %     1.post-extreme normal   2.non-extreme and non-post-extreme  3.
    %     skimed by effective month
    norm_m=((array_tmp>nmon_pstex)|(isnan(array_tmp))).*eff_m_tmp;
	norm_m=reshape(norm_m,[s1 s2 nm ny]);
    
end