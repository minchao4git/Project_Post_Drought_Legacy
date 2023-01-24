function [ex_y1 ex_y0 ex_pst]=def_extreme(data_ref, ptype)
    
    global pos_map_s lgs_map_s;
    global pos_map lgs_map;
        
    switch ptype
        case 'SAT'
            lgs_map_tmp=lgs_map;
            pos_map_tmp=pos_map;
            dsid=2; % EVI
        case 'DGVM'
            lgs_map_tmp=lgs_map;
            pos_map_tmp=pos_map;
            dsid=1; % DGVM
        case 'FLUXNET'
            fprintf('This is an EC dataset\n');
        
            lgs_map_tmp=lgs_map_s;
            pos_map_tmp=pos_map_s;
            dsid=1;
        otherwise
            error('Must be a 3D or 4D dataset!!');
    end
    
	% Parameters definition
    [s1 s2 nm ny]=size(data_ref);
    extrm_thr=-1.5;
    norm_thr=-1.0;
    
    % The number of consecutive summer extremes
    n_cex=zeros(s1,s2,ny);
    
    preex_ind=zeros(s1,s2,ny);
    pstex_ind=zeros(s1,s2,ny);
    n_cex_max=nan(s1,s2);
    ind_ex=nan(s1,s2,ny);
    ind_norm=nan(s1,s2,ny);
    
    ex_m_ind=((data_ref)<=extrm_thr); % extreme months indicator
    nrm_m_ind=((data_ref)>=norm_thr); % normal month indicator

    % Calculate the number of consecutive extreme year
    n=0;
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
            
            for y=1:ny
                
                ex_mi=ex_m_ind(i,j,1:nmon_chk,y);
                
                if nansum(ex_mi)>0 
                    ind_ex(i,j,y)=1;

                    if y>0 && nansum(ind_ex(i,j,y))==1
                        n_cex(i,j,y)=1;
                    end
                    if y>1 && nansum(ind_ex(i,j,(y-1:y)))==2
                        n_cex(i,j,y)=2;
                    end
                    if y>2 && nansum(ind_ex(i,j,(y-2:y)))==3
                        n_cex(i,j,y)=3;
                    end
                    if y>3 && nansum(ind_ex(i,j,(y-3:y)))==4
                        n_cex(i,j,y)=4;
                    end
                    if y>4 && nansum(ind_ex(i,j,(y-4:y)))==5
                        n_cex(i,j,y)=5;
                    end
                    if y>5 && nansum(ind_ex(i,j,(y-5:y)))==6
                        n_cex(i,j,y)=6;
                    end
                end

                norm_mi=nrm_m_ind(i,j,1:nmon_chk,y);
                if nansum(norm_mi)==nmon_chk
                    ind_norm(i,j,y)=1;
                    
                    if y>1 && n_cex(i,j,y)==0 && n_cex(i,j,y-1)>=1
                        pstex_ind(i,j,y)=1;
                    end
                end
            end
            ind_norm(isnan(ind_norm))=0;

            % mark the pre-extreme/normal years
            preex_ind(i,j,:)=~(pstex_ind(i,j,:)) & ~(n_cex(i,j,:)) & ind_norm(i,j,:);

            % Maximum number of consecutive extreme year
            n_cex_max(i,j)=nanmax(squeeze(n_cex(i,j,:)));
        end
    end

    % Then for the first year of extreme (non-consecutive extreme year), SPEI vs. EVI
    % FIRST extreme years indicator
    ex_y1=(n_cex==1);
    % month level
    tmp1=repmat(ex_y1,[1 1 1 4]);
    ex_y1_m=permute(tmp1,[1 2 4 3]);

    % normal years indicator
    ex_y0=(preex_ind==1);tmp1=repmat(ex_y0,[1 1 1 4]);
    ex_y0_m=permute(tmp1,[1 2 4 3]);

    % FIRST post-extreme years indicator
    ex_pst=(pstex_ind==1);tmp1=repmat(ex_pst,[1 1 1 4]); 
    ex_pst_m=permute(tmp1,[1 2 4 3]);
    
end