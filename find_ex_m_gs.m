function [ex_m_gs]=find_ex_m_gs(ex_m, ex_pst_m)
    % Identify extremes occured during 1st or 2nd GS
    
    [s1 s2 nm ny]=size(ex_m);

    ex_m_gs=nan(s1,s2,nm*ny);

    for i=1:s1
        for j=1:s2
            ex_ts=ex_m(i,j,:,:);
            ex_ts=ex_ts(:);

            ex_pst_ts=ex_pst_m(i,j,:,:);
            ex_pst_ts=ex_pst_ts(:);

            ex_ts(isnan(ex_ts))=0;
            ex_pst_ts(isnan(ex_pst_ts))=0;
            
            found_1stGS_m=0;
            m=2;
            
            while m<=nm*ny % NOTE: loop vaiable can not be changed in a for loop, but while loop does not have such a problem!
                % 1st GS
                break_flag=0;
                if ex_ts(m-1) && ex_pst_ts(m)>0 
                    ex_m_gs(i,j,m)=1; % 1st GS starts here
                    found_1stGS_m=m;
                    
                    for m_gs=(m+1):(m+12)
                        if ex_m_gs(i,j,m_gs-1)==1 && ex_pst_ts(m_gs)>0 % loop within the 1st GS
                            ex_m_gs(i,j,m_gs)=1;  % extend the 1st GS
                        else
                            m=m_gs;
                            break_flag=1;
                            break;
                        end
                    end
                end

                % 2nd GS
                if ex_pst_ts(m)>0 && (m-found_1stGS_m)<12
                    ex_m_gs(i,j,m)=2; % 2nd GS starts here

                    for m_gs=(m+1):(m+12)
                        if ex_m_gs(i,j,m_gs-1)==2 && ex_pst_ts(m_gs)>0 % loop within the 2nd GS
                            ex_m_gs(i,j,m_gs)=2;  % extend the 2nd GS
                        else
                            m=m_gs;
                            break_flag=1;
                            break;
                        end
                    end
                end
                
                if break_flag==0
                    m=m+1;
                end
            end
        end
    end
    ex_m_gs=reshape(ex_m_gs,[s1 s2 nm ny]);
end
