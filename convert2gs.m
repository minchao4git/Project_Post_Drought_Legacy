function [DATA_out_gs]=convert2gs(datain, sos_map, lgs_map)

    [s1 s2 nm ny nv]=size(datain);
    adtmp_rs=reshape(datain,[s1 s2 nm*ny nv]);

    ny_gs=ny-1;
    DATA_out_gs=nan(s1,s2,12,ny_gs,nv);

    % Growing season definition based on calculated SOS
    for var=1:nv
        for i=1:s1
            for j=1:s2

                gssn_s=sos_map(i,j);
                m_lgs=lgs_map(i,j);

                if ~isnan(gssn_s) || m_lgs==12
                    for y=1:ny_gs

                        if ~isnan(gssn_s)
                            m_s=(y-1)*12+gssn_s;
                        else
                            m_s=(y-1)*12+1; % for the case of m_lgs==12, sos use Jan.
                        end
                        m_e=m_s+m_lgs-1;

                        % season 1 growing season dataset
                        DATA_out_gs(i,j,1:m_lgs,y,var)=squeeze(adtmp_rs(i,j,m_s:m_e,var));
                    end % year
                end
            end % lon/lat
        end % lon/lat
    end % data source
end