% Get months for growing season sub-periods
function [m_rng]=get_sub_ssn(m_s, m_p, m_l)

    % For each growing season month, as well as the defined growing period
    m_rng={[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12],[],[] [] [] []};

    if ~isnan(m_s)
        gs_mid=floor(m_p-m_s+1); % the end of the early growing season

        % To deal with the case growing period cross the calendar year
        if gs_mid<0
            gs_mid=gs_mid+12;
        end

        % Dynamical adjust the growing season sub-periods
        m_rng{13}=[1:gs_mid];        % Early GS
        m_rng{14}=[(gs_mid+1):m_l];  % Late GS
        m_rng{15}=[m_rng{13} m_rng{14}];  % all growing season months

        % Middle (middle) of GS
        if m_l >=3 % growing season length have more than 3 months, pick the middle three months

            a=[gs_mid-1 gs_mid gs_mid+1];
            a(a<0)=a(a<0)+12;    % To deal with the case that peak of the month in January
            a(a>12)=a(a>12)-12;  % To deal with the case that peak of the month in December
            m_rng{16}=a;
        else
            m_rng{16}=[gs_mid];
        end

        m_rng{17}=m_l;              % End of GS
    elseif m_l==12
        m_rng{15}=[1:12];  % all calendar months
    end
end