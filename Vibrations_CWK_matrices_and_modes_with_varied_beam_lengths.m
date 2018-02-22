%Properties (Moduli - GPA, density - kg/m^3)
E = 200*10^9; G = 75*10^9; rho = 7850;
% Beam A
w_a = 0.048; t_a = 6.30*10^-3; La=0.586;m_a = w_a*t_a*La*rho;
I_a = w_a * t_a^3 / 12; c_a = w_a * t_a^3 / 3;
% Beam B
w_b = 0.048; t_b = 6.41*10^-3; Lb=0.4; m_b = w_b*t_b*Lb*rho;
I_b = w_b * t_b^3 / 12; c_b = w_b * t_b^3 / 3;
% Lengths of both beams
L_a = [La 0.564 0.538];
L_b = [Lb 0.352];
% Mass elements
m_1 = m_a/4; m_2 = m_1; m_3 = m_1; m_4 = m_a/8 + m_b/2; m_5 = m_b/4; 
m_6 = m_5;
% Mass Matrix
m = [m_1 0 0 0 0 0; 0 m_2 0 0 0 0; 0 0 m_3 0 0 0; 0 0 0 m_4 0 0; 
    0 0 0 0 m_5 0; 0 0 0 0 0 m_6];

% Loop through calculations for each length variation
for ii = 1:3
    
    for jj = 1:2
        
        % Set lengths
        l_a = L_a(ii)
        %% i
        l_b = L_b(jj)
    
        % Flexibility Influence Coefficients
        a_11 = (l_a/4)^3/(3*E*I_a); a_21 = ((l_a/4)^2 * (3*l_a/2 - l_a/4))/(6*E*I_a);
        a_12 = a_21; a_31 = ((l_a/4)^2 * (3*3*l_a/4 - l_a/4))/(6*E*I_a); a_13=a_31;
        a_41 =  ((l_a/4)^2 * (3*l_a - l_a/4))/(6*E*I_a); a_14=a_41; a_51 = a_14; 
        a_61 = a_14; a_15 = a_51; a_16 = a_61; a_22 =(l_a/2)^3/(3*E*I_a);
        a_32 = ((l_a/2)^2 * (3*3*l_a/4 - l_a/2))/(6*E*I_a); a_23=a_32;
        a_42 = ((l_a/2)^2 * (3*l_a - l_a/2))/(6*E*I_a); a_24=a_42;
        a_52 = a_42; a_62 = a_42; a_25 = a_52; a_26 = a_62; 
        a_33 =(3*l_a/4)^3/(3*E*I_a); 
        a_43 = ((3*l_a/4)^2 * (3*l_a - 3*l_a/4))/(6*E*I_a);
        a_34 = a_43; a_53 = a_43; a_35 = a_53; a_63 = a_43; a_36 = a_63;
        a_44 = (l_a)^3/(3*E*I_a); a_54 = a_44; a_64 = a_44; a_45 = a_54; 
        a_46 = a_64;
        a_55 = (l_a)^3/(3*E*I_a)+l_a*l_b^2/(4*G*c_a)+l_b^3/(24*E*I_b); a_66 = a_55;
        a_65 = l_a^3/(3*E*I_a)-l_a*l_b^2/(4*G*c_a); a_56 = a_65;

        % Flexibility Influence Matrix 
        a = [a_11 a_12 a_13 a_14 a_15 a_16; a_21 a_22 a_23 a_24 a_25 a_26;
            a_31 a_32 a_33 a_34 a_35 a_36; a_41 a_42 a_43 a_44 a_45 a_46;
            a_51 a_52 a_53 a_54 a_55 a_56; a_61 a_62 a_63 a_64 a_65 a_66]

        % Dynamical Matrix
        D = a*m

        % Solve for eigenvalues/vectors
        [R,L] = eig(D);
        [LS,n] = sort(diag(L));
        R = R(:,n);

        % Natural Frequencies
        tau = sqrt(1./LS)
        f = tau./(2*pi)

        % Define node coordinates
        x1 = [0 0 0 0 0];
        x2 = [0.20 0 -0.20];
        y1 = [0 0.141 0.282 0.423 0.564];
        y2 = [0.564 0.564 0.564];

        % Loop through mode shapes
        for kk=1:6
            max_z = max(abs(R(:,kk)));
            % Normalise amplitudes
            z = R(:,kk).'/max_z;
            z1 = [0 z(1:4)];
            z2 = [z(5) z(4) z(6)];
            % Create reference shape
            z1_stat = [0 0 0 0 0]
            z2_stat = [0 0 0]
            % Plot mode shape
            figure(((ii-1)*2+jj)*6+kk)
            plot3(x1, y1, z1, 'g', x2, y2, z2, 'g')
            hold on
            plot3(x1, y1, z1_stat, 'r', x2, y2, z2_stat, 'r')
            grid on
            title(['Mode Shape ' num2str(kk) ' - La=' num2str(l_a) ', Lb =' num2str(l_b)])
            xlabel('x (m)')
            ylabel('y (m)')
            zlabel('Relative z Displacement')
            hgexport(figure(((ii-1)*2+jj)*6+kk), ['ModeShape' num2str(kk) 'la' num2str(l_a) 'lb' num2str(l_b) '.jpg'], hgexport('factorystyle'), 'Format', 'jpeg');
        end
    end
end
