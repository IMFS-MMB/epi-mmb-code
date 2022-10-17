function [out, outstruct] = compute_profits(in)
    
    % "in" is T by 2 (no inner concavity) or T by 3 (with inner concavity)

    global N T q q_m1 rho_r_s rho_r_a w delta_n delta_h delta_s delta_l ...
        A beta alpha1 alpha2 gamma rho_d initial_infected varphi Pi_pq Pi_pn ...
        h_over_n_plus_h n_nd m_nd flag_internalize_q Pi_q flag_planner ...
        delta_d planner_obj_n
    

    % pre-allocates some vectors
    n_vec = zeros(T, 1);
%     m_vec = zeros(T, 1);
%     z_vec = zeros(T, 1);
    
    %---------------------------------------------------------------------%
    %             Distribution of workers across health states            %
    %---------------------------------------------------------------------%    

    % pre-allocates some variables
    if nargout>1
        n_vec = zeros(T, 1);
        l_vec = zeros(T, 1);
        a_vec = zeros(T, 1);
        s_vec = zeros(T, 1);
        r_s_vec = zeros(T, 1);
        r_a_vec = zeros(T, 1);
        d_vec = zeros(T, 1);
        nn_vec = zeros(T, 1);
        mn_vec = zeros(T, 1);
        nm_vec = zeros(T, 1);
        mm_vec = zeros(T, 1);
        nr_vec = zeros(T, 1);
        mr_vec = zeros(T, 1);
        production_vec = zeros(T, 1);
        profit_vec = zeros(T, 1);
        h_vec = zeros(T, 1);
        infected_workplace_vec = zeros(T, 1);
        p_vec = zeros(T, 1);
        n_tilde_vec = zeros(T, 1);
        m_tilde_vec = zeros(T, 1);
    end
    
    % starts profit
    disc_profits = 0;
    disc_planner_obj = 0;
    
    % loop to compute profits
    for t = 1:T
        % computes state variables
        if t==1
            i_tm1 = initial_infected * N; % initial conditions
            i_n_tm1 = i_tm1 * n_nd / (n_nd + m_nd);
            i_m_tm1 = i_tm1 * m_nd / (n_nd + m_nd);
            z_tm1 = 0;
            d_tm1 = 0;
            a_tm1 = 0;
            a_n_tm1 = 0;
            a_m_tm1 = 0;
            r_s_tm1 = 0;
            r_a_tm1 = 0;
            r_a_n_tm1 = 0;
            r_a_m_tm1 = 0;
            nn_tm1 = n_nd; % I think I don't need to store all of this...
            nm_tm1 = 0;
            mn_tm1 = 0;
            mm_tm1 = m_nd;%until here
            q_tm1 = q_m1; % check this (use formula for q)            
            
%             n_tm1 = n_nd; % initial conditions
%             m_tm1 = m_nd;
%             n_tilde_tm1 = (1 - h_over_n_plus_h) * initial_infected;
%             m_tilde_tm1 = h_over_n_plus_h * initial_infected;
%             s_tm1 = 0;
%             a_tm1 = 0;
%             r_s_tm1 = 0;
%             r_a_tm1 = 0;
%             d_tm1 = 0;
        else
            i_tm1 = i_t;
            i_n_tm1 = i_n_t;
            i_m_tm1 = i_m_t;
            z_tm1 = z_t;
            d_tm1 = d_t;
            a_tm1 = a_t;
            a_n_tm1 = a_n_t;
            a_m_tm1 = a_m_t;
            r_s_tm1 = r_s_t;
            r_a_tm1 = r_a_t;
            r_a_n_tm1 = r_a_n_t;
            r_a_m_tm1 = r_a_m_t;
            nn_tm1 = nn_t; % I think I don't need to store all of this...
            nm_tm1 = nm_t;
            mn_tm1 = mn_t;
            mm_tm1 = mm_t; % ... until here
            q_tm1 = q(t-1);               
            
            
            
%             n_tm1 = n_t;
%             m_tm1 = m_t;
%             n_tilde_tm1 = n_tilde_t;
%             m_tilde_tm1 = m_tilde_t;
%             s_tm1 = s_t;
%             a_tm1 = a_t;
%             r_s_tm1 = r_s_t;
%             r_a_tm1 = r_a_t;
%             d_tm1 = d_t;
        end
        
%         x_nn_vec = in(t, 1);
%         x_nm_vec = in(t, 2);        
        
        x_nn_t = in(t, 1);
%         x_nn_t = 0.5 + 0.5 * sin(in(t, 1));
        x_mn_t = 1 - x_nn_t;

        x_nm_t = in(t, 2);
%         x_nm_t = 0.5 + 0.5 * sin(in(t, 2));
        x_mm_t = 1 - x_nm_t;

        x_nr_t = 1;
        x_mr_t = 1 - x_nr_t;

        % choices in levels
        nn_t = x_nn_t * (nn_tm1 + nm_tm1 - varphi * i_n_tm1);
        mn_t = x_mn_t * (nn_tm1 + nm_tm1 - varphi * i_n_tm1);

        nm_t = x_nm_t * (mn_tm1 + mm_tm1 - varphi * i_m_tm1);
        mm_t = x_mm_t * (mn_tm1 + mm_tm1 - varphi * i_m_tm1);
        
        if isnan(nn_t)
            disp('stop here')
        end
        
        if t < T
            r_s_t = r_s_tm1 + (1-rho_d) * rho_r_s * z_tm1;
        else
            r_s_t = r_s_tm1 + z_tm1; % cure
        end

        nr_t = x_nr_t * r_s_t;
        mr_t = x_mr_t * r_s_t;     
        
%         r_s_t = r_s_t_hat - fr_t + h_t * k_r_s_t;

        nr_t = x_nr_t * r_s_t;
        mr_t = x_mr_t * r_s_t;

        % infection probability of on-site worker
        p_tm1 = min(Pi_pq * q_tm1 + Pi_pn * (i_n_tm1 + a_n_tm1), 1);

        % susceptible and exposed in t-1
        s_n_tm1 = nn_tm1 + nm_tm1 - i_n_tm1 - a_n_tm1 - r_a_n_tm1;
        e_n_tm1 = s_n_tm1 * p_tm1;

        if s_n_tm1 < 0
            disp('erro')
        end
        
        s_m_tm1 = mn_tm1 + mm_tm1 - i_m_tm1 - a_m_tm1 - r_a_m_tm1;
        e_m_tm1 = s_m_tm1 * q_tm1;
        
        if s_m_tm1 < 0
            disp('erro')
        end

        % uncertain in t before hiring/firing
        s_hat_n_t = s_n_tm1 * (1-p_tm1);
        i_hat_n_t = e_n_tm1;
        a_hat_n_t = (1-rho_r_a) * a_n_tm1 + (1-varphi) * i_n_tm1;
        r_a_hat_n_t = r_a_n_tm1 + rho_r_a * a_n_tm1;

        s_hat_m_t = s_m_tm1 * (1-q_tm1);
        i_hat_m_t = e_m_tm1;
        a_hat_m_t = (1-rho_r_a) * a_m_tm1 + (1-varphi) * i_m_tm1;
        r_a_hat_m_t = r_a_m_tm1 + rho_r_a * a_m_tm1;
        
        % some measures of workers by health states in t
        d_t = d_tm1 + rho_d * z_tm1;
        r_a_t = r_a_tm1 + rho_r_a * a_tm1;% - f_t_r_a + h_t * k_r_a_t;
        a_t = (1-rho_r_a) * a_tm1 + (1-varphi) * i_tm1;% - f_t_a + h_t * k_a_t;
        
        if t < T
            z_t = (1-rho_d) * (1-rho_r_s) * z_tm1 + varphi * i_tm1;
        else
            z_t = 0; % cure
        end
        
        if z_t < 0
            disp('stop here')
        end

        % incubated in t
        if nn_t + mn_t < 1e-10
            prop_nn_t = 0;
            prop_mn_t = 0;
        else
            prop_nn_t = nn_t / (nn_t + mn_t);
            prop_mn_t = mn_t / (nn_t + mn_t);
        end
        
        if nm_t + mm_t < 1e-10
            prop_nm_t = 0;
            prop_mm_t = 0;
        else
            prop_nm_t = nm_t / (nm_t + mm_t);
            prop_mm_t = mm_t / (nm_t + mm_t);
        end
        
%         if nh_t + mh_t < 1e-10
            prop_nh_t = 0;
            prop_mh_t = 0;
%         else
%             prop_nh_t = nh_t / (nh_t + mh_t);
%             prop_mh_t = mh_t / (nh_t + mh_t);
%         end
        
        i_nn_t = (e_n_tm1) * prop_nn_t;
        i_nm_t = (e_m_tm1) * prop_nm_t;

        i_n_t = i_nn_t + i_nm_t;
        
        if isnan(i_n_t)
            disp('stop here')
        end

        i_mn_t = (e_n_tm1) * prop_mn_t;
        i_mm_t = (e_m_tm1) * prop_mm_t;

        i_m_t = i_mn_t + i_mm_t;

        i_t = i_n_t + i_m_t;
        
        if i_t < 0
            disp('stop here')
        end

        % asymptomatics
%         if nn_t + mn_t < 1e-10
%             prop_nn_t = 0;
%             prop_mn_t = 0;
%         else
%             prop_nn_t = nn_t / (nn_t + mn_t);
%             prop_mn_t = mn_t / (nn_t + mn_t);
%         end
%         
%         if nm_t + mm_t < 1e-10
%             prop_nm_t = 0;
%             prop_mm_t = 0;
%         else
%             prop_nm_t = nm_t / (nm_t + mm_t);
%             prop_mm_t = mm_t / (nm_t + mm_t);
%         end
        
        a_nn_t = (a_hat_n_t) * prop_nn_t;
        a_mn_t = (a_hat_n_t) * prop_mn_t;
        a_nm_t = (a_hat_m_t) * prop_nm_t;
        a_mm_t = (a_hat_m_t) * prop_mm_t;

        a_n_t = a_nn_t + a_nm_t;
        a_m_t = a_mn_t + a_mm_t;

        % recovered asymptomatics
        r_a_nn_t = (r_a_hat_n_t) * prop_nn_t;
        r_a_mn_t = (r_a_hat_n_t) * prop_mn_t;
        r_a_nm_t = (r_a_hat_m_t) * prop_nm_t;
        r_a_mm_t = (r_a_hat_m_t) * prop_mm_t;

        r_a_n_t = r_a_nn_t + r_a_nm_t;
        r_a_m_t = r_a_mn_t + r_a_mm_t;       
        
        n_t = nn_t + nm_t;
        m_t = mn_t + mm_t;        
        
        n_tilde_t = i_nn_t + i_mn_t;
        m_tilde_t = i_nm_t + i_mm_t;
        
        
        % choices
        n = nn_t + nm_t + nr_t;
        m = mn_t + mm_t + mr_t;        
        
        if n < 0
            n = 0;
        end
        
%         % p and q
%         if t==1
%             q_tm1 = q_m1;
%         else
%             q_tm1 = q(t-1);
%         end        
        
%         % some state variables
%         r_s_t = r_s_tm1 + (1 - rho_d) * rho_r_s * s_tm1;
%         r_a_t = r_a_tm1 + rho_r_a * a_tm1;
%         s_t = (1 - rho_d) * (1 - rho_r_s) * s_tm1 + varphi * (n_tilde_tm1 + m_tilde_tm1);
%         a_t = (1 - rho_r_a) * a_tm1 + (1 - varphi) * (n_tilde_tm1 + m_tilde_tm1);
%         d_t = d_tm1 + rho_d * s_tm1;
%         
%         % choices
%         xnn_t = in(t, 1);
%         xnm_t = in(t, 2);
% %         xnn_t = sin(in(t, 1)) * 0.5 + 0.5;
% %         xnm_t = sin(in(t, 2)) * 0.5 + 0.5;
%         
%         nn_t =      xnn_t  * (n_tm1 - varphi * n_tilde_tm1);
%         mn_t = (1 - xnn_t) * (n_tm1 - varphi * n_tilde_tm1);
%         nm_t =      xnm_t  * (m_tm1 - varphi * m_tilde_tm1);
%         mm_t = (1 - xnm_t) * (m_tm1 - varphi * m_tilde_tm1);
%         
%         if alpha1==1
%             xnr_t = 1;
%         else
%             xnr_t = in(t, 3);
% %             xnr_t = sin(in(t, 3)) * 0.5 + 0.5;
%         end
%         
%         nr_t =      xnr_t  * r_s_t;
%         mr_t = (1 - xnr_t) * r_s_t;
%         
%         % p and q
%         if t==1
%             q_tm1 = q_m1;
%         else
%             q_tm1 = q(t-1);
%         end
%         
%         a_n_tm1 = n_tm1 * (a_tm1/(N - r_s_tm1 - s_tm1 - d_tm1));
% %         a_m_tm1 = m_tm1 * (a_tm1/(N - r_s_tm1 - s_tm1 - d_tm1));
%         
%         p_tm1 = Pi_pq * q_tm1 + Pi_pn * (a_n_tm1 + n_tilde_tm1);
%         p_tm1 = min(p_tm1, 1);
%         
%         % other state variables
%         n_t = nn_t + nm_t;
%         m_t = mn_t + mm_t;
%         n_tilde_t = (nn_t * p_tm1 + nm_t * q_tm1) * (1 - (r_a_t + a_t)/(N - r_s_t - s_t - d_t));
%         m_tilde_t = (mn_t * p_tm1 + mm_t * q_tm1) * (1 - (r_a_t + a_t)/(N - r_s_t - s_t - d_t));
        
        % prop. of infected in the workplace
        if n_tilde_t + m_tilde_t > 1e-6
            infected_workplace = (nn_t + mn_t) * (p_tm1 - q_tm1 * Pi_pq) * ...
                (1 - (r_a_t + a_t)/(N - r_s_t - z_t - d_t));
%         if p_tm1 ~= q_tm1
%             zz = 1;
%         end
            infected_workplace = infected_workplace/(n_tilde_t + m_tilde_t);
        else
            infected_workplace = NaN;
        end
%         
%         % choices
%         n = nn_t + nm_t + nr_t;
%         m = mn_t + mm_t + mr_t;
        
        % splitting m among h and l
        if delta_l < delta_h
            if alpha1==1 % no inner concavity
                h = 1/gamma * (gamma * A(t) * alpha2/(w*(delta_h - delta_l)))^(1/(1-alpha2)) - n/gamma;
                l = m - h;

                if h < 0
                    h = 0;
                    l = m;
                elseif h > m
                    h = m;
                    l = 0;
                end
            else % there is inner concavity
                
                if m < 1e-8
                    h = 0;
                    l = 0;
                else
                    
                    int_h = linspace(1e-5, m, 100);
                    trdoff = ((n^alpha1 + gamma * int_h.^alpha1).^(alpha2 - 1)).*(int_h.^(alpha1 - 1)) - ... 
                               ((delta_h - delta_l) * w)/(A(t)*alpha1*alpha2*gamma);
                    if trdoff(end) > 0
                        h = m;
                        l = 0;
                    elseif trdoff(1) < 0
                        h = 0;
                        l = m;     
                    else
                        [min_trdoff, i_min_trdoff] = min(trdoff.^2);
                        h = int_h(i_min_trdoff);
                        l = m - h;     
                        if ~isreal(h) || l < 0 || min_trdoff > 1e-5
                            Trdoff = @(h) ((n^alpha1 + gamma * h^alpha1)^(alpha2 - 1))*(h^(alpha1 - 1)) - ... 
                               ((delta_h - delta_l) * w)/(A(t)*alpha1*alpha2*gamma);
                              h = fzero(Trdoff, [1e-5, m]); %Matlab's solver for h 
                              l = m -  h; 
                          if l < 0 || l > m
                            error('out of bounds')
                          end                         
                        end
                    end
                    
                end
            end
        else
            l = 0;
            h = m;
        end
        
        if n < 0
            n = 0;
        end
        
        % if q is internalized, computes it
        if flag_internalize_q==1
            q(t) = Pi_q * (a_t + n_tilde_t + m_tilde_t);
        end
        
        % production in t
        production_t = A(t) * (n^alpha1 + gamma * h^alpha1)^alpha2;
        
        % profit in t
        profit_t = production_t - ...
            delta_n(t) * w * n - delta_h * w * h - ...
            delta_l * w * l - delta_s * w * z_t;
        
        % discounted profit
        disc_profits = disc_profits + beta^(t-1) * profit_t;
        
        % planner's objective function
        if planner_obj_n==1 || planner_obj_n==2
            planner_obj_t = production_t - delta_d * (d_t - d_tm1);
        elseif planner_obj_n==3
            planner_obj_t = production_t; 
        end
        
        disc_planner_obj = disc_planner_obj + beta^(t-1) * planner_obj_t;
        
        % stores some output variables if asked to
        if nargout>1
            % calculates p_t
            a_n_t = (nn_t + nm_t) * (a_t/(N - r_s_t - z_t - d_t));
            p_t = Pi_pq * q(t) + Pi_pn * (a_n_t + n_tilde_t);
            p_t = min(p_t, 1);
            
            % stores variables in vectors
            n_vec(t) = n;
            h_vec(t) = h;
            l_vec(t) = l;
            a_vec(t) = a_t;
            s_vec(t) = z_t;
            r_s_vec(t) = r_s_t;
            r_a_vec(t) = r_a_t; 
            d_vec(t) = d_t;
            nn_vec(t) = nn_t;
            mn_vec(t) = mn_t;
            nm_vec(t) = nm_t;
            mm_vec(t) = mm_t;
            nr_vec(t) = nr_t;
            mr_vec(t) = mr_t;
            profit_vec(t) = profit_t;
            infected_workplace_vec(t) = infected_workplace;
            p_vec(t) = p_t;
            n_tilde_vec(t) = n_tilde_t;
            m_tilde_vec(t) = m_tilde_t;
            production_vec(t) = production_t;
        end
        
    end
    
    % profit in the no-disease world
    delta_n_nd = 1;
    production_nd = A(T) * (n_nd^alpha1 + gamma * m_nd^alpha1)^alpha2;
    disc_profit_nd =  production_nd - delta_n_nd * w * n_nd - delta_h * w * m_nd;
    disc_profit_nd = disc_profit_nd/(1-beta);
    
    % same for planner
    disc_planner_obj_nd = production_nd/(1-beta);
    
    % variables in t=T+1
    n_Tp1 = (N - d_t - z_t) * n_nd/N;
    m_Tp1 = (N - d_t - z_t) * m_nd/N;
    production_Tp1 = A(T) * (n_Tp1^alpha1 + gamma * m_Tp1^alpha1)^alpha2;
    profit_Tp1 = production_Tp1 - delta_n_nd * w * n_Tp1 - delta_h * w * m_Tp1;
    
    % adds future discounted profits after T
    disc_profits = disc_profits + beta^T * profit_Tp1 / (1-beta);
    
    % same for planner
    if planner_obj_n==1
        disc_planner_obj = disc_planner_obj + beta^T * production_Tp1 / (1-beta);
    elseif planner_obj_n==3
        disc_planner_obj = disc_planner_obj + beta^T * production_Tp1 / (1-beta) - ...
            delta_d * d_t;
    end
    
%     disc_planner_obj = disc_planner_obj + beta^T * profit_Tp1 / (1-beta);
    
%     disc_planner_obj = disc_planner_obj + beta^T * production_Tp1 / (1-beta);
%     disc_planner_obj = delta_d * (- d_t) + (1 - delta_d) * disc_planner_obj;
    
    if flag_planner==0
        out = (disc_profits/disc_profit_nd - 1)*1e+7;
    else
        out = (disc_planner_obj/disc_planner_obj_nd - 1)*1e+7;
    end
    
    % we are looking for a minimum
    out = - out;
    
    % stores output variables in a structdigits(digitsOld) if asked to
    if nargout>1
        % some variables
        incubation = n_tilde_vec + m_tilde_vec;
        exposed = [incubation(2:T); incubation(T)];
        susceptible = N - s_vec - a_vec - r_s_vec - r_a_vec - d_vec - ...
            incubation - exposed;
        
        % saves variables to a struct
        outstruct = struct( ...
            'q', q, ...
            'n', n_vec, ...
            'h', h_vec, ...
            'l', l_vec, ...
            'a', a_vec, ...
            's', s_vec, ...
            'r_s', r_s_vec, ...
            'r_a', r_a_vec, ...
            'd', d_vec, ...
            'nn', nn_vec, ...
            'mn', mn_vec, ...
            'nm', nm_vec, ...
            'mm', mm_vec, ...
            'nr', nr_vec, ...
            'mr', mr_vec, ...
            'profit', profit_vec, ...
            'infected_workplace', infected_workplace_vec, ...
            'p', p_vec, ...
            'n_tilde', n_tilde_vec, ...
            'm_tilde', m_tilde_vec, ...
            'production', production_vec, ...
            'incubation', incubation, ...
            'exposed', exposed, ...
            'susceptible', susceptible, ...
            'A', A, ...
            'disc_profits', disc_profits, ...
            'solution', in);
    end
    
end






