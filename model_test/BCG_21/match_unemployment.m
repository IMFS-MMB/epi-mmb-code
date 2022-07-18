function match_unemployment(share2_add_lockdown)

%march
share_working_while_lockdown(3,1:15) = linspace(0,1,15)*(share1_add_lockdown(1)+share1_lockdown);

share_working_while_lockdown(3,16:75) = 0.4-linspace(0,1,60)*0.075/(l2_ss/(l1_ss+l2_ss));

%april
share_working_while_lockdown(3,16:45) = linspace(0,1,30)*(share1_add_lockdown(2)+share1_lockdown);

%may
share_working_while_lockdown(3,46:75) = linspace(0,1,30)*(share1_add_lockdown(2)+share1_lockdown);

%june
share_working_while_lockdown(3,76:105) = linspace(0,1,30)*(share1_add_lockdown(2)+share1_lockdown);

%july
share_working_while_lockdown(3,106:135) = linspace(0,1,30)*(share1_add_lockdown(2)+share1_lockdown);

%September
share_working_while_lockdown(3,136:165) = share1_lockdown+share1_add_lockdown - linspace(0,1,60)*share1_add_lockdown;

%October
share_working_while_lockdown(3,166:195) = share1_lockdown+share1_add_lockdown - linspace(0,1,60)*share1_add_lockdown;


  

  l1_path_lockdown_day = -max(res_.work(2,:)'-share_working_while_lockdown(2,:)',0).*(Nmat(2,:)'-Dmat(2,:)')...
        -min(res_.work(2,:)',share_working_while_lockdown(2,:)')*(1-share_working_while_infected).*Rmat(2,:)'...
        -(1-res_.work(2,:)')*(1-share_working_while_infected).*Rmat(2,:)'...  
        -Dmat(2,:)';

  l1_path_lockdown = day2month([zeros(15,1);l1_path_lockdown_day(1:end-15)]);

  l2_path_lockdown_day = -max(res_.work(3,:)'-share_working_while_lockdown(3,:)',0).*(Nmat(3,:)'-Dmat(3,:)')...
        -min(res_.work(3,:)',share_working_while_lockdown(3,:)')*(1-share_working_while_infected).*Rmat(3,:)'...
        -(1-res_.work(3,:)')*(1-share_working_while_infected).*Rmat(3,:)'...  
        -Dmat(3,:)';

  l2_path_lockdown = day2month([zeros(15,1);l2_path_lockdown_day(1:end-15)]);

  
  
  
    
    