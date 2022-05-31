%UScontacts

%% pullting in data:

UScontactsother
contacts_.other = (Contact_other_agg_3+Contact_other_agg_3')/2;
UScontactsathome
contacts_.home = (Contact_home_agg_3+Contact_home_agg_3')/2;
UScontactsatwork
contacts_.work = (Contact_work_agg_3+Contact_work_agg_3')/2;
UScontactsatschool
contacts_.school = (Contact_school_agg_3+Contact_school_agg_3')/2;

% cleaning up contact matrices:
cleaning = 1;
if cleaning == 1
    aux1 = contacts_.school(1,1);  %about 7.0
    contacts_.school = 0*contacts_.school;  
    contacts_.school(1,1) = aux1; 
    aux1 = contacts_.work(2,2);  %about 5.4
    contacts_.work = 0*contacts_.work;  
    contacts_.work(2,2) = aux1; 
end

%% contacts_.other = (Contact_other_agg_3+Contact_other_agg_3')/2;
contacts3_.other = contacts_.other;
contacts3_.home = contacts_.home;
contacts3_.school = contacts_.school;
contacts3_.work = contacts_.work;

Ninit_agg_3 = [sum(popsharebyage(1:4));...
               sum(popsharebyage(5:13));...
               sum(popsharebyage(14:16))];

%% four population groups (two within the middle-aged)
% We split the middle-aged into two groups to distinguish the sector of
% employment

% other:
contacts4_.other(1,1) = contacts_.other(1,1);
contacts4_.other(1,4) = contacts_.other(1,3);
contacts4_.other(4,1) = contacts_.other(3,1);
contacts4_.other(4,4) = contacts_.other(3,3);

contacts4_.other(1,2) = contacts_.other(1,2)/2;
contacts4_.other(1,3) = contacts_.other(1,2)/2;
contacts4_.other(4,2) = contacts_.other(3,2)/2;
contacts4_.other(4,3) = contacts_.other(3,2)/2;

contacts4_.other(2,1) = contacts_.other(2,1);
contacts4_.other(3,1) = contacts_.other(2,1);
contacts4_.other(2,4) = contacts_.other(2,3);
contacts4_.other(3,4) = contacts_.other(2,3);

contacts4_.other(2,2) = contacts_.other(2,2)/2; 
contacts4_.other(2,3) = contacts_.other(2,2)/2; 
contacts4_.other(3,2) = contacts_.other(2,2)/2; 
contacts4_.other(3,3) = contacts_.other(2,2)/2; 

% home:
contacts4_.home(1,1) = contacts_.home(1,1);
contacts4_.home(1,4) = contacts_.home(1,3);
contacts4_.home(4,1) = contacts_.home(3,1);
contacts4_.home(4,4) = contacts_.home(3,3);

contacts4_.home(1,2) = contacts_.home(1,2)/2;
contacts4_.home(1,3) = contacts_.home(1,2)/2;
contacts4_.home(4,2) = contacts_.home(3,2)/2;
contacts4_.home(4,3) = contacts_.home(3,2)/2;

contacts4_.home(2,1) = contacts_.home(2,1);
contacts4_.home(3,1) = contacts_.home(2,1);
contacts4_.home(2,4) = contacts_.home(2,3);
contacts4_.home(3,4) = contacts_.home(2,3);

contacts4_.home(2,2) = contacts_.home(2,2)/2; 
contacts4_.home(2,3) = contacts_.home(2,2)/2; 
contacts4_.home(3,2) = contacts_.home(2,2)/2; 
contacts4_.home(3,3) = contacts_.home(2,2)/2; 

% school:
contacts4_.school(1,1) = contacts_.school(1,1);
contacts4_.school(1,4) = contacts_.school(1,3);
contacts4_.school(4,1) = contacts_.school(3,1);
contacts4_.school(4,4) = contacts_.school(3,3);

contacts4_.school(1,2) = contacts_.school(1,2)/2;
contacts4_.school(1,3) = contacts_.school(1,2)/2;
contacts4_.school(4,2) = contacts_.school(3,2)/2;
contacts4_.school(4,3) = contacts_.school(3,2)/2;

contacts4_.school(2,1) = contacts_.school(2,1);
contacts4_.school(3,1) = contacts_.school(2,1);
contacts4_.school(2,4) = contacts_.school(2,3);
contacts4_.school(3,4) = contacts_.school(2,3);

contacts4_.school(2,2) = contacts_.school(2,2)/2; 
contacts4_.school(2,3) = contacts_.school(2,2)/2; 
contacts4_.school(3,2) = contacts_.school(2,2)/2; 
contacts4_.school(3,3) = contacts_.school(2,2)/2;

% work:
contacts4_.work(1,1) = contacts_.work(1,1);
contacts4_.work(1,4) = contacts_.work(1,3);
contacts4_.work(4,1) = contacts_.work(3,1);
contacts4_.work(4,4) = contacts_.work(3,3);

contacts4_.work(1,2) = contacts_.work(1,2)/2;
contacts4_.work(1,3) = contacts_.work(1,2)/2;
contacts4_.work(4,2) = contacts_.work(3,2)/2;
contacts4_.work(4,3) = contacts_.work(3,2)/2;

contacts4_.work(2,1) = contacts_.work(2,1);
contacts4_.work(3,1) = contacts_.work(2,1);
contacts4_.work(2,4) = contacts_.work(2,3);
contacts4_.work(3,4) = contacts_.work(2,3);

contacts4_.work(2,2) = contacts_.work(2,2)/2; 
contacts4_.work(2,3) = contacts_.work(2,2)/2; 
contacts4_.work(3,2) = contacts_.work(2,2)/2; 
contacts4_.work(3,3) = contacts_.work(2,2)/2;

Ninit_agg_4 = [sum(popsharebyage(1:4));...
               sum(popsharebyage(5:13))*(0.25/(0.25+0.40));...
               sum(popsharebyage(5:13))*(0.40/(0.25+0.40));...
               sum(popsharebyage(14:16))];

%% aggregate population groups 
% 
contacts1_.other = Contact_other_agg_1;
contacts1_.home = Contact_home_agg_1;
contacts1_.school = Contact_school_agg_1;
contacts1_.work = Contact_work_agg_1;

Ninit_agg_1 = 1;

% %reaggregaging
% Ninit_agg_4 = [Ninit_agg_3(1,1);Ninit_agg_3(2,1)/2;Ninit_agg_3(2,1)/2;Ninit_agg_3(3,1)]';
% Contact_home_agg_3 = ...
% [sum(sum(contacts4_.other(1:1,1:1),2).*Ninit_agg_4(1:1)')/sum(Ninit_agg_4(1:1)),...
%  sum(sum(contacts4_.other(1:1,2:3),2).*Ninit_agg_4(1:1)')/sum(Ninit_agg_4(1:1)),...
%  sum(sum(contacts4_.other(1:1,4:4),2).*Ninit_agg_4(1:1)')/sum(Ninit_agg_4(1:1));...
%  sum(sum(contacts4_.other(2:3,1:1),2).*Ninit_agg_4(2:3)')/sum(Ninit_agg_4(2:3)),...
%  sum(sum(contacts4_.other(2:3,2:3),2).*Ninit_agg_4(2:3)')/sum(Ninit_agg_4(2:3)),...
%  sum(sum(contacts4_.other(2:3,4:4),2).*Ninit_agg_4(2:3)')/sum(Ninit_agg_4(2:3));...
%  sum(sum(contacts4_.other(4:4,1:1),2).*Ninit_agg_4(4:4)')/sum(Ninit_agg_4(4:4)),...
%  sum(sum(contacts4_.other(4:4,2:3),2).*Ninit_agg_4(4:4)')/sum(Ninit_agg_4(4:4)),...
%  sum(sum(contacts4_.other(4:4,4:4),2).*Ninit_agg_4(4:4)')/sum(Ninit_agg_4(4:4))];

clear contacts_


