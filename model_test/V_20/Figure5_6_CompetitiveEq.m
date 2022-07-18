% Code to reproduce figures 5 and 6 in Policy During an Epidemic With
% Super-Spreaders, Van Vlokhoven (2020)
% Evolution of epidemic in the competitive equilibirum


%clear all

%% parameters
T = 100;    %(1 period is a week), so solve for T weeks
dt = 1/2;  %time steps

beta = 0.96^(dt/52);     %dicounting
  
phi = [0.07 0.2];       %weight on social consumption in utility for each type of agent
nphi = length(phi);     %number of type of agents
g = [0.7 0.3];          %density phi
g = g/sum(g);
E_phi = phi*g';         %expectation phi

%epsilon = 0.02;     %fraction infected initially       // initial epi shock
helper=load('inf_ini.mat');
epsilon=helper.helper;

alpha = 0.7;        %returns to scale production function
As = (E_phi/(1-E_phi))^(1-alpha);   %productivity social sector: set A_s such that P_s is 1 in economy with no infections (Ar=L=1)

%disease parameters
p = 0.07;

lambda_R = 0.4975*dt;  %recovery rate
lambda_D = 0.0025*dt;  %death rate

% create grids for disease states
dS = 0.025;
S_grid = 0:dS:1;     
dI = 0.03;   
I_grid = 0:dI:0.5; 

T_grid = 1:dt:T;

nS = length(S_grid);    %make sure nS  = nR
nI = length(I_grid);
nT = length(T_grid);

xs_grid = 0.3:(0.7/300):1.05;   %grid for share of phi spent on social good

[I_mat,S_mat,XS_mat] = meshgrid(I_grid,S_grid,xs_grid);
I_tmp = I_mat(:,:,1);
S_tmp = S_mat(:,:,1);

alive = min(1-(1-S_mat - I_mat)*lambda_D/(lambda_D+lambda_R),1);    %mass of agents alive depending on state
Itilde_mat = I_mat./alive;  %fraction of alive agents that are infected

%% solve steady state (pre-disease onset) and value of life
B_ss = ((As).^(alpha/(alpha-1))+As)./((1+(As).^(1/(alpha-1))).^alpha);  %Ps=1,L=1,Ar=1
cr_ss = (1-phi)*B_ss;
cs_ss = phi*B_ss;
CS_ss = cs_ss*g';

% set parameters to get following degree distribution in steady state
degree_target = [10 25]*dt;
d_CS = 1; %parameter that governs to what extent degree depends on aggregate consumption (in paper set to 1)
a = (degree_target(2)-degree_target(1))/((CS_ss.^(d_CS))*B_ss*(phi(2)-phi(1)));
b = degree_target(1) - a*(CS_ss.^(d_CS))*cs_ss(1);
%degree = b + a * C^d_CS * c
degree_ss = b+a*cs_ss*(CS_ss.^(d_CS));
degree_ss_exp = degree_ss*g';

value_life = 16*(52/dt)*B_ss;    %B_ss is income per capita per week. 10 qaly (= 1 million dollar) corresponds to 16 times GDP per capita

u_lowerbar = (1-beta)*value_life/B_ss - log(B_ss);

u_ss = (1-phi).*log(cr_ss./(1-phi)) + phi.*log(cs_ss./phi) + u_lowerbar;

V_ss = u_ss/(1-beta); %value stady state when u_lowerbar = 0

%% value function iteration during epidemic
%for simulations
S_phi = NaN(nT,nphi);
S_phi(1,:) = (1-epsilon*degree_ss/degree_ss_exp); %reflects that those with a higher degree are more likely to be infected initially
I_phi = NaN(nT,nphi);
I_phi(1,:) = epsilon*degree_ss/degree_ss_exp;
cs_eq = NaN(nT,nphi);   %consumption of social good over time for each type
eps1 = 1;
iter1 = 0;
gamma = 1;    %update rate time variables (homothopy)
gammaP = 0.2;  

% initial guess (taken from SIR model)
%f is probability that a given social interaction is with an infectious agent
f = [0.0244946492271106 0.0331573335731152 0.0445267930518307 0.0591588599698761 0.0774947770486490 0.0996601619496799 0.125190369922288 0.152751650910325 0.180016907078723 0.203904392937005 0.221278107209198 0.229900038116812 0.229136161388631 0.220003342443696 0.204602552845394 0.185355149837881 0.164426522692483 0.143459523143270 0.123545356922253 0.105313281878629 0.0890530632515690 0.0748272008879687 0.0625584940648112 0.0520923990632809 0.0432386757376051 0.0357976867771286 0.0295760059334871 0.0243949019980018 0.0200942545652154 0.0165336574862144 0.0135918790738595 0.0111654364031940 0.00916676062503583 0.00752224348413933 0.00617033336002323 0.00505977121706789 0.00414800819810430 0.00339981716273339 0.00278609346377806 0.00228283114134780 0.00187025650627964 0.00153209983594102 0.00125498636283789 0.00102792910705062 0.000841907895504827 0.000689520828446821 0.000564696325754989 0.000462455621730230 0.000378717134732643 0.000310135504962639 0.000253969274310280 0.000207972190082396 0.000170303967312078 0.000139457061250440 0.000114196601077970 9.35111350421904e-05 7.65722515517821e-05 6.27014837662583e-05 5.13431886005571e-05 4.20423247906674e-05 3.44262471801887e-05 2.81897927827011e-05 2.30830643825485e-05 1.89014243962592e-05 1.54772995278410e-05 1.26734688073544e-05 1.03775667072133e-05 8.49758149702453e-06 6.95816873007169e-06 5.69763231986370e-06 4.66545235166833e-06 3.82026063982203e-06 3.12818295556864e-06 2.56148152509113e-06 2.09744341798946e-06 1.71747029220148e-06 1.40633302533186e-06 1.15156136671343e-06 9.42944153282516e-07 7.72120061822058e-07 6.32242497533874e-07 5.17705189421250e-07 4.23917495513319e-07 3.47120412945607e-07 2.84235919102396e-07 2.32743605761984e-07 1.90579661973392e-07 1.56054157052692e-07 1.27783308504088e-07 1.04634020232112e-07 8.56784681815034e-08 7.01569132269259e-08 5.74472508726273e-08 4.70400773267508e-08 3.85182726068159e-08 3.15402823725570e-08 2.58264284253550e-08 2.11476992011992e-08 1.73165709814269e-08 1.41794919361374e-08 1.16107277611685e-08 9.50732224198809e-09 7.78497076399378e-09 6.37464137952888e-09 5.21980800395560e-09 4.27418484723523e-09 3.49987127697576e-09 2.86583275899908e-09 2.34665699146929e-09 1.92153537838655e-09 1.57342901985826e-09 1.28838579199784e-09 1.05498114497908e-09 8.63860206402860e-10 7.07362837454764e-10 5.79216614083447e-10 4.74285427868389e-10 3.88363630485361e-10 3.18007471067930e-10 2.60397070464610e-10 2.13223400311782e-10 1.74595737036933e-10 1.42965881543139e-10 1.17066107295249e-10 9.58583497633035e-11 7.84926007330567e-11 6.42728399250317e-11 5.26291384593596e-11 4.30948160715754e-11 3.52877365391979e-11 2.88949916386566e-11 2.36603597646274e-11 1.93740365524796e-11 1.58642259065564e-11 1.29902544021908e-11 1.06369330862709e-11 8.70994069697832e-12 7.13204326186486e-12 5.83999855552864e-12 4.78202134736472e-12 3.91570784636544e-12 3.20633615459171e-12 2.62547461138673e-12 2.14984225068352e-12 1.76037569846554e-12 1.44146511157389e-12 1.18032853424171e-12 9.66499596527833e-13 7.91408021571326e-13 6.48036128372433e-13 5.30637562710211e-13 4.34506982915003e-13 3.55791469487430e-13 2.91336099849943e-13 2.38557498857553e-13 1.95340296964513e-13 1.59952346083945e-13 1.30975294986911e-13 1.07247741698677e-13 8.78186844367515e-14 7.19094054014648e-14 5.88822597190739e-14 4.82151185963476e-14 3.94804423666983e-14 3.23281446742777e-14 2.64715609914892e-14 2.16759590872435e-14 1.77491309448246e-14 1.45336890528608e-14 1.19007583042728e-14 9.74481067412403e-15 7.97943564994737e-15 6.53387689313672e-15 5.35019632058158e-15 4.38095194276347e-15 3.58729638592341e-15 2.93741988695308e-15 2.40527535614998e-15 1.96953440827399e-15 1.61273251956659e-15 1.32056904857373e-15 1.08133406556444e-15 8.85439017833256e-16 7.25032419923147e-16 5.93685165609687e-16 4.86132848931585e-16 3.98064766478786e-16 3.25951144136961e-16 2.66901663525792e-16 2.18549617862063e-16 1.78957054207487e-16 1.46537100196782e-16 1.19990362096513e-16 9.82528450250347e-17 8.04533079727574e-17 6.58783444093666e-17 5.39437889065878e-17 4.41713037522049e-17 3.61692071453898e-17]';
L = [1 0.999975000000000 0.999941091824183 0.999895452731155 0.999834641482844 0.999754682652989 0.999651344439439 0.999520692117010 0.999359942787708 0.999168502562926 0.998948859208205 0.998706868717057 0.998451112124182 0.998191441812663 0.997937282957856 0.997696324726366 0.997473912628023 0.997273069187451 0.997094884648053 0.996939032531298 0.996804257862736 0.996688772308601 0.996590543817816 0.996507493248839 0.996437618277919 0.996379064517968 0.996330160363800 0.996289428107328 0.996255580341479 0.996227507881133 0.996204263355619 0.996185043147783 0.996169169330508 0.996156072560517 0.996145276434462 0.996136383521136 0.996129063103765 0.996123040560702 0.996118088255840 0.996114017784113 0.996110673410818 0.996107926548287 0.996105671124533 0.996103819712656 0.996102300304918 0.996101053630280 0.996100030928179 0.996099192103963 0.996098504202708 0.996097940147982 0.996097477700627 0.996097098599977 0.996096787856072 0.996096533166736 0.996096324437780 0.996096153388283 0.996096013226024 0.996095898380709 0.996095804284779 0.996095727193372 0.996095664036489 0.996095612297636 0.996095569914221 0.996095535195826 0.996095506757162 0.996095483463065 0.996095464383397 0.996095448756060 0.996095435956668 0.996095425473689 0.996095416888063 0.996095409856505 0.996095404097805 0.996095399381618 0.996095395519266 0.996095392356206 0.996095389765857 0.996095387644545 0.996095385907357 0.996095384484747 0.996095383319760 0.996095382365750 0.996095381584513 0.996095380944764 0.996095380420882 0.996095379991883 0.996095379640584 0.996095379352914 0.996095379117348 0.996095378924450 0.996095378766492 0.996095378637146 0.996095378531229 0.996095378444498 0.996095378373477 0.996095378315321 0.996095378267700 0.996095378228705 0.996095378196774 0.996095378170627 0.996095378149216 0.996095378131684 0.996095378117328 0.996095378105573 0.996095378095947 0.996095378088065 0.996095378081611 0.996095378076325 0.996095378071998 0.996095378068454 0.996095378065553 0.996095378063177 0.996095378061231 0.996095378059638 0.996095378058333 0.996095378057265 0.996095378056390 0.996095378055674 0.996095378055088 0.996095378054607 0.996095378054214 0.996095378053892 0.996095378053629 0.996095378053413 0.996095378053236 0.996095378053091 0.996095378052972 0.996095378052875 0.996095378052796 0.996095378052731 0.996095378052677 0.996095378052634 0.996095378052598 0.996095378052569 0.996095378052545 0.996095378052525 0.996095378052509 0.996095378052496 0.996095378052485 0.996095378052476 0.996095378052469 0.996095378052463 0.996095378052458 0.996095378052454 0.996095378052451 0.996095378052448 0.996095378052446 0.996095378052444 0.996095378052443 0.996095378052442 0.996095378052441 0.996095378052440 0.996095378052439 0.996095378052439 0.996095378052439 0.996095378052438 0.996095378052438 0.996095378052438 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052437 0.996095378052436 0.996095378052436 0.996095378052437 0.996095378052436 0.996095378052436 0.996095378052436 0.996095378052436 0.996095378052436 0.996095378052436 0.996095378052436]';
Ps = 1-3*f; %guess for prices social good
CS = repmat(CS_ss,nT,1);    %guess for aggregate consumption social good

%initialization
L1 = L;
Ps1 = Ps;
f1 = f;
CS1 = CS;

epsL = 1;
epsf = 1;
epsP = 1;
epsCS = 1;
a0 = a;  
while (epsL > 1e-3 || epsf > 1e-4 || epsP > 1e-3 || epsCS > 1e-3) && iter1 < 100      
    iter1 = iter1 + 1;
    
    %update guesses
    L = gamma*L1 + (1-gamma)*L;
    Ps = gammaP*Ps1 + (1-gammaP)*Ps;
    f = gamma*f1 + (1-gamma)*f;
    CS = gammaP*CS1 + (1-gammaP)*CS;
    a = a0*CS.^d_CS;
    
    %income per capita per period
    B = (((As*Ps).^(alpha/(alpha-1))+As*Ps)./((1+(As*Ps).^(1/(alpha-1))).^alpha)).*(L.^(alpha-1));   
     
    V = NaN(nS,nI,nT,nphi);
    cs_fin = NaN(nS,nI,nT,nphi);
    cs_fin(:,:,nT,1)= repmat(phi(1)*B(end)/Ps(end),nS,nI,1);
    cs_fin(:,:,nT,2)= repmat(phi(2)*B(end)/Ps(end),nS,nI,1);
    % assume that after last period everyone that is I moves to R
    % in last period, consume steady state
    V(:,:,nT,1) = ( (1-phi(1))*log((B(end)-Ps(end)*cs_fin(:,:,nT,1))/(1-phi(1))) + phi(1)*log(cs_fin(:,:,nT,1)/phi(1)) + u_lowerbar)/(1-beta); 
    V(:,:,nT,2) = ( (1-phi(2))*log((B(end)-Ps(end)*cs_fin(:,:,nT,2))/(1-phi(2))) + phi(2)*log(cs_fin(:,:,nT,2)/phi(2)) + u_lowerbar)/(1-beta); 

    for i_phi = 1:nphi   
        for t = nT-1:-1:1   %backward induction to solve for consumption      
            cs = XS_mat*phi(i_phi)*B(t)/Ps(t); %consumption social good
            cr = (1-XS_mat*phi(i_phi))*B(t);    %consumption regular good
            degree = b+a(t)*cs;
            S_out = p*f(t)*degree.*S_mat;
            Sprime = max(S_mat-S_out,S_grid(1));    %susceptible next period
            Iprime = I_mat*(1-lambda_R-lambda_D)+S_out; %infected next period

            %flow utility
            U = phi(i_phi)*log(cs) + (1-phi(i_phi))*log(cr) - phi(i_phi)*log(phi(i_phi)) - (1-phi(i_phi))*log(1-phi(i_phi)) + u_lowerbar;
            %value next period depending on probability being in each state
            Wprime = interp2(I_tmp,S_tmp,V(:,:,t+1,i_phi),Iprime,Sprime,'linear',-inf);

            %optimize over consumption
            Wmatrix = U + beta*(1-lambda_D*Itilde_mat).*Wprime;     %
            [Wtemp,cs_ind] = max(Wmatrix,[],3);
            V(:,:,t,i_phi) = Wtemp;
            cs_fin(:,:,t,i_phi) = xs_grid(cs_ind)*phi(i_phi)*B(t)/Ps(t);
        end

        %solve for time path S and I
        for t=1:nT-1    %iterate forward

            cs_eq(t,i_phi) = interp2(I_tmp,S_tmp,cs_fin(:,:,t,i_phi),I_phi(t,i_phi),S_phi(t,i_phi),'linear');

            S_out = p*f(t)*(b+a(t)*cs_eq(t,i_phi)).*S_phi(t,i_phi);
            S_phi(t+1,i_phi) = S_phi(t,i_phi) - S_out; 
            I_phi(t+1,i_phi) = I_phi(t,i_phi)*(1-lambda_R-lambda_D)+S_out;

        end
        cs_eq(nT,i_phi) = interp2(I_tmp,S_tmp,cs_fin(:,:,nT,i_phi),I_phi(nT,i_phi),S_phi(nT,i_phi),'linear');


    end
    degree_eq = b+repmat(a,1,nphi).*cs_eq;


    frac_alive = 1-(1-S_phi-I_phi)*lambda_D/(lambda_R+lambda_D);

    L1 = frac_alive*g';         %mass of agents alive (and labor supply)
    f1 = ((degree_eq.*I_phi)*g')./((degree_eq.*frac_alive)*g'); %probability that a person you meet is infected 

    cs_demand = (cs_eq.*frac_alive)*g';
    Ps1 = ((max((As^(1/alpha)) * L1./(cs_demand.^(1/alpha)) - 1,0.00001)).^(alpha-1))/As;    %Price at which social market clears. make sure that term to the power is positive

    CS1 = cs_demand;

    eps1 = norm(L1-L) + norm(f1-f) + norm(Ps1-Ps);
    epsL = norm(L1-L);
    epsf = norm(f1-f);
    epsP = norm(Ps1-Ps);
    epsCS = norm(CS1-CS);
    fprintf('iter big loop=%d epsL=%1.8f epsf=%1.8f epsP=%1.8f epsCS=%1.8f\n',iter1,epsL,epsf,epsP,epsCS) 
end


B_fin = B;

I_overall = I_phi*g';

value1 = interp2(I_tmp,S_tmp,V(:,:,1,1),I_phi(1,1),S_phi(1,1),'linear');
value2 = interp2(I_tmp,S_tmp,V(:,:,1,2),I_phi(1,2),S_phi(1,2),'linear');
        
value_competitive = g(1) * value1 +  g(2) * value2;

%% create plots

nT2 = nT-100;
%{
%%% figure 5 in paper 
figure(1)
yyaxis left
plot(T_grid(1:nT2),I_overall(1:nT2),'linewidth',2,'Color','b')
hold on
plot(T_grid(1:nT2),I_phi((1:nT2),2),'linewidth',2,'Color','r')
plot(T_grid(1:nT2),I_phi((1:nT2),1),'linewidth',2,'Color','g')
plot(T_grid(1:nT2),f(1:nT2),'linewidth',2,'Color','m')
ylabel('Fraction of Initial Population','FontSize',24)
yyaxis right
plot(T_grid(1:nT2),1-frac_alive((1:nT2),:)*g','linewidth',2,'Color','k')
ylabel('Deaths as Fraction of Initial Population','FontSize',24)
legend('Infected','Infected high degree','Infected low degree','Prob meeting is infected','Death (right axis)')
% legend('Infected','Death (right axis)')
xlabel('Time (weeks)','FontSize',24)
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ax.YAxis(1).Limits = ([0, 0.3]);
ax.YAxis(2).Limits = ([0, 0.005]);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 13)
hold off

%%% income per capita over time
figure(2)
plot(T_grid(1:nT2),B_fin(1:nT2)/B_ss,'linewidth',2,'Color','k')
hold on
ylabel('Income per capita relative to steady state','FontSize',24)
xlabel('Time (weeks)','FontSize',24)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 13)
hold off

%%% Consumption social good
figure(3)
plot(T_grid(1:nT2),cs_eq((1:nT2),1)/cs_ss(1),'linewidth',2)
hold on
plot(T_grid(1:nT2),cs_eq((1:nT2),2)/cs_ss(2),'--','linewidth',2)
legend('Low degree','High degree')
ylabel('Cons social good relative to steady state','FontSize',24)
xlabel('Time (weeks)','FontSize',24)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 13)
ax = gca;
ax.YAxis.Limits = ([0.4, 1.1]);
hold off

%%% Degree
figure(4)
plot(T_grid(1:nT2),degree_eq((1:nT2),1)/degree_ss(1),'linewidth',2)
hold on
plot(T_grid(1:nT2),degree_eq((1:nT2),2)/degree_ss(2),'--','linewidth',2)
legend('Low degree','High degree')
ylabel('Degree relative to steady state','FontSize',24)
xlabel('Time (weeks)','FontSize',24)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 13)
ax = gca;
ax.YAxis.Limits = ([0.4, 1.1]);
hold off
%}






