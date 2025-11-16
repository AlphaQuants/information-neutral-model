global dt;
dt = 0.02;

global dk
dk=0.005;
global S0;
S0 =100;

global CallStrikeInsiderTrader;
CallStrikeInsiderTrader = 100;

global PutStrikeInsiderTrader;
PutStrikeInsiderTrader = 0;

global KOBarrier;
KOBarrier = 1E+8;


global CallStrikeLargeTrader;
CallStrikeLargeTrader = 110;
global PutStrikeLargeTrader;
PutStrikeLargeTrader = 0;

global sigma;
sigma = 0.25;

global Xmin;
Smin = 10;
global Xmax;

Xmid= log(100);
Xmin = log(Smin);
ndiscX = 300;

global dX;
dX = (Xmid-Xmin)/ndiscX*3.0/2.0;
dS = 1;
Xmax= Xmin +dX*ndiscX;
Smax =exp(Xmax);
global Xgrid;
Xgrid=Xmin:dX:Xmax;

global kgrid kmin;
kmin=0.0;
kgrid=kmin:dk:1.0;
global T;
T= 0.33;
t0 =0;

% how much time to settle?
% lets assume 1 month

 
 global dtau;
 dtau = 1.0/100.0;

% the daily volume average is of percentage of notional
global DVA ADV;
DVA= 5;
ADV= DVA;
global tmpCostAlpha;
tmpCostAlpha = 2;
% lambda how much move in percentage

global lambda;
lambda = 0.2;%0.1/DVA;



global mu;
mu = 0.02;

global r;
r= 0.01;

global phi;
phi = 20.0;

global nu;
nu =  0.1;
global lambda_term;

global tau;
tau = 1.0/24.0;
lambda_term = nu;

global rho;
rho =0.5;

global kappaStar;
global varrho;
var = rho + phi ;

lambda_values = 0 : 0.1 : 2.0;       % 21 values
nL = numel(lambda_values);

res_store      = zeros(1, nL);      % stores res(201,57)
insider_store  = zeros(1, nL);      % stores InsiderOption(201,57)
bs_store       = zeros(1, nL);      % if you want BSOption as well

for idx = 1:nL
    lambda = lambda_values(idx);    % update the global lambda
    assignin('base','lambda',lambda);   % set global lambda

    [res, InsiderOption, BSOption] = PDESolver3();

    % store the values you want
    res_store(idx)     = res(201,57)
    insider_store(idx) = InsiderOption(201,57)
    bs_store(idx)      = BSOption(201,57);  % optional
end

% After the loop:
disp(res_store);
disp(insider_store);


%[res,InsiderOption,BSOption]=PDESolver3();


[fvC, resgridC] = calcGreeks2(exp(Xgrid),CallStrikeLargeTrader, r*100, 0, sigma*100, T, 'c', 1);
[fvP, resgridP] = calcGreeks2(exp(Xgrid),PutStrikeLargeTrader, r*100, 0, sigma*100, T, 'put', 1);
plot( exp(Xgrid), resgridC.delta)
hold on
plot(exp(Xgrid),kappaStar(:,1)*dk);
%plot(exp(Xgrid),kappaStar(:,round(size(kgrid,2)/2))*dk+ kgrid(round(size(kgrid,2)/2)));
hold off
legend('Black Scholes Delta','Large Trader Delta rate \phi\chi_0  when V_0=0', 'Location','northwest')
ylabel('Delta to be executed')
xlabel('Spot (S_0)')
set(gca,'XLim',[40 160])
set(gca,'XTick',(40:20:160))

%exp(Xgrid),InsiderOption(:,kgridp)+res(:,kgridp),...
kgridp=134;
plot(exp(Xgrid),BSOption(:,kgridp) ,exp(Xgrid),InsiderOption(:,kgridp),...
    exp(Xgrid),max(exp(Xgrid)-CallStrikeInsiderTrader,0), '--' ,...
'LineWidth',1)

legend({'Black Scholes Price','Insider Price', 'Insider Price + Large Trader Cost Function' , 'Call Payoff'},...
     'Location','northwest')
 set(gca,'XLim',[40 160])
set(gca,'XTick',(40:20:160))

xlabel('Spot (S_0)')


meshc(kgrid,exp(Xgrid), TerminalCost());
xlabel('Current Shares (V_T)')
ylabel('Spot (S_T)')
zlabel('Cost Function J(T,.,.)')

plot(exp(Xgrid(1:200)),BSOption(1:200,180));
hold on
plot(exp(Xgrid(1:200)),InsiderOption(1:200,180));
xlabel('Spot (S_0)')
ylabel('Option price (O(0,S_0))')
legend('Black Scholes price','Insider trader price', 'Location','northwest')

BSOption(136,1)