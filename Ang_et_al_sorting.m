%% This document performs the same estimating process as Ang-et-al-2006 %%
% To reproduce table 1 with concentration and sparsity  we can do the following:
% Either regress the excess returns on market factor and NC or NS
% (separately or together)
% or just regress on NC and NS 

clear;
clc;

%% For now just regress on NS and NC and Mrk and sort %%

%monthly return for common stock from 1976-2018
load('raw_CS_return.mat'); 

% monthly cap for common stock from 1976-2018
load('raw_CS_cap.mat'); 

%market factor
load('mrkt_96_2013.mat');  

%risk free rate
load('RF_96_2013.mat'); 

%loading factors (concentration and sparsity)
load('Concentration2.mat'); 
load('Concentration1.mat');

%load VIX
load('VIXmonthly.mat') 
%inno_VIX = diff(VIX);

CONCENTRATION = VIXmonthly;

%% for implied volatility, lets only look at 1996 to 2013 total of 216 months %%
t1 = datetime(1975,12,1);
t = t1 + calmonths(1:516);
t = t';

% PERMNO header for each col
permno = raw_CS_return(1,:); 

%common stock return to Jan 2014-one month ahead
cs_rtn = raw_CS_return(242:458,:); 

%common stock market cap to Jan 2014-one month ahead
cs_cap = raw_CS_cap(242:458,:); 

%% Test for NC in innovation %%
inno_NC = diff(CONCENTRATION);

%common stock return
cs_rtn = raw_CS_return(243:457,:); 

%common stock market cap
cs_cap = raw_CS_cap(243:457,:); 
mrkt_96_2013 = mrkt_96_2013(2:end,:);
RF_96_2013 = RF_96_2013(2:end,:);
CONCENTRATION = inno_NC;
%%
% 
% R = cs_rtn;
% P = permno;
% m = size(R,1)/24; %number of periods biennial
% for i=1:m
%     R = cs_rtn;
%     P = pemno;
%     bi_cs_rtn = R(1:24,:);
%     num_NaN = isnan(bi_cs_rtn); %find NaN in each col
%     total_num_NaN = sum(num_NaN); %get the total number of NaN for two years
%     NaN_rmv = find(total_num_NaN >= 3);%index of col that have more than 3 NaN and should be removed
%     bi_cs_rtn(:,NaN_rmv) = [];
%     fillmissing(bi_cs_rtn,'linear'); %fill the missing by linear interpolation
%     P(:,NaN_rmv) = [];
%     P_R = [P;bi_cs_rtn]; %return and their permno as header
%     P_R_CELL{i} = P_R; %the cell that contain biennial data
%     R(1:24,:) = [];% remove the first two years in each loop
% end


%% Rolling window analysis: at each period t we construct portfolios by %%
% regressing values form t-1 to t-w - same as Herskovic

% window length
w = 24; 
k = size(cs_rtn,1);

% number of intervals/number of portfolios
N = k - w ; 

for i=1:N
P = permno; %global set of permno
R = cs_rtn; %global set of return from 96 to 2013
C = cs_cap; %global market cap
P_t = permno;
R_t = R(w+1,:); %retun on period t
C_t = C(w+1,:);
num_NaN = isnan(R_t); %find NaN in each col of R_t
NaN_rmv = find(num_NaN == 1); %firms that dont have return at t
R_t(:,NaN_rmv) = [];%returns of firms that return exist at t
P_t(:,NaN_rmv) = []; %permnos of firms that return exist at t
C_t(:,NaN_rmv) = [];

num_NaN = isnan(C_t); %do the same circle for Market caps
NaN_rmv = find(num_NaN == 1); 
C_t(:,NaN_rmv) = [];
P_t(:,NaN_rmv) = [];
R_t(:,NaN_rmv) = [];

% We should look at the same firms in P_t and find their returns fron t-1 to
% t-w by searching P_t in R and consider those that have at least 21 points 
size_P_t = size(P_t,2); %how many firms to look for
P_t_1 = P_t;
r_t_1 = zeros(w,size_P_t); %return series from t-w to t-1
c_t_1 = zeros(w,size_P_t);
for j=1:size_P_t
    pr = P_t(1,j); %pick the j permno in P_t
    idx = find(P == pr); %find the index in global permno that permno is pr
    r = cs_rtn(1:w,idx); %return of that pr
    c = cs_cap(1:w,idx); %the same for cap
    r_t_1(:,j) = r; 
    c_t_1(:,j) = c;
end

% check for NaN in r_t_1
num_NaN_t_1 = isnan(r_t_1); 
total_num_NaN_t_1 = sum(num_NaN_t_1);
NaN_rmv = find(total_num_NaN_t_1 >= 3);

% remove them
r_t_1(:,NaN_rmv) = []; 
P_t_1(:,NaN_rmv) = [];
r_t_1_noNaN = fillmissing(r_t_1,'linear');
R_t_1{i} = [P_t_1;r_t_1_noNaN];


% we need to go back to R_t and P_t and remove R(t) and their permno that
% were not available in the window frame
r_t = zeros(1,size(r_t_1,2));
cap = zeros(1,size(r_t_1,2));
for k=1:size(r_t_1,2)
    pr2 = P_t_1(:,k);
    idx2 = find(P_t == pr2); %find the index in permno that permno is pr
%     cap_permno = C_t(1,:);
%     idx3 = find(cap_permno == pr2);
    cap(1,k) = C_t(1,idx2);
    r_t(1,k) = R_t(:,idx2); 
end

R_t_cell{i} = [P_t_1;r_t;cap];


cs_rtn(1,:) = [];
cs_cap(1,:) = [];
end
cs_rtn = raw_CS_return(242:457,:);
cs_cap = raw_CS_cap(242:457,:); %common stock market cap

%% creating window frames for RF and Market factor and Concentration %%
rf = RF_96_2013;
mrkt = mrkt_96_2013;

MARKET = zeros(w,N);
RF = zeros(w,N);
NC = zeros(w,N);

for i=1:N
MARKET(:,i) = mrkt(1:w,1);
RF(:,i) = rf(1:w,1);
NC(:,i) = CONCENTRATION(1:w,1);
mrkt(1,:) = [];
rf(1,:) = [];
CONCENTRATION(1,:) = [];
end

%% Running regressions for windows %%
% we need R_t_1 & R_t_cell , MARKET and RF and a column of one each time
% and also the factors. For now: concentration with 2 versions
% dont forget to adjust the dimension since you would regress on both level
% and innovation

%beta1 = zeros(1,N);
%beta2 = zeros(1,N);
%beta1_CI = cell(1,N);
%beta2_CI = cell(1,N); No CI for now, I can add them later

% declaring the independent variables
for i=1:N
    depends = R_t_1{i}(2:end,:); %set of y
    permno = R_t_1{i}(1,:);
    beta1 = zeros(1,size(depends,2));
    beta2 = zeros(1,size(depends,2));
    x1 = MARKET(:,i);
    x2 = NC(:,i);
    rfree = RF(:,i);
    for j=1:size(depends,2)
    y = depends(:,j) - rfree; %excess return for permno j in R_t_1 i;
    intercept = ones(w,1);
    X = [intercept, x1, x2];
    [b,bint] = regress(y,X);
    beta1(:,j) = b(2,1);
    beta2(:,j) =  b(3,1);
    end
    BETA_NC{i} = [permno;beta2]; %storing the preformation betas of the factor concentartion
end

% Sorting BETA_NC
for i=1:N
    bnc = BETA_NC{i}(2,:);
    permno = BETA_NC{i}(1,:);
    IndexNC = zeros(1,size(bnc,2));
    sorted_betanc = zeros(1,size(bnc,2));
    [B,I] = sort(bnc);
    IndexNC(1,:) = I';
    sorted_betanc(1,:) = B';
    permn = zeros(1,size(IndexNC,2)); %to selec the permno as well
    for j=1:size(IndexNC,2)
        ind = IndexNC(:,j);
        prm = permno(1,ind);
        permn(1,j) = prm;
    end
    sortedNC{i} = [permn;IndexNC;sorted_betanc];  %%yYOU HAVE WRITE THE CORRECT PERMNO
end

% In sortedNC, first row corresponds to permno, 2nd row is index of the firm,
% 3rd is the sorted return
% creating quintile portfolios for each period N=192
qnt = 5;
for i=1:N


% total sorted firms for the window
totalsorted = sortedNC{i}; 

% dimension of each quitile (num of firms in each)
q1_num = round(size(totalsorted,2)/qnt);
q5_num = q1_num; 
q2_num = q1_num;
q4_num = q1_num;
q3_num = size(totalsorted,2)-q1_num-q2_num-q4_num-q5_num; %the rest in the middle
quintile_sizes = [q1_num,q2_num,q3_num,q4_num,q5_num];
QUINTILES{i} = mat2cell(totalsorted,3,quintile_sizes); %in each cell there is one quintile portfolio

end

%% Creating portfolios in period (t) from QUINTILES %%
% we look at permnos in each QUINTILES and return their corresponding return
% and market cap in period t in R_t_cell. so we can creat the value-weighted
% returns.

PRTFL = cell(N,qnt); %total portfolios
for i=1:N
    all_quantiles = QUINTILES{1,i};
    Rt = R_t_cell{i};
    for j=1:qnt
        quantile = QUINTILES{1,i}{1,j};
        Rt_return_cap = zeros(3,size(quantile,2));
        for k=1:size(quantile,2)
            pr = quantile(1,k);
            retnIndex = find(Rt(1,:) == pr);
            retrn = Rt(2,retnIndex);
            caps = Rt(3,retnIndex);
            Rt_return_cap(2,k) = retrn; %store its return in R_t
            Rt_return_cap(3,k) = caps; %store its market cap
            Rt_return_cap(1,k) = pr; %store the prmno
            
        end
        PRTFL{i,j} = Rt_return_cap;
    end
    
    
end

% calculate value weighted returns for each portfolio
prtfl_returns = zeros(N,qnt);
for i=1:N
    for j=1:qnt
        
        prtfl = PRTFL{i,j};
        shares = prtfl(3,:)./sum(prtfl(3,:));
        prtfl_returns(i,j) = sum(prtfl(2,:).*shares);
        
    end
end

%% test %%
[h,p,ci,stats]=ttest2(prtfl_returns(:,1),prtfl_returns(:,5))
