% sum rate 
%% OPTIMAL VS SUBOPTIMAL NOMA ASR.


% sum rate 


clear; close all;
format compact

 
K_V = [1 2 4:4:32]; 
Rate_V_OMA1= [0.7077   0.8221    0.8490    0.8791    0.9361    0.9689    1.166    1.3189    1.5380    1.7550];
Rate_V_SNOMA1 = [0.8352    1.0019    1.2124    1.7019    1.8946    1.9757   1.9991    2.0187    2.1323    2.1787];

Rate_V_NOMA2 = [0.8352    1.1068    1.4207    1.7414    1.9473    2.0778    2.1812    2.2598    2.3284    2.3871];

Rate_V_SNOMA2= [0.8352    1.1759    1.4732    1.7642    1.9977    2.1823    2.3008    2.4175    2.4746   2.5040];



figure(1)
         x = K_V;
        Y = [Rate_V_OMA1; Rate_V_SNOMA1; Rate_V_NOMA2;Rate_V_SNOMA2];
        plot(x,Y,'MarkerSize',12,'LineWidth',2);
        legend('OMA','Suboptimal','Optimal','Proposed','')
        xlabel('number of users');
        ylabel('average sum rate in bits/channel use');





%% OPTIMAL/SUBOPTIMAL VS VARYING USERS
% sum rate for 
% LS + Rayleigh fading


% sum rate for 
% LS + Rayleigh fading


clear; 

format compact

K_V = [1 2 4:4:32 40:8:64];

%Nr = 64;  Nt = 1;  P_dB = 0 dB;
Rate_V_Capacity_SUB = [2.0616    3.5490    5.6884    9.0317   11.8165   14.2726   16.4451   18.3394   20.1641   21.9163   25.1423   28.0123   30.7569   33.2954];
Rate_V_Capacity_OMA = [2.0616    3.1997    4.7201    6.8501    8.4214    9.6442   10.6390   11.5598   12.4160   13.1636   14.4782   15.5561   16.5490   17.3647];

Rate_V_Capacity_PR = [4.3793    7.4572   12.2419   19.8583   26.2659   31.9470   36.9913   41.6205   45.9480   50.0451   57.6386   64.4988   70.9077   76.9589];
Rate_V_Capacity_OP = [4.3793    7.2702   11.6411   18.2427   23.4094   27.6697   31.2928   34.6301   37.7190   40.4581   45.4734   49.7621   53.6306   56.9793];


figure(2)


        x = K_V;
        Y = [Rate_V_Capacity_SUB;Rate_V_Capacity_OMA; Rate_V_Capacity_PR;Rate_V_Capacity_OP];
        plot(x,Y,'MarkerSize',12,'LineWidth',2);
        legend('suboptimal','Oma','Proposed','Optimal','')
        xlabel('number of users');
        ylabel('average sum rate in bits/channel use');






%% end

%% ASR CHANNEL INTRA/INTERFERANCE.


clear;  
format compact;

%K_V = [1 2 4:4:32 40:8:64];
K_V = [1 2 4:4:32];
Rate_V_NOMA = zeros(1,length(K_V));
Rate_V_OMA = zeros(1,length(K_V));

N0 = 1;
p_sum = 10.^(10/10);

display_num = 1e1;
Channel_Num = 1e3;


PL_Factor = -3.76;
gTh = (35/289)^PL_Factor;
PL_Avrg = 55.7484;  %d0 = 35
m_log = 5.45579259206826;

ls_avrg = 0;
rng('default');
for knum_ix = 1:length(K_V)
    User_Num = K_V(knum_ix);
    snr_ind = p_sum/User_Num/N0;
    snr_sum = p_sum/N0;

    for h_ix = 1:Channel_Num

        if mod(h_ix,display_num) == 1
            User_Num,h_ix,tic
        end
        
        G_All = zeros(1,User_Num);
        for k_ix = 1:User_Num
            
            % the fading with PL > gTh is discarded
            Ploss1 = Inf;
            while Ploss1 > gTh
                a = rand(1,1); % uniform distribution with unit interval.
                a = sqrt(a) * sqrt(3)/2;

                sign = (-1)^randi([0,1]);
                b = rand(1,1);
                b = sign*a*b/sqrt(3);
                clear sign;

                Ploss1 = norm([a,b])^PL_Factor;
                clear a b                 
            end
            lognormal = 10.^(randn(1,1)*0.8);
            g_LS = Ploss1*lognormal/PL_Avrg/m_log;
            
            Rayleigh_Gain = abs( randn(1,1) + 1i*randn(1,1) )^2 /2; 
            G_All(k_ix) = g_LS *Rayleigh_Gain;
          end        
        ls_avrg = ls_avrg + sum(G_All);
        
        rate_noma = real( log2(1 + sum(G_All)*snr_ind) );
        rate_oma = sum( real( log2(1 + G_All*snr_sum) ) )/User_Num;

        Rate_V_NOMA(knum_ix) =  Rate_V_NOMA(knum_ix) + rate_noma;
        Rate_V_OMA(knum_ix) =  Rate_V_OMA(knum_ix) + rate_oma;

        if mod(h_ix,display_num) == 0
            toc
        end
    end
end

Rate_V_NOMA = Rate_V_NOMA/Channel_Num;
Rate_V_OMA = Rate_V_OMA/Channel_Num;
ls_avrg = ls_avrg/Channel_Num/sum(K_V);

figure(3);
plot(K_V,Rate_V_NOMA,'b-*','MarkerSize',12,'LineWidth',2);hold ; grid 
plot(K_V,Rate_V_OMA,'r-o','MarkerSize',12,'LineWidth',2);hold ; grid ;
xlabel('users');
ylabel('average sum rate in bits');
legend('NOMA','OMA');


%% GROUPING ALGORITHM.
%users vs resource blocks
%Grouping algorithm

clc; clear variables; 

N = 10^5;
Pt = 30;                        %Max BS Tx power (dBm)
pt = (10^-3)*db2pow(Pt);        %Max BS Tx power (Linear scale)
No = -114;                      %Noise power (dBm)
no = (10^-3)*db2pow(No);        %Noise power (Linear scale)


r = 0.5:0.5:10;                 %Far user target rate range (R*)

df = 1000; dn = 500;            %Distances

eta = 4;                        %Path loss exponent

p1 = zeros(1,length(r));
p2 = zeros(1,length(r));
pa1 = zeros(1,length(r));
pa2 = zeros(1,length(r));

af = 0.75; an = 0.25;       %Fixed PA (for comparison)



hf = sqrt(df^-eta)*(randn(1,N) + 1i*randn(1,N))/sqrt(2);
hn = sqrt(dn^-eta)*(randn(1,N) + 1i*randn(1,N))/sqrt(2);

g1 = (abs(hf)).^2;
g2 = (abs(hn)).^2;

for u = 1:length(r)
    epsilon = (2^(r(u)))-1;         %Target SINR for far user

        %BASIC FAIR PA%
%        aaf = min(1,epsilon*(no + pt*g1)./(pt*g1*(1+epsilon)));
%        aan = 1 - aaf;

     %IMPROVED FAIR PA%
     aaf = epsilon*(no + pt*g1)./(pt*g1*(1+epsilon));
     aaf(aaf>1) = 0;
     aan = 1 - aaf;



    gamma_f = pt*af*g1./(pt*g1*an + no);
    gamma_nf = pt*af*g2./(pt*g2*an + no);
    gamma_n = pt*g2*an/no;

    gamm_f = pt*aaf.*g1./(pt*g1.*aan + no);
    gamm_nf = pt*aaf.*g2./(pt*g2.*aan + no);
    gamm_n = pt*g2.*aan/no;


    Cf = log2(1 + gamma_f);
    Cnf = log2(1 + gamma_nf);
    Cn = log2(1 + gamma_n);

    Ca_f = log2(1 + gamm_f);
    Ca_nf = log2(1 + gamm_nf);
    Ca_n = log2(1 + gamm_n);

    for k = 1:N
        if Cf(k) < r(u)
            p1(u) = p1(u) + 1;
        end
        if (Cnf(k)<r(u))||(Cn(k) < r(u))
            p2(u) = p2(u) + 1;
        end

        if Ca_f(k) < r(u)
            pa1(u) = pa1(u) + 1;
        end
        if aaf(k) ~= 0
            if (Ca_n(k) < r(u)) || (Ca_nf(k) < r(u))
                pa2(u) = pa2(u) + 1;
            end
        else
            if Ca_n(k) < r(u)
                 pa2(u) = pa2(u) + 1;
            end
        end
    end
end
pout1 = p1/N;
pout2 = p2/N;
pouta1 = pa1/N;
pouta2 = pa2/N;



figure(4);
plot(r,pout1,'--+r','linewidth',2); hold ; grid ;
plot(r,pout2,'--ob','linewidth',2);
plot(r,pouta1,'r','linewidth',2);
plot(r,pouta2,'b','linewidth',2);
plot(r,pouta2,'b','linewidth',2)
xlabel('Target Resource block');
ylabel('Outage probability of users');
xlim([r(1) r(end)]);
legend('K+1 user (fixed PA)','K user (fixed PA)','K+1 user (fair PA)','K user (fair PA)');



%% end

%% power allocation variables algorithm decording order.

clc; clear variables; 

N = 10^5;
df = 5000; dn = 1000; %Distances

eta = 4;
hf = sqrt(df^-eta)*(randn(1,N) + 1i*randn(1,N))/sqrt(2);
hn = sqrt(dn^-eta)*(randn(1,N) + 1i*randn(1,N))/sqrt(2);

gf = (abs(hf)).^2;
gn = (abs(hn)).^2;

R1 = 1;  %Target rate bps/Hz

epsilon = (2^(R1))-1; %Target SINR

%Transmit power
Pt1 = 0:30;
pt11 = (10^-3)*db2pow(Pt1);

%Noise power
No = -114;
no = (10^-3)*db2pow(No);

b1 = 0.75; b2 = 0.25; %Fixed PA for comparison
for u = 1:length(pt11)
    a1 = epsilon*(no + pt11(u)*gf)./(pt11(u)*gf*(1+epsilon));
    a1(a1>1) = 0;
    a2 = 1 - a1;
    
    %Sum rate of fair PA
    C1 = log2(1 + pt11(u)*a1.*gf./(pt11(u)*a2.*gf + no));
    C2 = log2(1 + pt11(u)*a2.*gn/no);
    
    C_sum(u) = mean(C1+C2);
    
    %Sum rate of fixed PA
    C1f = log2(1 + pt11(u)*b1.*gf./(pt11(u)*b2.*gf + no));
    C2f = log2(1 + pt11(u)*b2.*gn/no);
    C_sumf(u) = mean(C1f+C2f);
end

figure(5)
plot(Pt1,C_sum,'linewidth',1.5); 
hold ; 
grid ;
plot(Pt1,C_sumf,'linewidth',1.5);
legend('P.optimal','Optimal')
xlabel('Transmit power (dBm)');
ylabel('Sum rate (bps/Hz)');

%% END

%% DECORDING ORDER

clc; 
clear variables;


load decordingworkspace1.mat
%Create QPSKModulator and QPSKDemodulator objects
QPSKmod = comm.QPSKModulator('BitInput',true); 
QPSKdemod = comm.QPSKDemodulator('BitOutput',true); 

%Perform QPSK modulation
xmod1 = step(QPSKmod, x1);
xmod2 = step(QPSKmod, x2);
xmod3 = step(QPSKmod, x3);

%Do super position coding
x = sqrt(a1)*xmod1 + sqrt(a2)*xmod2 + sqrt(a3)*xmod3;

for u = 1:length(Pt2)
	
%Received signals
   y1 = sqrt(pt22(u))*x.*h1 + n1;	
   y2 = sqrt(pt22(u))*x.*h2 + n2;	
   y3 = sqrt(pt22(u))*x.*h3 + n3;
   
%Perform equalization
   eq1 = y1./h1;
   eq2 = y2./h2;
   eq3 = y3./h3;
   
%Decode at user K-2 (Direct decoding)
   dec1 = step(QPSKdemod, eq1);
   
%Decode at user K-1
   dec12 = step(QPSKdemod, eq2);		%Direct demodulation to get U1's data
   dec12_remod = step(QPSKmod, dec12);		%Remodulation of U1's data
   rem2 = eq2 - sqrt(a1*pt22(u))*dec12_remod;	%SIC to remove U1's data
   dec2 = step(QPSKdemod, rem2);		%Direct demodulation of remaining signal
   
%Decode at user K
   dec13 = step(QPSKdemod, eq3);		%Direct demodulation to get U1's data
   dec13_remod = step(QPSKmod, dec13);		%Remodulation of U1's data
   rem31 = eq3 - sqrt(a1*pt22(u))*dec12_remod;	%SIC to remove U1's data
   dec23 = step(QPSKdemod, rem31);		%Direct demodulation of remaining signal to get U2's data
   dec23_remod = step(QPSKmod, dec23);		%Remodulation of U2's data
   rem3 = rem31 - sqrt(a2*pt22(u))*dec23_remod;	%SIC to remove U2's data
   dec3 = step(QPSKdemod, rem3);		%Demodulate remaining signal to get U3's data
   
%BER calculation
   ber1(u) = biterr(dec1, x1)/N;
   ber2(u) = biterr(dec2, x2)/N;
   ber3(u) = biterr(dec3, x3)/N;
end

figure(6)
semilogy(Pt2, ber1, '-o', 'linewidth', 2);
hold; 
grid;
semilogy(Pt2, ber2, '-o', 'linewidth', 2);
semilogy(Pt2, ber3, '-o', 'linewidth', 2);

xlabel('Transmit power (dBm)');
ylabel('BER');
legend('P.OPTIMAL', 'OPTIMAL', 'OMA');



%% END
