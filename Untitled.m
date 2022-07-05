clc;
% close all;
clear all;

f = 9.5*10^6;
c = 3*10^8;
lambda = c/f;

N =6;
d =[0.45*lambda  0.6*lambda]
amp = [0 1];
ang = [0 360];

fitness_condition = 1;

Main_beam_dir = 0;
nulls_locations=[40];

% nulls_locations = [0 20 22 24 23 26 28 29 30 32 34 36 38 39 40 42 44 46 48  50 52 54 20+90 22+90 24+90 23+90 26+90 28+90 29+90 30+90 32+90 34+90 36+90 38+90 39+90 40+90]; ;
% nulls_locations=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 28 29 30 31 32 32 34 35 36 37 38 39 40 41 42 43 44 45 46  47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74];
%  nulls_locations=0:1:80;
% for k=1:179;
%     nulls_locations (1)=1;
%    nulls_locations (k+1)=nulls_locations (k)+nulls_locations (1);
% end

beamwidth = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dist_val = 2;
S = 500;
Itt_max = N+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = 1;
phi_p = 0.9;
phi_g = 0.9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% CODE ADDDED %%%%%%%%%%%%%%%%%

%%% I am computing individual distance and storing them in vectors now.
%%% dist_val defines that how many realizations you need.
% d_elem_1 = (d(2)-d(1)).*rand(1,dist_val) + d(1);
% d_elem_2 = (d(2)-d(1)).*rand(1,dist_val) + d(1);
% d_elem_3 = (d(2)-d(1)).*rand(1,dist_val) + d(1);
% d_elem_4 = (d(2)-d(1)).*rand(1,dist_val) + d(1);
% d_elem_5 = (d(2)-d(1)).*rand(1,dist_val) + d(1);
% d_elem_6 = (d(2)-d(1)).*rand(1,dist_val) + d(1);
% d_elem_7 = (d(2)-d(1)).*rand(1,dist_val) + d(1);
%
%
% % keyboard;
% ii = 1;                       %% I have initialized ii to be 1 and will increment it in the below for loops.
% tic
% for d_1 = 1:length(d_elem_1)
%     for d_2 = 1:length(d_elem_2)
%         for d_3 = 1:length(d_elem_3)
%             for d_4 = 1:length(d_elem_4)
%                 for d_5 = 1:length(d_elem_5)
%                     for d_6 = 1:length(d_elem_6)
%                         for d_7 = 1:length(d_elem_7)
%
% %% Now i am constructing the distance_bw_elem vector %%
% %% NO CHANGE IN REST PART OF CODE %%
%
% distance_bw_elem = [d_elem_1(d_1) d_elem_2(d_2) d_elem_3(d_3) d_elem_4(d_4) d_elem_5(d_5) d_elem_6(d_6) d_elem_7(d_7)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:S
    
    distance_bw_elem1 =  (d(2)-d(1)).*rand(1,N-1) + d(1);
    distance_bw_elem2 =  (d(2)-d(1)).*rand(1,N-1) + d(1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    distance_bw_elem= distance_bw_elem2-distance_bw_elem1;
    distance(1) = 0;
    
    
    
    for indx_d = 1:length(distance_bw_elem)
        
        distance(indx_d+1) = distance_bw_elem(indx_d) + distance(indx_d);
        
    end
    
    amplitude = (amp(2)-amp(1)).*rand(1,N) + amp(1);
    phase    = (ang(2)-ang(1)).*rand(1,N) + ang(1);
    
    
    x_d(ii,:) = distance_bw_elem;
    x_amp(ii,:) = amplitude;
    x_ang(ii,:) = phase;
    
    p_d(ii,:)  =  x_d(ii,:);
    p_amp(ii,:) = x_amp(ii,:);
    p_ang(ii,:) = x_ang(ii,:);
    
    [AF fitness_p(ii)]=compute_fitness(N,distance,amplitude,phase,lambda,fitness_condition,Main_beam_dir,nulls_locations,beamwidth);
    %[AF fitness_p(ii)]=compute_fitness(N,distance,amplitude,phase,lambda,fitness_condition,Main_beam_dir,nulls_locations,beamwidth);
    
    if (ii==1)
        
        g_d = p_d(1,:);
        g_amp = p_amp(1,:);
        g_ang = p_ang(1,:);
        fitness_g = fitness_p(1);
        
    elseif (fitness_p(ii)  < fitness_g)
        
        g_d = p_d(ii,:);
        g_amp = p_amp(ii,:);
        g_ang = p_ang(ii,:);
        fitness_g = fitness_p(ii);
        
    end
    
    v_d(ii,:) = -2*(abs(max(distance_bw_elem) - min(distance_bw_elem)))*rand(1,length(distance_bw_elem)) +  (abs(max(distance_bw_elem) - min(distance_bw_elem))) ;
    v_amp(ii,:) = -2*(abs(max(amplitude) - min(amplitude)))*rand(1,length(amplitude)) + (abs(max(amplitude) - min(amplitude)));
    v_ang(ii,:) = -2*(abs(max(phase) - min(phase)))*rand(1,length(phase)) + (abs(max(phase) - min(phase)));
    
    clear distance;
    clear AF;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    ii = ii+1;
    %
    %                         end
    %                     end
    %                 end
    %             end
    %         end
    %     end
    % end
    % %toc;
    % %keyboard
    %
    % total_time = toc*Itt_max
    % fprintf('For the given number of Itterations it will take %d minutes \n',total_time/60)
    % choice = input('Enter 1 to continue or 0 to stop the program here');
    %
    % if (choice == 1)
    %
    %     S = length(x_d)
    %S = 100
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear ii

for outer_itt = 1:Itt_max
    
    for ii= 1:S
        
        r_p = rand(1);
        r_g = rand(1);
        
        v_d(ii,:) = w*v_d(ii,:) + phi_p*r_p*(p_d(ii,:) -x_d(ii,:)) + phi_g*r_g*(g_d-x_d(ii,:));
        v_amp(ii,:) = w*v_amp(ii,:) + phi_p*r_p*(p_amp(ii,:) -x_amp(ii,:)) + phi_g*r_g*(g_amp-x_amp(ii,:));
        v_ang(ii,:) = w*v_ang(ii,:) + phi_p*r_p*(p_ang(ii,:) -x_ang(ii,:)) + phi_g*r_g*(g_ang-x_ang(ii,:));
        
        x_d(ii,:) = x_d(ii,:) + v_d(ii,:);
        x_amp(ii,:) = x_amp(ii,:) + v_amp(ii,:);
        x_ang(ii,:) = x_ang(ii,:) + v_ang(ii,:);
        
        
        %d_up_limit = x_d(ii,:)>d(2);
        %d_dw_limit = x_d(ii,:)<d(1);
        
        %x_d(ii,d_up_limit) = d(2);
        %x_d(ii,d_dw_limit) = d(1);
        
        %amp_up_limit = x_amp(ii,:)>amp(2);
        %amp_dw_limit = x_amp(ii,:)<amp(1);
        
        %x_amp(ii,amp_up_limit) = amp(2);
        %x_amp(ii,amp_dw_limit) = amp(1);
        
        
        distance(1) = 0;
        
        for indx_d = 1:length(x_d(ii,:))
            distance(indx_d+1) = x_d(ii,indx_d) + distance(indx_d);
        end
        
        
        
        [AF fitness]=compute_fitness(N,distance,x_amp(ii,:),x_ang(ii,:),lambda,fitness_condition,Main_beam_dir,nulls_locations,beamwidth);
        
        if (fitness < fitness_p(ii))
            
            p_d(ii,:) = x_d(ii,:);
            p_amp(ii,:) = x_amp(ii,:);
            p_ang(ii,:) = x_ang(ii,:);
            fitness_p(ii) = fitness;
            
            if (fitness < fitness_g);
                
                g_d = p_d(ii,:);
                g_amp = p_amp(ii,:);
                g_ang = p_ang(ii,:);
                fitness_g = fitness;
                
            end
            
        end
        
        clear fitness;
        clear distance;
        clear r_p
        clear r_g;
        
    end
    
    fprintf('Itteration %d of %d Fitness value is %d \n',outer_itt, Itt_max,-10*log10(-fitness_g));
    fitness1(outer_itt,:)=-10*log10(-fitness_g);
    
end

distance(1) = 0;  %% The first element is at the distance 0

for indx_d = 1:length(g_d)
    distance(indx_d+1) = g_d(indx_d) + distance(indx_d);
end



[AF fitness_check]=compute_fitness(N,distance,g_amp,g_ang,lambda,fitness_condition,Main_beam_dir,nulls_locations,beamwidth);


AF_power = abs(AF).^2;                                      %% Computing power in the array factor expression
AF_plot = AF_power./max(AF_power);                          %% Normalize with respect to maximum power


PATTERN  = 10*log10(AF_plot)
[Max_PATTERN Indx_Max ] = max(PATTERN)% max value in pattern is in Max PATTERN and Indx_Max contains its index.
PATTERN_HPBW = PATTERN>(Max_PATTERN-3); % look at angles at which PATTERN IS ABOVE -3 dB of max value
fprintf('HALF POWER BEAM WIDTH IS %d',sum(PATTERN_HPBW));

Angle_first_null=0;
AF_Processed = 10*log10(AF_power);
for ang_index = 89:-1:2
    if AF_Processed(ang_index)>AF_Processed(ang_index+1)
        Angle_first_null = ang_index;
        break;
    end
end
SLL = max(AF_Processed) - max(AF_Processed(1:Angle_first_null+1))

%   keyboard

it=1:Itt_max;
figure(1),clf
plot(it,-fitness1)
phi = 0:1:179;
figure
if (fitness_condition==0)
    
    plot(phi,10*log10(AF_plot));                            %% plotting the power in array factor in dB salce
    hold on
    axis([0 180 -110 0]);
    
elseif (fitness_condition==1)
    
    plot(phi,10*log10(AF_plot));                            %% plotting the power in array factor in dB salce
    hold on
    axis([0 180 -60 0]);
    
elseif (fitness_condition==2)
    figure(2),clf
    plot(phi,10*log10(AF_plot));                            %% plotting the power in array factor in dB salce
    hold on
%     axis([0 180 -80 0]);
    
elseif (fitness_condition==3)
    
    plot(phi,10*log10(AF_plot));                            %% plotting the power in array factor in dB salce
    hold on
    axis([0 180 -60 0]);
    
elseif (fitness_condition==4)
    
    plot(phi,10*log10(AF_plot));                            %% plotting the power in array factor in dB salce
    hold on
    axis([0 180 -60 0]);
    
elseif (fitness_condition==5)
    
    plot(phi,10*log10(AF_plot));                            %% plotting the power in array factor in dB salce
    hold on
    axis([0 180 -60 0]);
end

% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% else
%     fprintf('You have stopped the program \n');
%
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g_d = g_d';
g_ang = g_ang';
g_amp = g_amp'
a=20*log10(AF_plot);
b=(a/2)';
%  SLL=(-(AF_power(Main_beam_dir+1) / sum(AF_power([1: Main_beam_dir-beamwidth  Main_beam_dir+beamwidth:180]))));
% SLL=SLL'




