function [AF, fitness]=compute_fitness(N,distance,amplitude,phase,lambda,fitness_condition,Main_beam_dir,nulls_locations,beamwidth )


k  = 2*pi/lambda;
phi = 0:1:179;               %% phi angle range over which Array factor is found
AF = zeros(1,length(phi));

%    distance =[0 0.5 1 1.5 2 2.5 3 3.5].*lambda;            %%%% distance fixed for N=8 %%%
%    distance =[0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5].*lambda;     %%%% distance fixed for N=16 %%%
%    distance =[0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 ].*lambda;     %%%% distance fixed for N=12 %%%
   distance =[0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 10.5 11 11.5 12 12.5 13 13.5 14 14.5 15 15.5].*lambda;  %%%% distance fixed for N=32 %%%
%   distance =[0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 10.5 11 11.5 12 12.5 13 13.5 14 14.5 15 15.5 16 16.5 17 17.5 18 18.5 19 19.5 20 20.5 21 21.5 22 22.5 23 23.5 24 24.5 25 25.5 26 26.5 27 27.5 28 28.5 29 29.5 30 30.5 31 31.5].*lambda;  %%%% distance fixed for N=64 %%%
%distance = (0:0.5:31.5).*lambda

% amplitude=ones(1,N);
%       amplitude = [0.5 0.4 0.5 0.6 0.7 0.8 0.9 1 .9 .8 .7 .6 .5 0.4 0.5 0.4];


%  distance=[  0     0.7421    0.5824+0.7421    0.3730+0.5824+0.7421    0.3841+0.3730+0.5824+0.7421    0.3880+0.3841+0.3730+0.5824+0.7421    0.3333+ 0.3880+0.3841+0.3730+0.5824+0.7421    0.4211+0.3333+ 0.3880+0.3841+0.3730+0.5824+0.7421    0.3878+0.4211+0.3333+ 0.3880+0.3841+0.3730+0.5824+0.7421    0.4835+0.3878+0.4211+0.3333+ 0.3880+0.3841+0.3730+0.5824+0.7421    0.4622+0.4835+0.3878+0.4211+0.3333+ 0.3880+0.3841+0.3730+0.5824+0.7421    0.6250+0.4622+0.4835+0.3878+0.4211+0.3333+ 0.3880+0.3841+0.3730+0.5824+0.7421];
% distance=[0 0.7476    0.7785+0.7476    0.6227+0.7785+0.7476    0.5390+0.6227+0.7785+0.7476    0.5375+0.5390+0.6227+0.7785+0.7476    0.4389+0.5375+0.5390+0.6227+0.7785+0.7476    0.5082+0.4389+0.5375+0.5390+0.6227+0.7785+0.7476    0.4193+0.5082+0.4389+0.5375+0.5390+0.6227+0.7785+0.7476    0.4794+0.4193+0.5082+0.4389+0.5375+0.5390+0.6227+0.7785+0.7476    0.4820+0.4794+0.4193+0.5082+0.4389+0.5375+0.5390+0.6227+0.7785+0.7476    0.5664+0.4820+0.4794+0.4193+0.5082+0.4389+0.5375+0.5390+0.6227+0.7785+0.7476    0.4457+ 0.5664+0.4820+0.4794+0.4193+0.5082+0.4389+0.5375+0.5390+0.6227+0.7785+0.7476    0.7128+0.4457+ 0.5664+0.4820+0.4794+0.4193+0.5082+0.4389+0.5375+0.5390+0.6227+0.7785+0.7476    0.7741+0.7128+0.4457+ 0.5664+0.4820+0.4794+0.4193+0.5082+0.4389+0.5375+0.5390+0.6227+0.7785+0.7476    0.8910+0.7741+0.7128+0.4457+ 0.5664+0.4820+0.4794+0.4193+0.5082+0.4389+0.5375+0.5390+0.6227+0.7785+0.7476];
phase = 0.*ones(1,N);

for n = 1:N
    
    AF = AF + amplitude(n).*exp(j*phase(n)*pi/180).*exp(j.*k.*distance(n).*cosd(phi));
    %      AF = AF + amplitude(n).*exp(j*phase(n)*pi/180).*exp(j.*k.*distance.*cosd(phi));
end

AF_power = abs(AF).^2;



AF_plot = AF_power./max(AF_power);                          %% Normalize with respect to maximum power
PATTERN  = 10*log10(AF_plot);
[Max_PATTERN Indx_Max ] = max(PATTERN); % max value in pattern is in Max PATTERN and Indx_Max contains its index.
PATTERN_HPBW = PATTERN>(Max_PATTERN-3); % look at angles at which PATTERN IS ABOVE -3 dB of max value
HPBW=sum(PATTERN_HPBW);

Angle_first_null=0;
AF_Processed = 10*log10(AF_power);
for ang_index = 89:-1:2
    if AF_Processed(ang_index)>AF_Processed(ang_index+1)
        Angle_first_null = ang_index;
        break;
    end
end
SLL = max(AF_Processed) - max(AF_Processed(1:Angle_first_null+1));

if (fitness_condition ==0)  %%%%%%%%%%% my fitness %%%%%%%%%%%%
    fitness=max(AF_power(nulls_locations)/max(AF_power));
    %--------------------------------For Side Lobe Level-------------------
    if (fitness_condition ==4)  %%%%%%%%%%% my fitness %%%%%%%%%%%%
                                 Angle_first_null=0;
                                 AF_Processed = 10*log10(AF_power);
                                 for ang_index = 89:-1:2
                                     if AF_Processed(ang_index)>AF_Processed(ang_index+1)
                                         Angle_first_null = ang_index;
                                         break;
                                     end
                                 end
                                 SLL = max(AF_Processed) - max(AF_Processed(1:Angle_first_null+1));
                                 
                                 %-----------------------------Fitness function-----------------------
                            %     
                            %     
                            %     A=1;                   %weightin coeficent
                            %     B=.5;
                            %     for phi=0:1:90
                            %         
                            %         
                            %         fitness=(A*20*log(abs(max(AF_power)))-B*20*log(abs(max(SLL))))/(min(abs(AF_power(90))));
                            %     end
                            %     
                            %     for phi=90:1:180
                            %         
                            %         
                            %         fitness=(A*20*log(abs(max(AF_power)))-B*20*log(abs(max(SLL))))/(min(abs(AF_power(90))));
                            %     end
                                
                                %
                            elseif (fitness_condition ==1)  %%%%%%%%%%% my fitness 2 %%%%%%%%%%%%
                                
                                AF_plot = AF_power./max(AF_power);                          %% Normalize with respect to maximum power
                                PATTERN  = 10*log10(AF_plot);
                                [Max_PATTERN Indx_Max ] = max(PATTERN); % max value in pattern is in Max PATTERN and Indx_Max contains its index.
                                PATTERN_HPBW = PATTERN>(Max_PATTERN-3); % look at angles at which PATTERN IS ABOVE -3 dB of max value
                                HPBW=sum(PATTERN_HPBW);
                                
                                Angle_first_null=0;
                                AF_Processed = 10*log10(AF_power);
                                for ang_index = 89:-1:2
                                    if AF_Processed(ang_index)>AF_Processed(ang_index+1)
                                        Angle_first_null = ang_index;
                                        break;
                                    end
                                end
                                SLL = max(AF_Processed) - max(AF_Processed(1:Angle_first_null+1));
                                
                                if HPBW<=10&&SLL>16
                                    
                                    fitness=(max(HPBW))-max(SLL);
                                    
                                else
                                    fitness=(max(HPBW-10)-max(SLL-16));
                                end
                                
                                             fitness=fitness/max(abs(AF_power));
                                
                                
                                
                                
                            elseif (fitness_condition ==2)
                                     AF_plot = AF_power./max(AF_power);                          %% Normalize with respect to maximum power
                                         PATTERN  = 10*log10(AF_plot);
                                
                                   for phase1=1:1:60;
                                      if PATTERN(phase1)<-1;
                                knull = [50 60];
                                sidelobe1 = [1:nulls_locations(1)];
                                sidelobe2 = [nulls_locations(2):nullstep:180];
                                mainlobe = [nulls_locations(1):nulls_locations(2)];
                                     fitness = max(AF_power((nulls_locations+1)));
                                fitness = sum(AF_power(sidelobe1))/nullstep +
                                sum(AF_power(sidelobe2))/nullstep + sum(AF_power(mainlobe)) ;
                                fitness = max(mean(AF_power(sidelobe1))+ mean(AF_power(knull)));
                                itness = -20*log10((AF(Main_beam_dir))/max(AF(nulls_locations+1)));
                                else
                                     fitness = max(AF_power((nulls_locations1+1)));
                                end
                                fitness = max(AF_power((nulls_locations+1)));
                                     fitness=abs(fitness1);
                                     fitness=min(PATTERN)/max(PAF_power);
                                      end
                                             end
                                         end
                                   end
                                
                                     AF_plot = AF_power./max(AF_power);                          %% Normalize with respect to maximum power
                                         PATTERN  = 10*log10(AF_plot);
                                         [Max_PATTERN Indx_Max ] = max(PATTERN); % max value in pattern is in Max PATTERN and Indx_Max contains its index.
                                          PATTERN_HPBW = PATTERN>(Max_PATTERN-3); % look at angles at which PATTERN IS ABOVE -3 dB of max value
                                          HPBW=sum(PATTERN_HPBW);
                                
                                          Angle_first_null=0;
                                                 AF_Processed = 10*log10(AF_power);
                                                for ang_index = 89:-1:2
                                                     if AF_Processed(ang_index)>AF_Processed(ang_index+1)
                                                     Angle_first_null = ang_index;
                                                         break;
                                                     end
                                                end
                                                 SLL = max(AF_Processed) - max(AF_Processed(1:Angle_first_null+1));
                                
                                
                                   k=1;
                                %      for y=1:1:min(nulls_locations+1)
                                        while k<179;
                                        if phi(k)<min(nulls_locations+1)
                                
                                
                                          if HPBW<=10&&SLL>16
                                
                                            fitness=(max(HPBW))-max(SLL);
                                
                                          else
                                             fitness=(max(HPBW-10)-max(SLL-16));
                                          end
                                        end
                                if phi(k)==min(nulls_locations+1)
                                
                                     fitness = -20*log10((AF_power(Main_beam_dir))/max(AF_power(k+1)));
                                end
                                if phi(k)>min(nulls_locations+1)
                                
                                    if HPBW<=10&&SLL>16
                                
                                            fitness=(max(HPBW))-max(SLL);
                                
                                          else
                                             fitness=(max(HPBW-10)-max(SLL-16));
                                    end
                                end
                                
                                           k=k+1;
                                     end
                                
                                %            AF_plot = AF_power./max(AF_power);                          %% Normalize with respect to maximum power
                                %          PATTERN  = 10*log10(AF_plot);
                                %          [Max_PATTERN Indx_Max ] = max(PATTERN); % max value in pattern is in Max PATTERN and Indx_Max contains its index.
                                %           PATTERN_HPBW = PATTERN>(Max_PATTERN-3); % look at angles at which PATTERN IS ABOVE -3 dB of max value
                                %           HPBW=sum(PATTERN_HPBW);
                                %
                                %           Angle_first_null=0;
                                %                  AF_Processed = 10*log10(AF_power);
                                %                 for ang_index = 89:-1:2
                                %                      if AF_Processed(ang_index)>AF_Processed(ang_index+1)
                                %                      Angle_first_null = ang_index;
                                %                          break;
                                %                      end
                                %                 end
                                %                  SLL = max(AF_Processed) - max(AF_Processed(1:Angle_first_null+1));
                                %
                                %           if HPBW<=10&&SLL>16
                                %
                                %             fitness=(max(HPBW))-max(SLL);
                                %
                                %           else
                                %              fitness=(max(HPBW-10)-max(SLL-16));
                                %           end
                                
                                %
                                %
                                    elseif phi(k)==min(nulls_locations+1)
                                
                                       fitness=max(AF_power(nulls_locations+1));
                                
                                    elseif min(nulls_locations+1)<phi(k)<max(nulls_locations+1)
                                
                                
                                          fitness = -( sum(AF_power(Main_beam_dir-beamwidth+1:Main_beam_dir+beamwidth+1)) / sum(AF_power([1: Main_beam_dir-beamwidth  Main_beam_dir+beamwidth+2:180])));
                                
                                
                                    elseif phi(k)==max(nulls_locations+1)
                                          fitness=max(AF_power(nulls_locations+1));
                                
                                    elseif  phi(k)>max(nulls_locations+1)
                                
                                           AF_plot = AF_power./max(AF_power);                          %% Normalize with respect to maximum power
                                         PATTERN  = 10*log10(AF_plot);
                                         [Max_PATTERN Indx_Max ] = max(PATTERN); % max value in pattern is in Max PATTERN and Indx_Max contains its index.
                                          PATTERN_HPBW = PATTERN>(Max_PATTERN-3); % look at angles at which PATTERN IS ABOVE -3 dB of max value
                                          HPBW=sum(PATTERN_HPBW);
                                
                                          Angle_first_null=0;
                                                 AF_Processed = 10*log10(AF_power);
                                                for ang_index = 89:-1:2
                                                     if AF_Processed(ang_index)>AF_Processed(ang_index+1)
                                                     Angle_first_null = ang_index;
                                                         break;
                                                     end
                                                end
                                                 SLL = max(AF_Processed) - max(AF_Processed(1:Angle_first_null+1));
                                
                                          if HPBW<=10&&SLL>16
                                
                                            fitness=(max(HPBW))-max(SLL);
                                
                                          else
                                             fitness=(max(HPBW-10)-max(SLL-16));
                                          end
                                               k=k+1;
                                        end
                                     end
                            elseif (fitness_condition ==3 )  %%%%%%%%%%% my fitness 2 %%%%%%%%%%%%
                                
                                AF_plot = AF_power./max(AF_power);                          %% Normalize with respect to maximum power
                                PATTERN  = 10*log10(AF_plot);
                                [Max_PATTERN Indx_Max ] = max(PATTERN); % max value in pattern is in Max PATTERN and Indx_Max contains its index.
                                PATTERN_HPBW = PATTERN>(Max_PATTERN-3); % look at angles at which PATTERN IS ABOVE -3 dB of max value
                                HPBW=sum(PATTERN_HPBW);
                                
                                Angle_first_null=0;
                                AF_Processed = 10*log10(AF_power);
                                for ang_index = 89:-1:2
                                    if AF_Processed(ang_index)>AF_Processed(ang_index+1)
                                        Angle_first_null = ang_index;
                                        break;
                                    end
                                end
                                SLL = max(AF_Processed) - max(AF_Processed(1:Angle_first_null+1));
                                for phase=0:1:50
                                    if SLL>16
                                        
                                        fitness=min(AF_power);
                                        
                                    else
                                        fitness=min(AF_power)-16;
                                    end
                                end
end