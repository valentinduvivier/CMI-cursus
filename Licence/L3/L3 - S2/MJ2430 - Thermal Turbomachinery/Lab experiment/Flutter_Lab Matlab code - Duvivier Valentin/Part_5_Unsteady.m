%% To do this part : 
%   - be sure that the fils les_bugges.mat is loaded
%   - be sure to run the Part_1 if not already done 

%% Steady Part    
    Pression_B0_Unsteady = Flutter_B0([31:45],1);
    
%% Unsteady Part

    % Time vector
    time = x_46(:,1)';

    % We plot the corrected pressures for each 15 taps
    figure(30)
    for k = 1:size(p_corr_b0,2)
        plot(time, p_corr_b0(:,k));
        hold on
    end
    grid on
    hold off
    title('Unsteady corrected loading ; blade b0'); xlabel('time [s]'); ylabel('Pressure [kPa]');

%% Work on logical functions - Creation errors bands

        % We only keep the value ofthe pressure such that we don't take
        % into account the flutter part
            p1 = p_corr_b0(:,1);        p_no_flutter_1 = p1((0.5 < time) & (time <= 3.5));
            p2 = p_corr_b0(:,2);        p_no_flutter_2 = p2((0.5 < time) & (time <= 3.5));
            p3 = p_corr_b0(:,3);        p_no_flutter_3 = p3((0.5 < time) & (time <= 3.5));
            p4 = p_corr_b0(:,4);        p_no_flutter_4 = p4((0.5 < time) & (time <= 3.5));
            p5 = p_corr_b0(:,5);        p_no_flutter_5 = p5((0.5 < time) & (time <= 3.5));
            p6 = p_corr_b0(:,6);        p_no_flutter_6 = p6((0.5 < time) & (time <= 3.5));
            p7 = p_corr_b0(:,7);        p_no_flutter_7 = p7((0.5 < time) & (time <= 3.5));
            p8 = p_corr_b0(:,8);        p_no_flutter_8 = p8((0.5 < time) & (time <= 3.5));
            p9 = p_corr_b0(:,9);        p_no_flutter_9 = p9((0.5 < time) & (time <= 3.5));
            p10 = p_corr_b0(:,10);      p_no_flutter_10 = p10((0.5 < time) & (time <= 3.5));
            p11 = p_corr_b0(:,11);      p_no_flutter_11 = p11((0.5 < time) & (time <= 3.5));
            p12 = p_corr_b0(:,12);      p_no_flutter_12 = p12((0.5 < time) & (time <= 3.5));
            p13 = p_corr_b0(:,13);      p_no_flutter_13 = p13((0.5 < time) & (time <= 3.5));
            p14 = p_corr_b0(:,14);      p_no_flutter_14 = p14((0.5 < time) & (time <= 3.5));
            p15 = p_corr_b0(:,15);      p_no_flutter_15 = p15((0.5 < time) & (time <= 3.5));
            
            p_no_flutter = [p_no_flutter_1, p_no_flutter_2, p_no_flutter_3, p_no_flutter_4, p_no_flutter_5, p_no_flutter_6, p_no_flutter_7, p_no_flutter_8, p_no_flutter_9, p_no_flutter_10, p_no_flutter_11, p_no_flutter_12, p_no_flutter_13, p_no_flutter_14, p_no_flutter_15];
        
        %%
        % We plot for the values of pressure that are not in flutter
        % We check if everything went well
        figure(31)
        for k = 1:size(p_no_flutter,2)
            t(k) = plot(time((0.5<time)&(time<=3.5)), p_no_flutter(:,k));
            hold on
        end
        hold off
        legend([t],{'Point1', 'Point2', 'Point3', 'Point4', 'Point5', 'Point6', 'Point7', 'Point8', 'Point9', 'Point10', 'Point11', 'Point12', 'Point13', 'Point14', 'Point15'}, 'location', 'northwest');
        title('Unsteady pressure for each point for  0.5 < time < 3.5');
        xlabel('time'); ylabel('Pressure');
            
        %%
        % Mean value for the pressure
            
            f = zeros(1,size(p_no_flutter,2));   % Sum pressure

            for j = 1:size(p_no_flutter,2)
                r = 0;  % Nb iteration
                for i = 1:size(p_no_flutter,1)
                    f(1,j) = f(1,j) + p_no_flutter(i,j);
                    r = r + 1;
                end
            end
            
            mean_pressure = f/r;
            
            %%
            
        % We search for the max and min of the pressure
        % These values will give us the extremums that can be taken by our
        % variations in the no-flutter part and to hus be able to construct the
        % error bands
            min_p_no_flutter = min(p_no_flutter);
            max_p_no_flutter = max(p_no_flutter);
     
        % Interval of fluctuation :
            P_B0_Unsteady = Pression_B0_Unsteady(1);
                for k = 1:size(p_no_flutter,2)-1
                    P_B0_Unsteady = [P_B0_Unsteady Pression_B0_Unsteady(1)];
                end
        
            Interval_fluctuation = [P_B0_Unsteady - min_p_no_flutter; P_B0_Unsteady + max_p_no_flutter];
        
            % We reorganise the Interval of fluctuation so that the min are
            % above the max (see the matrix)
            for j = 1:15
                for i = 1
                    if Interval_fluctuation(i,j) > Interval_fluctuation(i+1,j)
                        c = Interval_fluctuation(i,j);
                        Interval_fluctuation(i,j) = Interval_fluctuation(i+1,j);
                        Interval_fluctuation(i+1,j) = c;
                    end
                end
            end
%%
        figure(32)
        plot(coordinates, Pression_B0_Unsteady);
        hold on
        for k = 1:size(coordinates,1)
            p(k) = plot([coordinates(k) coordinates(k)], Interval_fluctuation(:,k));
        end
        hold off
        legend([p],{'Point1', 'Point2', 'Point3', 'Point4', 'Point5', 'Point6', 'Point7', 'Point8', 'Point9', 'Point10', 'Point11', 'Point12', 'Point13', 'Point14', 'Point15'},'Location','northeast');
        xlabel('cooordinates'); ylabel('pressure [kPa]'); title('Shape of the loading while considering range of fluctuation due to unsteady loading');

%% Display of variation : example point 1

    % Multiplicator to adapt the variations
        multiplicator_positive =  1.5;
        multiplicator_negative =  0.67;


   % Offset : we move the variation around the considered pressure of
   % reference ( = steady one)
        for k = 1:size(p_corr_b0,1)
             p_corr_b0(k,1) =  p_corr_b0(k,1) + (Pression_B0_Unsteady(1) - mean_pressure(1));
        end
   
   % Adaptation of the range of uncertainty
       for k = 1:size(p_corr_b0,1)
           if (p_corr_b0(k,1) > Pression_B0_Unsteady(1)) 
               % What's above the mean_pressure increases
               p_corr_b0(k,1) = p_corr_b0(k,1)*multiplicator_positive;
           else if (p_corr_b0(k,1) < Pression_B0_Unsteady(1))
               % What's below the mean_pressure decreases
                p_corr_b0(k,1) = p_corr_b0(k,1)*multiplicator_negative;
               else
                % What's on the mean_pressure keep constant
                p_corr_b0(k,1) = p_corr_b0(k,1);
               end
           end
       end

    % Combination : we display the exemple results
        P_B0_Unsteady = Pression_B0_Unsteady(1);
        for k = 1:size(time,2)-1
            P_B0_Unsteady = [P_B0_Unsteady Pression_B0_Unsteady(1)];
        end
        
        figure(33)
        p1 = plot(time, P_B0_Unsteady, 'r-*');
        hold on
        p2 = plot(time, p_corr_b0(:,1), 'g-');
        hold off
        legend([p1 p2],{'Steady pressure','Unsteady variations'},'Location','northwest');
        grid on
        xlabel('time [s]'); ylabel('Pressure [kPa]'); title('Pressure point 1, blade 0 considering unsteady variations');

%     plot(coordinates, Pression_B0_a10_f15, 'r-');
%     plot(coordinates, Pression_B0_a10_f20, 'r-');
%     p2 = plot(coordinates, Pression_B0_a20_f10, 'g-');
%     plot(coordinates, Pression_B0_a20_f15, 'g-');
%     plot(coordinates, Pression_B0_a20_f20, 'g-');
%     p3 = plot(coordinates, Pression_B0_a30_f10, 'y*-');
%     plot(coordinates, Pression_B0_a30_f15, 'y*-');
%     plot(coordinates, Pression_B0_a30_f20, 'y*-');
% axis([0 3.5 90 110]);
