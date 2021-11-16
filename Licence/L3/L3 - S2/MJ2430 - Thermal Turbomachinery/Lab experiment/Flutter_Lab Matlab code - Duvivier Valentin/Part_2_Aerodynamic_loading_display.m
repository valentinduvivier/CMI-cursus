%% Aerodynamic loading display : B-2

% The pressure on each point along the distance gives the loading
% We consider below the blade B-2
    Pression_Bm2_a10_f10 = A_a10_f10([1:15],1);
    Pression_Bm2_a10_f15 = A_a10_f15([1:15],1);
    Pression_Bm2_a10_f20 = A_a10_f20([1:15],1);
    Pression_Bm2_a20_f10 = A_a20_f10([1:15],1);
    Pression_Bm2_a20_f15 = A_a20_f15([1:15],1);
    Pression_Bm2_a20_f20 = A_a20_f20([1:15],1);
    Pression_Bm2_a30_f10 = A_a30_f10([1:15],1);
    Pression_Bm2_a30_f15 = A_a30_f15([1:15],1);
    Pression_Bm2_a30_f20 = A_a30_f20([1:15],1);
   
    figure(10)
    p1 = plot(coordinates, Pression_Bm2_a10_f10, 'r-');
    hold on
    plot(coordinates, Pression_Bm2_a10_f15, 'r-');
    plot(coordinates, Pression_Bm2_a10_f20, 'r-');
    p2 = plot(coordinates, Pression_Bm2_a20_f10, 'g-');
    plot(coordinates, Pression_Bm2_a20_f15, 'g-');
    plot(coordinates, Pression_Bm2_a20_f20, 'g-');
    p3 = plot(coordinates, Pression_Bm2_a30_f10, 'y*-');
    plot(coordinates, Pression_Bm2_a30_f15, 'y*-');
    plot(coordinates, Pression_Bm2_a30_f20, 'y*-');
    hold off
    legend([p1 p2 p3],{'Angle = 10°','Angle = 20°','Angle = 30°'},'Location','northwest');
    grid on;
    title('Aerodynamic loading for blade b-2, each 15 sensors');

    
%% Aerodynamic loading display : B-1

% The pressure on each point along the distance gives the loading
% We consider below the blade B-1
    Pression_Bm1_a10_f10 = A_a10_f10([16:30],1);
    Pression_Bm1_a10_f15 = A_a10_f15([16:30],1);
    Pression_Bm1_a10_f20 = A_a10_f20([16:30],1);
    Pression_Bm1_a20_f10 = A_a20_f10([16:30],1);
    Pression_Bm1_a20_f15 = A_a20_f15([16:30],1);
    Pression_Bm1_a20_f20 = A_a20_f20([16:30],1);
    Pression_Bm1_a30_f10 = A_a30_f10([16:30],1);
    Pression_Bm1_a30_f15 = A_a30_f15([16:30],1);
    Pression_Bm1_a30_f20 = A_a30_f20([16:30],1);
   
    figure(11)
    p1 = plot(coordinates, Pression_Bm1_a10_f10, 'r-');
    hold on
    plot(coordinates, Pression_Bm1_a10_f15, 'r-');
    plot(coordinates, Pression_Bm1_a10_f20, 'r-');
    p2 = plot(coordinates, Pression_Bm1_a20_f10, 'g-');
    plot(coordinates, Pression_Bm1_a20_f15, 'g-');
    plot(coordinates, Pression_Bm1_a20_f20, 'g-');
    p3 = plot(coordinates, Pression_Bm1_a30_f10, 'y*-');
    plot(coordinates, Pression_Bm1_a30_f15, 'y*-');
    plot(coordinates, Pression_Bm1_a30_f20, 'y*-');
    hold off
    legend([p1 p2 p3],{'Angle = 10°','Angle = 20°','Angle = 30°'},'Location','northwest');
    grid on;
    title('Aerodynamic loading for blade b-1, each 15 sensors');

    
%% Aerodynamic loading display : B0

% The pressure on each point along the distance gives the loading
% We consider below the blade B0
    Pression_B0_a10_f10 = A_a10_f10([31:45],1);
    Pression_B0_a10_f15 = A_a10_f15([31:45],1);
    Pression_B0_a10_f20 = A_a10_f20([31:45],1);
    Pression_B0_a20_f10 = A_a20_f10([31:45],1);
    Pression_B0_a20_f15 = A_a20_f15([31:45],1);
    Pression_B0_a20_f20 = A_a20_f20([31:45],1);
    Pression_B0_a30_f10 = A_a30_f10([31:45],1);
    Pression_B0_a30_f15 = A_a30_f15([31:45],1);
    Pression_B0_a30_f20 = A_a30_f20([31:45],1);
    
    figure(12)
    p1 = plot(coordinates, Pression_B0_a10_f10, 'r-');
    hold on
    plot(coordinates, Pression_B0_a10_f15, 'r-');
    plot(coordinates, Pression_B0_a10_f20, 'r-');
    p2 = plot(coordinates, Pression_B0_a20_f10, 'g-');
    plot(coordinates, Pression_B0_a20_f15, 'g-');
    plot(coordinates, Pression_B0_a20_f20, 'g-');
    p3 = plot(coordinates, Pression_B0_a30_f10, 'y*-');
    plot(coordinates, Pression_B0_a30_f15, 'y*-');
    plot(coordinates, Pression_B0_a30_f20, 'y*-');
    hold off
    legend([p1 p2 p3],{'Angle = 10°','Angle = 20°','Angle = 30°'},'Location','northwest');
    grid on;
    title('Aerodynamic loading - blade b0, each 15 sensors');
    xlabel('normalised position x'); ylabel('Loading');

%% Aerodynamic loading display : B+1

% The pressure on each point along the distance gives the loading
% We consider below the blade B+1
    Pression_Bp1_a10_f10 = A_a10_f10([46:60],1);
    Pression_Bp1_a10_f15 = A_a10_f15([46:60],1);
    Pression_Bp1_a10_f20 = A_a10_f20([46:60],1);
    Pression_Bp1_a20_f10 = A_a20_f10([46:60],1);
    Pression_Bp1_a20_f15 = A_a20_f15([46:60],1);
    Pression_Bp1_a20_f20 = A_a20_f20([46:60],1);
    Pression_Bp1_a30_f10 = A_a30_f10([46:60],1);
    Pression_Bp1_a30_f15 = A_a30_f15([46:60],1);
    Pression_Bp1_a30_f20 = A_a30_f20([46:60],1);
    
    figure(13)
    p1 = plot(coordinates, Pression_Bp1_a10_f10, 'r-');
    hold on
    plot(coordinates, Pression_Bp1_a10_f15, 'r-');
    plot(coordinates, Pression_Bp1_a10_f20, 'r-');
    p2 = plot(coordinates, Pression_Bp1_a20_f10, 'g-');
    plot(coordinates, Pression_Bp1_a20_f15, 'g-');
    plot(coordinates, Pression_Bp1_a20_f20, 'g-');
    p3 = plot(coordinates, Pression_Bp1_a30_f10, 'y*-');
    plot(coordinates, Pression_Bp1_a30_f15, 'y*-');
    plot(coordinates, Pression_Bp1_a30_f20, 'y*-');
    hold off
    legend([p1 p2 p3],{'Angle = 10°','Angle = 20°','Angle = 30°'},'Location','northwest');
    grid on;
    title('Aerodynamic loading for blade b+1, each 15 sensors');

    
%% Aerodynamic loading display : B+2

% The pressure on each point along the distance gives the loading
% We consider below the blade B+2
    Pression_Bp2_a10_f10 = A_a10_f10([61:75],1);
    Pression_Bp2_a10_f15 = A_a10_f15([61:75],1);
    Pression_Bp2_a10_f20 = A_a10_f20([61:75],1);
    Pression_Bp2_a20_f10 = A_a20_f10([61:75],1);
    Pression_Bp2_a20_f15 = A_a20_f15([61:75],1);
    Pression_Bp2_a20_f20 = A_a20_f20([61:75],1);
    Pression_Bp2_a30_f10 = A_a30_f10([61:75],1);
    Pression_Bp2_a30_f15 = A_a30_f15([61:75],1);
    Pression_Bp2_a30_f20 = A_a30_f20([61:75],1);
   
    figure(14)
    p1 = plot(coordinates, Pression_Bp2_a10_f10, 'r-');
    hold on
    plot(coordinates, Pression_Bp2_a10_f15, 'r-');
    plot(coordinates, Pression_Bp2_a10_f20, 'r-');
    p2 = plot(coordinates, Pression_Bp2_a20_f10, 'g-');
    plot(coordinates, Pression_Bp2_a20_f15, 'g-');
    plot(coordinates, Pression_Bp2_a20_f20, 'g-');
    p3 = plot(coordinates, Pression_Bp2_a30_f10, 'y*-');
    plot(coordinates, Pression_Bp2_a30_f15, 'y*-');
    plot(coordinates, Pression_Bp2_a30_f20, 'y*-');
    hold off
    legend([p1 p2 p3],{'Angle = 10°','Angle = 20°','Angle = 30°'},'Location','northwest');
    grid on;
    title('Aerodynamic loading for blade b+2, each 15 sensors');

    
    
%% Display with subplot : B-2 & B-1 & B0 & B+1 & B+2    

figure(15)

subplot(2,3,1)
    p1 = plot(coordinates, Pression_Bm2_a10_f10, 'r-');
    hold on
    plot(coordinates, Pression_Bm2_a10_f15, 'r-');
    plot(coordinates, Pression_Bm2_a10_f20, 'r-');
    p2 = plot(coordinates, Pression_Bm2_a20_f10, 'g-');
    plot(coordinates, Pression_Bm2_a20_f15, 'g-');
    plot(coordinates, Pression_Bm2_a20_f20, 'g-');
    p3 = plot(coordinates, Pression_Bm2_a30_f10, 'y*-');
    plot(coordinates, Pression_Bm2_a30_f15, 'y*-');
    plot(coordinates, Pression_Bm2_a30_f20, 'y*-');
    hold off
    legend([p1 p2 p3],{'Angle = 10°','Angle = 20°','Angle = 30°'},'Location','northwest');
    grid on;
    title('Aerodynamic loading - blade b-2');
    xlabel('normalised position x'); ylabel('Loading');
    axis([-0.6 0.6 98 101.5]);
    
 subplot(2,3,2)
    p1 = plot(coordinates, Pression_Bm1_a10_f10, 'r-');
    hold on
    plot(coordinates, Pression_Bm1_a10_f15, 'r-');
    plot(coordinates, Pression_Bm1_a10_f20, 'r-');
    p2 = plot(coordinates, Pression_Bm1_a20_f10, 'g-');
    plot(coordinates, Pression_Bm1_a20_f15, 'g-');
    plot(coordinates, Pression_Bm1_a20_f20, 'g-');
    p3 = plot(coordinates, Pression_Bm1_a30_f10, 'y*-');
    plot(coordinates, Pression_Bm1_a30_f15, 'y*-');
    plot(coordinates, Pression_Bm1_a30_f20, 'y*-');
    hold off
    legend([p1 p2 p3],{'Angle = 10°','Angle = 20°','Angle = 30°'},'Location','northwest');
    grid on;
    title('Aerodynamic loading - blade b-1');
    xlabel('normalised position x'); ylabel('Loading');
    axis([-0.6 0.6 98 101.5]);
    
subplot(2,3,3)
    p1 = plot(coordinates, Pression_B0_a10_f10, 'r-');
    hold on
    plot(coordinates, Pression_B0_a10_f15, 'r-');
    plot(coordinates, Pression_B0_a10_f20, 'r-');
    p2 = plot(coordinates, Pression_B0_a20_f10, 'g-');
    plot(coordinates, Pression_B0_a20_f15, 'g-');
    plot(coordinates, Pression_B0_a20_f20, 'g-');
    p3 = plot(coordinates, Pression_B0_a30_f10, 'y*-');
    plot(coordinates, Pression_B0_a30_f15, 'y*-');
    plot(coordinates, Pression_B0_a30_f20, 'y*-');
    hold off
    legend([p1 p2 p3],{'Angle = 10°','Angle = 20°','Angle = 30°'},'Location','northwest');
    grid on;
    title('Aerodynamic loading - blade b0');
    xlabel('normalised position x'); ylabel('Loading');
    axis([-0.6 0.6 98 101.5]);
    
subplot(2,3,4)
    p1 = plot(coordinates, Pression_Bp1_a10_f10, 'r-');
    hold on
    plot(coordinates, Pression_Bp1_a10_f15, 'r-');
    plot(coordinates, Pression_Bp1_a10_f20, 'r-');
    p2 = plot(coordinates, Pression_Bp1_a20_f10, 'g-');
    plot(coordinates, Pression_Bp1_a20_f15, 'g-');
    plot(coordinates, Pression_Bp1_a20_f20, 'g-');
    p3 = plot(coordinates, Pression_Bp1_a30_f10, 'y*-');
    plot(coordinates, Pression_Bp1_a30_f15, 'y*-');
    plot(coordinates, Pression_Bp1_a30_f20, 'y*-');
    hold off
    legend([p1 p2 p3],{'Angle = 10°','Angle = 20°','Angle = 30°'},'Location','northwest');
    grid on;
    title('Aerodynamic loading - blade b+1');
    xlabel('normalised position x'); ylabel('Loading');
    axis([-0.6 0.6 98 101.5]);

subplot(2,3,5)
    p1 = plot(coordinates, Pression_Bp2_a10_f10, 'r-');
    hold on
    plot(coordinates, Pression_Bp2_a10_f15, 'r-');
    plot(coordinates, Pression_Bp2_a10_f20, 'r-');
    p2 = plot(coordinates, Pression_Bp2_a20_f10, 'g-');
    plot(coordinates, Pression_Bp2_a20_f15, 'g-');
    plot(coordinates, Pression_Bp2_a20_f20, 'g-');
    p3 = plot(coordinates, Pression_Bp2_a30_f10, 'y*-');
    plot(coordinates, Pression_Bp2_a30_f15, 'y*-');
    plot(coordinates, Pression_Bp2_a30_f20, 'y*-');
    hold off
    legend([p1 p2 p3],{'Angle = 10°','Angle = 20°','Angle = 30°'},'Location','northwest');
    grid on;
    title('Aerodynamic loading - blade b+2');
    xlabel('normalised position x'); ylabel('Loading');
    axis([-0.6 0.6 98 101.5]);

    