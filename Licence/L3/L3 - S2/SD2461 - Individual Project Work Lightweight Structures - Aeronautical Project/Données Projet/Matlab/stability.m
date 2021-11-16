%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Stability          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_wmac=x_wac-MAC_w/4;    % x coordinate of MAC leading edge
hn=(x_wac-x_wmac)/MAC_w+V_h*a_h/a_w;
h=(x_cg-x_wmac)/MAC_w;
SM=hn-h;

figure(1)
plot3(0,-x_cg,0,'*')            % plot cg in graph
hold on
plot3(0,-hn*MAC_w,0,'*')        % plot neutral point in graph