figure(50)
    subplot(1,3,1);
    k = round(N_t/20);
    k_prime = round(iter_t/20);
    p_1 = plot(x, U_2(2*k + 1, 2:end-1), 'r', x, U_3(2*k_prime + 1, 2:end-1), 'b', 'linewidth', 2);
    axis([0 10 0.5 3]);
%   xlabel('Position on the mesh [m]'); 
    ylabel('Wave''s height [m]');     title(['t = ',num2str(round(k*T/N_t,2)),' s']);
    legend([p_1(1), p_1(2)], {'Roe (a = H/5)', 'Roe (a = 2H)'}, 'location', 'northwest');
    
    subplot(1,3,2);
    k = round(N_t/15);
    k_prime = round(iter_t/15);
    p_2 = plot(x, U_2(2*k + 1, 2:end-1), 'r', x, U_3(2*k_prime + 1, 2:end-1), 'b', 'linewidth', 2);
    axis([0 10 0.5 3]);
    xlabel('Position on the mesh [m]'); % ylabel('Wave''s height [m]');
    title(['t = ',num2str(round(k*T/N_t,2)),' s']);
    legend([p_2(1), p_2(2)], {'Roe (a = H/5)', 'Roe (a = 2H)'}, 'location', 'northwest');
    
    subplot(1,3,3);
    k = round(N_t/7);
    k_prime = round(iter_t/7);
    p_3 = plot(x, U_2(2*k + 1, 2:end-1), 'r', x, U_3(2*k_prime + 1, 2:end-1), 'b', 'linewidth', 2);
    axis([0 10 0.5 3]);
%   xlabel('Position on the mesh [m]'); ylabel('Wave''s height [m]');
    title(['t = ',num2str(round(k*T/N_t,2)),' s']);
    legend([p_3(1), p_3(2)], {'Roe (a = H/5)', 'Roe (a = 2H)'}, 'location', 'northeast');
    