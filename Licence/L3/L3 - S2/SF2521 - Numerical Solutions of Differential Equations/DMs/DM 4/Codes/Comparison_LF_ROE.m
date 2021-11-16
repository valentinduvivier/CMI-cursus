figure(30)
    subplot(1,3,1);
    d = 20;
    k = round(N_t/d);
    k_prime = round(iter_t/d);
    p_1 = plot(x, U_2(2*k + 1, 2:end-1), 'r-', x, U(2*k_prime + 1, 2:end-1), 'b--', 'linewidth', 2);
    axis([0 10 0.9 1.4]);
    ylabel('Wave''s height [m]'); title(['t = ',num2str(round(k*T/N_t,2)),' s']);
    legend([p_1(1), p_1(2)], {'Roe', 'Lax-Friedriech'}, 'location', 'northwest');
    
    subplot(1,3,2);
    d = 10;
    k = round(N_t/d);
    k_prime = round(iter_t/d);
    p_2 = plot(x, U_2(2*k + 1, 2:end-1), 'r-', x, U(2*k_prime + 1, 2:end-1), 'b--', 'linewidth', 2);
    axis([0 10 0.9 1.4]);
    xlabel('Position on the mesh [m]'); title(['t = ',num2str(round(k*T/N_t,2)),' s'], 'linewidth', 5);
    legend([p_2(1), p_2(2)], {'Roe', 'Lax-Friedriech'}, 'location', 'northwest');
    
    subplot(1,3,3);
    d = 7;
    k = round(N_t/d);
    k_prime = round(iter_t/d);
    p_3 = plot(x, U_2(2*k + 1, 2:end-1), 'r-', x, U(2*k_prime + 1, 2:end-1), 'b--', 'linewidth', 2);
    axis([0 10 0.9 1.4]);
    title(['t = ',num2str(round(k*T/N_t,2)),' s']);
    legend([p_3(1), p_3(2)], {'Roe', 'Lax-Friedriech'}, 'location', 'northwest');
    