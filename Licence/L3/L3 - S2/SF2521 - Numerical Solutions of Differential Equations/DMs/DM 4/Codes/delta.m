figure(130)

        k = N_t;
        color_1 = 1/15; color_2 = 1/50; color_3 = 1/15;
       
%         plot(x, U_delta(2, 2:end-1)./U_delta(1, 2:end-1), 'b');
        hold on;
        plot(x, U_delta(4, 2:end-1)./U_delta(3, 2:end-1), 'm');
        plot(x, U_delta(6, 2:end-1)./U_delta(5, 2:end-1), 'r');
        plot(x, U_delta(8, 2:end-1)./U_delta(7, 2:end-1), 'b');
        plot(x, U_delta(10, 2:end-1)./U_delta(9, 2:end-1), 'c');
        
        plot(x, U_delta(12, 2:end-1)./U_delta(11, 2:end-1), 'g');
        plot(x, U_delta(14, 2:end-1)./U_delta(13, 2:end-1), 'y');
        plot(x, U_delta(16, 2:end-1)./U_delta(15, 2:end-1), 'b');
        plot(x, U_delta(18, 2:end-1)./U_delta(17, 2:end-1), 'g');
        plot(x, U_delta(20, 2:end-1)./U_delta(19, 2:end-1), 'r');
        
        plot(x, U_delta(22, 2:end-1)./U_delta(21, 2:end-1), 'k');
%         plot(x, U_delta(24, 2:end-1)./U_delta(23, 2:end-1), 'color', 12*[color_1 color_2 color_3]);
%         plot(x, U_delta(26, 2:end-1)./U_delta(25, 2:end-1), 'color', 13*[color_1 color_2 color_3]);
%         plot(x, U_delta(28, 2:end-1)./U_delta(27, 2:end-1), 'color', 14*[color_1 color_2 color_3]);
%         plot(x, U_delta(30, 2:end-1)./U_delta(29, 2:end-1), 'color', 15*[color_1 color_2 color_3]);
        axis([4 6 2 3.5]);
%         axis([0 10 1.5 4]);
        xlabel('Position on the mesh [m]'); ylabel('Wave''s speed [m.s^{-1}]'); title(['t = ',num2str(round(k*T/N_t,2)),' s']);
        grid on;
        legend({['\delta = ',num2str(delta(1)),''], ['\delta = ',num2str(delta(2)),''], ['\delta = ',num2str(delta(3)),''], ...
            ['\delta = ',num2str(delta(4)),''], ['\delta = ',num2str(delta(5)),''], ['\delta = ',num2str(delta(6)),''], ['\delta = ',num2str(delta(7)),''], ...
            ['\delta = ',num2str(delta(8)),''], ['\delta = ',num2str(delta(9)),''], ['\delta = ',num2str(delta(10)),'']}, 'location', 'northwest');
