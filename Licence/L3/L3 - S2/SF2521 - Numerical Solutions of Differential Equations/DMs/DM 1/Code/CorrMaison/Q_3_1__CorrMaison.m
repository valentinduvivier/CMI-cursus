
%% Q_3.1 

% To decoment if you want to display the result of the question Q_3.1
        
        t_1 = round(N_t/100); t_2 = round(N_t/50); t_3 = round(N_t/30); t_4 = round(N_t/20);
        
        figure(2)
        HM1 = imagesc(reshape(FF(:,t_1),M,N)); 
        title(['Heat repartition, t = ',num2str(Dt*t_1),'s for S1']);
        xlabel('x'); ylabel('y');
        
        figure(3)
        HM1 = imagesc(reshape(FF(:,t_2),M,N)); 
        title(['Heat repartition, t = ',num2str(Dt*t_2),'s for S1']);
        xlabel('x'); ylabel('y');
        
        figure(4)
        HM1 = imagesc(reshape(FF(:,t_3),M,N)); 
        title(['Heat repartition, t = ',num2str(Dt*t_3),'s for S1']);
        xlabel('x'); ylabel('y');
        
        figure(5)
        HM1 = imagesc(reshape(FF(:,t_4),M,N)); 
        title(['Heat repartition, t = ',num2str(Dt*t_4),'s for S1']);
        xlabel('x'); ylabel('y');


        % Drawing

        figure(6)
            for k = 1:N_t
                HM1 = surf(reshape(FF(:,k),M,N));
%                 axis([0 M 0 N 0 1]);
                xlabel('x position'); ylabel('y position'); zlabel('Heat (°C)'); title(['Heat over the mesh in the case of source S2 at t = ',num2str(round(Dt*k,2)),' s']);
                drawnow
            end
  