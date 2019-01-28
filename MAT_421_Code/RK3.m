function RK3= RK3(time_step,y0)

   t0=0;    %start time
   tf=1;    %final Time
   dt=time_step;    %time step
   y=y0;    %initial condition
   ysol=[y];% approximate solution vector
   time=[t0];   %time vector
   while t0<=tf %iterator for RK-3
    temp_t=t0;  
    k1=f(temp_t,y);
    k2=f(temp_t+(dt/2),y+(0.5*dt*k1));
    k3=f(temp_t+dt,y+dt*(2*k2-k1));
    y=y+(dt*(k1+4*k2+k3)/6);
    ysol=[ysol y];
    t0=t0+dt;
    time=[time t0];
    
   end
    exactsol=g(time);
    hold on
    plot(time,ysol,'r-','MarkerSize',24,'LineWidth',2);
    title('RK-3 method for yprime=y');
    xlabel('Time','FontSize',24);
    ylabel('Solution','FontSize',24);
    legend('RK-3');
    hold off
end    
    
function yprime=f(t,y)
    yprime=y;
end

function exact=g(t)
    exact=exp(t);
end