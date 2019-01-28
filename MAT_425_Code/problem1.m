%%Problem 1
function problem1(time,y0)
    hold on
    for c=1:size(time)
        
        dt=time(1,c);
        RK3(dt,y0);
    end
    hold off
end
