maxlev = 7;
p=[ 0 1; 0 0];
%t=[ 0 1/4 3/4 1; 0 -1/4 1/4 0] 
t=[ 0 1/4 1/2 3/4 1; 0 -1/4 0 1/4 0]
figure(1)
set(gcf,'color','white')
pp = plot(p(1,:),p(2,:),'b');
axis equal
hold on
%%plot(.2*t(1,:),.2*t(2,:)+.7,'r');
hold off
axis off%%
%%pause(5)
for lev = 1:maxlev
   q = koch(lev,p,t);
   set(pp,'xdata',q(1,:),'ydata',q(2,:))
            %%pause(2)
end