%% Solve the non-linear BVP
%%-------------------------
%%   y?? = y^3-yy?
%%   y(1)=1/2, y(2)=1/3
f = @(x,z) [z(2);z(1)^3-z(1)*z(2)]; 
a = 1; alpha = 1/2;
b=2;beta =1/3;
%% Parameters
s0 = 0;
s1 = 1;
N = 1000;
nbIter = 5;
%% save all?
saveAllfunctions = 1;
%% Initialization
[z,x] = RK4(f,[a,b],[alpha s0],1/N);
y_b_km1= z(1,end);
s_k  = s1;
s_km1= s0;

if (saveAllfunctions==1)
  stockY = zeros(nbIter+1,N+1);
  stockS = zeros(1,nbIter+1);
  stockY(1,:) = z(1,:);
  stockS(1)   = s0;
end
%% error
yExact = @(x) 1./(1+x);
stockError = zeros(1,nbIter+1);
stockError(1) = max(z(1,:)-yExact(x));
%%---------------------------%% %%----- loop ---%% %%---------------------------%% 
for k=1:nbIter
  %% Solve a new IVP
[z,x] = RK4(f,a,[alpha,s_k],b,1/N);
y_b_k = z(1,end);
%% the new s (secant method)
s = s_k - (y_b_k - beta)/(y_b_k-y_b_km1)*(s_k-s_km1); %% save
if (saveAllfunctions==1)
    stockY(k+1,:) = z(1,:); 
    stockS(k+1) = s_k;
end
stockError(k+1) = max(abs(z(1,:)-yExact(x))); %% update
s_km1 = s_k;
s_k = s;
 y_b_km1 = y_b_k;
end
%%---------------------------%%
%%---- plot 
plot(x,stockY(3,:),'linewidth',5,x,yExact(x),'linewidth',5)
xlabel('x')
ylabel('y')
legend('y numeric','y exact','location','northeast')
title('y''=y^3-yy, y(1)=1/2,y(2)=1/3')