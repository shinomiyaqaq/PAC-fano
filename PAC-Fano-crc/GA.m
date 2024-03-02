%高斯近似法 GA(0.7499,8)
function u = GA(sigma,N)
n=log2(N);%n次运算 列
u(1:N,1)=2/(sigma^2);
for j=1:N/2
    for i=1:n
    u(2*j-1,i+1)=phi_ni(1-(1-phi(u(j,i)))^2);
    u(2*j,i+1)=u(j,i)*2;
    end
end
u=u(:,log2(N)+1);
end
% u = bitrevorder(u);
 
function y=phi(x) %φ函数
    if x<10 && x>0
        y=exp(-0.4527*(x^0.86)+0.0218);
    else if x>=10
        y=sqrt(pi/x)*exp(-x/4)*(1-10/(7*x));    
    end
    end
end
 
%phi逆函数
function x = phi_ni(x)
%部分用闭合表达式，部分用数值解法，速度进一步提升！
    x0 = 0.0388;
    x1 = x0 - (phi(x0) - x)/derivative_phi(x0);
    delta = abs(x1 - x0);
    epsilon = 1e-3;
    
    while(delta >= epsilon)
        x0 = x1;
        x1 = x1 - (phi(x1) - x)/derivative_phi(x1);
        %当x1过大，放宽epsilon
        if x1 > 1e2
            epsilon = 10;
        end       
            delta = abs(x1 - x0);
    end
    x = x1;
end

%derivative_phi()函数
function dx = derivative_phi(x)%导数
if (x >= 0)&&(x <= 10)
%     dx = -0.4527*0.86*x^(-0.14)*phi(x);
    dx=-(194661*exp(109/5000 - (4527*x^(43/50))/10000))/(500000*x^(7/50));
else
%     dx = exp(-x/4)*sqrt(pi/x)*(-1/2/x*(1 - 10/7/x) - 1/4*(1 - 10/7/x) + 10/7/x/x);
    dx=(pi^(1/2)*exp(-x/4)*(10/(7*x) - 1)*(1/x)^(1/2))/4 + (10*pi^(1/2)*exp(-x/4)*(1/x)^(1/2))/(7*x^2) + (pi^(1/2)*exp(-x/4)*(10/(7*x) - 1))/(2*x^2*(1/x)^(1/2));
end
end
