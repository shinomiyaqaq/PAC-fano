function x = polarencode(u)%极化码编码
N = length(u);
n = log2(N);
x = u;
m =1; 
for d = n-1:-1:0
    for i = 1:2*m:N
        a = x(i:i+m-1); %第一部分
        b = x(i+m:i+2*m-1); % 第二部分
        x(i:i+2*m-1) = [mod(a+b,2) b]; %结合
    end
    m = m*2;%基底
end
end