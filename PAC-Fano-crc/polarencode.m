function x = polarencode(u)%���������
N = length(u);
n = log2(N);
x = u;
m =1; 
for d = n-1:-1:0
    for i = 1:2*m:N
        a = x(i:i+m-1); %��һ����
        b = x(i+m:i+2*m-1); % �ڶ�����
        x(i:i+2*m-1) = [mod(a+b,2) b]; %���
    end
    m = m*2;%����
end
end