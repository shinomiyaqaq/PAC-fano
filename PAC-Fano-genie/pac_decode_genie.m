function [phi,VPC,data]= ...
    pac_decode_genie(r,delta,maxcycles,poly,RP,sigma,E0)
N = size(r,2);                          % 码字长度
L=2;%两个译码器
% T=[0,-4,-6,-8];
T=[0,0];                                  % threshold 初始阈值

gama=zeros(L,N+2); %路径度量
gama(:,1)= -inf;%初始化
gama(:,2)=0;
V=zeros(L,N);%暂存结果
P=zeros(L,N);%暂存比特
VPC=zeros(1,L);%记录译码时延

l=1;%第一个译码器
[VPC(l),gama(l,:),V(l,:),T(l),P(l,:),b(l)] = individual_decoder(r,delta,maxcycles,poly,RP,sigma,E0,l, gama(l,:),V(l,:),T(l),P(l,:) );%译码完成
[~,index] = sort(gama(l,3:N+2));
T(l+1)= floor(gama(index(1))/3/delta) *3*delta;%选出第二个译码器的初始T

% clc
% floor(gama(index(1))/delta) *delta     %Fano-genie
% floor(gama(index(1))/2/delta) *2*delta  %Fano-genie2
% floor(gama(index(1))/3/delta) *3*delta    %Fano-genie3
% floor(gama(index(N))/3/delta)   %%Fano-max

l=2;%第二个译码器
[VPC(l),gama(l,:),V(l,:),T(l),P(l,:),b(l)] = individual_decoder(r,delta,maxcycles,poly,RP,sigma,E0,l, gama(l,:),V(l,:),T(l),P(l,:) );
phi=sum(VPC(l));%复杂度计算
VPC=VPC(l);
data=V(l,:);%译码结果

end

