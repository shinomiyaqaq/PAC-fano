function [phi,VPC,data]= ...
    pac_decode_genie(r,delta,maxcycles,poly,RP,sigma,E0)
N = size(r,2);                          % ���ֳ���
L=2;%����������
% T=[0,-4,-6,-8];
T=[0,0];                                  % threshold ��ʼ��ֵ

gama=zeros(L,N+2); %·������
gama(:,1)= -inf;%��ʼ��
gama(:,2)=0;
V=zeros(L,N);%�ݴ���
P=zeros(L,N);%�ݴ����
VPC=zeros(1,L);%��¼����ʱ��

l=1;%��һ��������
[VPC(l),gama(l,:),V(l,:),T(l),P(l,:),b(l)] = individual_decoder(r,delta,maxcycles,poly,RP,sigma,E0,l, gama(l,:),V(l,:),T(l),P(l,:) );%�������
[~,index] = sort(gama(l,3:N+2));
T(l+1)= floor(gama(index(1))/3/delta) *3*delta;%ѡ���ڶ����������ĳ�ʼT

% clc
% floor(gama(index(1))/delta) *delta     %Fano-genie
% floor(gama(index(1))/2/delta) *2*delta  %Fano-genie2
% floor(gama(index(1))/3/delta) *3*delta    %Fano-genie3
% floor(gama(index(N))/3/delta)   %%Fano-max

l=2;%�ڶ���������
[VPC(l),gama(l,:),V(l,:),T(l),P(l,:),b(l)] = individual_decoder(r,delta,maxcycles,poly,RP,sigma,E0,l, gama(l,:),V(l,:),T(l),P(l,:) );
phi=sum(VPC(l));%���Ӷȼ���
VPC=VPC(l);
data=V(l,:);%������

end

