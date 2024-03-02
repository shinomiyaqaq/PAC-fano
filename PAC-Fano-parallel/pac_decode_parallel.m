function [phi,VPC,data]= ...
    pac_decode_parallel(r,delta,maxcycles,poly,RP,sigma,E0)
N = size(r,2);
L=4;%独立译码器个数
% T=[0,-4,-6,-8];
T=[0,0,0,0];               %初始阈值
a=zeros(1,L);
i=ones(1,L);
b=zeros(1,L);

gama=zeros(L,N+2);%路径度量
gama(:,1)= -inf;
gama(:,2)=0;
V=zeros(L,N);%预估比特
P=zeros(L,N);%寄存器存储比特
a(1)=1;%开启第一个译码器
% T(1)=0;
VPC=zeros(1,L);%记录译码时间

while (sum(a)~=0)%开始译码
    for l=1:L
        if a(l)==1
            [VPC(l),gama(l,:),V(l,:),i(l),T(l),P(l,:),b(l)] = individual_decoder(r,delta,maxcycles,poly,RP,sigma,E0,l, gama(l,:),V(l,:),i(l),T(l),P(l,:) );%并行译码
        end
        
        if b(l)==true && l<=L-1
            a(l+1)=1;
            %             T(l+1)=min(gama(l,3:N+2));
            T(l+1)= floor(min(gama(l,3:N+2))/delta) *delta;%更新阈值
            gama(l+1,:)=gama(l,:);%路径度量
            V(l+1,:)=V(l,:);%预估比特
            P(l+1,:)=P(l,:);
            i(l+1)=i(l);%比特位置
            b(l)=0;
        end
        
        if i(l)==N
            a(l:L)=0;
        end
        
%         if T(l)<T(l+1)%早停策略
%             a(l:l)=0;
%         end
        
        for x=1:L
            if i(x)==N+1
                
                                for j=1:L-1
                                    [VPC(j),gama(j,:),V(j,:),i(j),T(j),P(j,:),b(j)] = individual_decoder2(r,delta,maxcycles,poly,RP,sigma,E0,l, gama(j,:),V(j,:),i(j),T(j),P(j,:) );%第1到L-1个也译码一下
                                end
                
%                 [~,index] = sort(gama(:,N+2),'descend'); %排序
                [~,index] = sort(VPC,'descend'); %排序
                for m=1:L
                    g=index(m);
                    if gama(g,N+2)~=0
                        
                        data=V(g,:);
                        phi=sum(VPC);
                        VPC=VPC(g);
                        return%输出结果
                    end
                end
            end
        end
        
    end
end
end

