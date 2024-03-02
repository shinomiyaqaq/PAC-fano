function [phi,VPC,data]= ...
    pac_decode_parallel(r,delta,maxcycles,poly,RP,sigma,E0)
N = size(r,2);
L=4;%��������������
% T=[0,-4,-6,-8];
T=[0,0,0,0];               %��ʼ��ֵ
a=zeros(1,L);
i=ones(1,L);
b=zeros(1,L);

gama=zeros(L,N+2);%·������
gama(:,1)= -inf;
gama(:,2)=0;
V=zeros(L,N);%Ԥ������
P=zeros(L,N);%�Ĵ����洢����
a(1)=1;%������һ��������
% T(1)=0;
VPC=zeros(1,L);%��¼����ʱ��

while (sum(a)~=0)%��ʼ����
    for l=1:L
        if a(l)==1
            [VPC(l),gama(l,:),V(l,:),i(l),T(l),P(l,:),b(l)] = individual_decoder(r,delta,maxcycles,poly,RP,sigma,E0,l, gama(l,:),V(l,:),i(l),T(l),P(l,:) );%��������
        end
        
        if b(l)==true && l<=L-1
            a(l+1)=1;
            %             T(l+1)=min(gama(l,3:N+2));
            T(l+1)= floor(min(gama(l,3:N+2))/delta) *delta;%������ֵ
            gama(l+1,:)=gama(l,:);%·������
            V(l+1,:)=V(l,:);%Ԥ������
            P(l+1,:)=P(l,:);
            i(l+1)=i(l);%����λ��
            b(l)=0;
        end
        
        if i(l)==N
            a(l:L)=0;
        end
        
%         if T(l)<T(l+1)%��ͣ����
%             a(l:l)=0;
%         end
        
        for x=1:L
            if i(x)==N+1
                
                                for j=1:L-1
                                    [VPC(j),gama(j,:),V(j,:),i(j),T(j),P(j,:),b(j)] = individual_decoder2(r,delta,maxcycles,poly,RP,sigma,E0,l, gama(j,:),V(j,:),i(j),T(j),P(j,:) );%��1��L-1��Ҳ����һ��
                                end
                
%                 [~,index] = sort(gama(:,N+2),'descend'); %����
                [~,index] = sort(VPC,'descend'); %����
                for m=1:L
                    g=index(m);
                    if gama(g,N+2)~=0
                        
                        data=V(g,:);
                        phi=sum(VPC);
                        VPC=VPC(g);
                        return%������
                    end
                end
            end
        end
        
    end
end
end

