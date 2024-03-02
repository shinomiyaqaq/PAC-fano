function [metric,cycles,data,VPC] = ...
    pac_decode(r,delta,maxcycles,poly,RP,sigma,E0)

visits = zeros(1,length(RP));%��¼�ڵ�ķ��ʴ���
VPC = 0;
polyb = dec2bin(base2dec(num2str(poly), 8))-'0';
polyL = polyb == 1;
c = length(polyb) - 1;      %���쳤��


N = size(r,2);                          %���ֳ���
i = 1;                                  % ��ǰ�ڵ�
T = 0;                                  % ��ֵ
M = zeros(N+2,1);                       % ·��������
M(1) = -inf;                            % ���ڵ㷵��
M(2) = 0;                               % ���ڵ�·������

t = zeros(N,1);                         % �ڷ��صĹ������ж��˻ص��ĸ��ڵ�
inpseq = zeros(N,1);                    %��ǰ����
parseq = zeros(N,1);                    %��������ʱ���
tm = zeros(2,1);                        %��ǰ·�������ıȽ�
tb = zeros(2,1);                        %  ��ǰ·��
parity = zeros(2,1);                    % ��ʱ�洢��·������

reg = zeros(1,c + N);                   %��λ�Ĵ��� ����Ľ��

maxcycles = maxcycles * N;              %�������
dobreakA = 0;

% ·����������
n = log2(N);%����
llr = zeros(1,N);       % ����
L = zeros(n+1,N);       % �ڵ�洢�Ķ�����Ȼ��
ucap = zeros(n+1,N);    %  �о����
ns = zeros(1,2*N-1);    % �ڵ�״̬

f = @(a,b) (1-2*(a<0)).*(1-2*(b<0)).*min(abs(a),abs(b));%f����
g = @(a,b,c) b+(1-2*c).*a;%g����
L(1,:) = -2*(1/sigma^2)*r;      % ������Ȼ��
node = 0; depth = 0;    % �Ӹ��ڵ�㿪ʼ
% done = 0;             % ��������

cycles = 0;%����

while cycles <= maxcycles % ѭ��A
    t(i) = 1; %��ʼ��
    if(dobreakA)
        break;
    end
    while(1) %ѭ��B
        cycles = cycles + 1;
        if i - 1 < node  %��0��ʼ
            ucap(n+1,:) = parseq;
            for I = (node+1):-1:i+1 % Ҷ�ڵ�
                npos = (2^depth-1) + I;
                
                ns(npos) = 0;
                tempU = npos;
                while (mod(tempU,2) == 0)
                    tempU = tempU/2;
                    ns(tempU) = 0;
                end
                tempL = floor(tempU/2); %�����ƶ�
                ns(tempL) = 1;
            end
            depth = n; %������ȴﵽn
            node = (i-1);
            
            llr(node+1) = L(n+1,node+1); % LLR
            %%%%%%%%%%%%%%%%%%%%%%%%%End%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            stop = i;
            ucap(n+1,:) = parseq;
            done = 0;               %��������
            while(done == 0)
                if depth == n  % �ж��Ƿ���Ҷ�ڵ�
                    llr(node+1) = L(n+1,node+1);
                    if node == (stop-1)
                        done = 1;
                    else
                        node = floor(node/2); depth = depth - 1;
                    end
                else
                    % ��Ҷ�ڵ�
                    npos = (2^depth-1) + node + 1;  %�ڵ��λ��
                    if ns(npos) == 0 %ǰ����֧
                        temp = 2^(n-depth);
                        Ln = L(depth+1,temp*node+1:temp*(node+1));  %LLR
                        a = Ln(1:temp/2); b = Ln(temp/2+1:end);
                        node = node*2; depth = depth+1; %��һ���ڵ�
                        temp = temp/2;
                        L(depth+1,temp*node+1:temp*(node+1)) = f(a,b);%f���� ����LLR
                        ns(npos) = 1;
                    else
                        if ns(npos) == 1 %ǰ����֧
                            temp = 2^(n-depth);
                            Ln = L(depth+1,temp*node+1:temp*(node+1));  % LLR
                            a = Ln(1:temp/2); b = Ln(temp/2+1:end);
                            lnode = 2*node; ldepth = depth+1;
                            ltemp = temp/2;
                            ucapn = ucap(ldepth+1,ltemp*lnode+1:ltemp*(lnode+1));
                            node = node*2+1; depth = depth+1;
                            temp = temp/2;
                            L(depth+1,temp*node+1:temp*(node+1)) = g(a,b,ucapn);%g���� ����LLR
                            ns(npos) = 2;
                        else %ǰ�����ڵ�
                            temp = 2^(n-depth);
                            lnode = 2*node; rnode = 2*node+1; cdepth = depth+1; %���ҷ�֧
                            ctemp = temp/2;
                            ucapl = ucap(cdepth+1,ctemp*lnode+1:ctemp*(lnode+1));
                            ucapr = ucap(cdepth+1,ctemp*rnode+1:ctemp*(rnode+1));
                            ucap(depth+1,temp*node+1:temp*(node+1)) = [mod(ucapl+ucapr,2) ucapr];
                            node = floor(node/2); depth = depth - 1;
                            ns(npos) = 3;
                        end
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%End%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        ri = llr(i);
        state = reg(1:c);
        
        if ~RP(i) % �������
            visits(i) = visits(i) + 1;
            u0 = genparity_Rone_logical(0,state,polyL); %�������
            m = calc_met_gaussian(ri,u0,E0(i));
            M(i+2) = M(i+1) + m;
            inpseq(i) = 0;
            parseq(i) = u0;
        else      %  ��Ϣ����
            visits(i) = visits(i) + 1;
            u0 = genparity_Rone_logical(0,state,polyL);
            m0 = calc_met_gaussian(ri,u0,E0(i)); % ��֧����
            
            u1 = genparity_Rone_logical(1,state,polyL);
            m1 = calc_met_gaussian(ri,u1,E0(i));
            
            if m1 > m0
                tm(1) = m1;
                tb(1) = 1;
                parity(1) = u1;
                tm(2) = m0;
                tb(2) = 0;
                parity(2) = u0;
            else
                tm(1) = m0;
                tb(1) = 0;
                parity(1) = u0;
                tm(2) = m1;
                tb(2) = 1;
                parity(2) = u1;
            end
            
            M(i+2) = M(i+1) + tm(t(i));
            inpseq(i) = tb(t(i));
            parseq(i) = parity(t(i));
        end
        
        if M(i+2) >= T % �ʵ�������ֵ
            if M(i+1) < T + delta
                while T + delta <= M(i+2)
                    T = T + delta;
                end
            end
            reg = [inpseq(i) reg(1,1:end-1)];
            i = i + 1;
            VPC = VPC + 1;
            if i == N + 1 % �������
                dobreakA = 1;
            end
            break         %ǰ��ѭ��A
        else  %  ·������û������ֵ
            dobreak = 0;
            while(1) % ѭ��C
                if M(i) <  T
                    T = T - delta;
                    dobreak = 1;
                    break
                else
                    %  ���ݿ����µ������� �ж�
                    i = i - 1;
                    reg = [reg(2:end) 0];
                    t(i) = t(i) + 1;
                    if ~RP(i)    % �Ƿ�֧ ����
                        continue
                    end
                    if t(i) == 2
                    end
                    
                    if t(i) == 3 % ���Ը���Ľڵ� ����
                        continue;       % ת��ѭ��C
                    end
                    break               % ǰ��ѭ��B
                end
            end %ѭ��C����
            if(dobreak)
                break;                  %����ѭ��Bǰ��ѭ��A
            end
        end
        
    end % ѭ��B����
end
data = inpseq;%������
metric = M;%·������
end


