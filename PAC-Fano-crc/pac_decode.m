function [metric,cycles,data,VPC] = ...
    pac_decode(r,delta,maxcycles,poly,RP,sigma,E0)

visits = zeros(1,length(RP));%记录节点的访问次数
VPC = 0;
polyb = dec2bin(base2dec(num2str(poly), 8))-'0';
polyL = polyb == 1;
c = length(polyb) - 1;      %构造长度


N = size(r,2);                          %码字长度
i = 1;                                  % 当前节点
T = 0;                                  % 阈值
M = zeros(N+2,1);                       % 路径度量表
M(1) = -inf;                            % 根节点返回
M(2) = 0;                               % 根节点路径度量

t = zeros(N,1);                         % 在返回的过程中判断退回到哪个节点
inpseq = zeros(N,1);                    %当前输入
parseq = zeros(N,1);                    %卷积后的暂时结果
tm = zeros(2,1);                        %当前路径度量的比较
tb = zeros(2,1);                        %  当前路径
parity = zeros(2,1);                    % 暂时存储的路径度量

reg = zeros(1,c + N);                   %移位寄存器 卷积的结果

maxcycles = maxcycles * N;              %最大轮数
dobreakA = 0;

% 路径度量参数
n = log2(N);%阶数
llr = zeros(1,N);       % 输入
L = zeros(n+1,N);       % 节点存储的对数似然比
ucap = zeros(n+1,N);    %  判决结果
ns = zeros(1,2*N-1);    % 节点状态

f = @(a,b) (1-2*(a<0)).*(1-2*(b<0)).*min(abs(a),abs(b));%f函数
g = @(a,b,c) b+(1-2*c).*a;%g函数
L(1,:) = -2*(1/sigma^2)*r;      % 对数似然比
node = 0; depth = 0;    % 从根节点点开始
% done = 0;             % 结束译码

cycles = 0;%轮数

while cycles <= maxcycles % 循环A
    t(i) = 1; %初始化
    if(dobreakA)
        break;
    end
    while(1) %循环B
        cycles = cycles + 1;
        if i - 1 < node  %从0开始
            ucap(n+1,:) = parseq;
            for I = (node+1):-1:i+1 % 叶节点
                npos = (2^depth-1) + I;
                
                ns(npos) = 0;
                tempU = npos;
                while (mod(tempU,2) == 0)
                    tempU = tempU/2;
                    ns(tempU) = 0;
                end
                tempL = floor(tempU/2); %向右移动
                ns(tempL) = 1;
            end
            depth = n; %译码深度达到n
            node = (i-1);
            
            llr(node+1) = L(n+1,node+1); % LLR
            %%%%%%%%%%%%%%%%%%%%%%%%%End%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            stop = i;
            ucap(n+1,:) = parseq;
            done = 0;               %结束译码
            while(done == 0)
                if depth == n  % 判断是否是叶节点
                    llr(node+1) = L(n+1,node+1);
                    if node == (stop-1)
                        done = 1;
                    else
                        node = floor(node/2); depth = depth - 1;
                    end
                else
                    % 非叶节点
                    npos = (2^depth-1) + node + 1;  %节点的位置
                    if ns(npos) == 0 %前往左支
                        temp = 2^(n-depth);
                        Ln = L(depth+1,temp*node+1:temp*(node+1));  %LLR
                        a = Ln(1:temp/2); b = Ln(temp/2+1:end);
                        node = node*2; depth = depth+1; %下一个节点
                        temp = temp/2;
                        L(depth+1,temp*node+1:temp*(node+1)) = f(a,b);%f函数 更新LLR
                        ns(npos) = 1;
                    else
                        if ns(npos) == 1 %前往右支
                            temp = 2^(n-depth);
                            Ln = L(depth+1,temp*node+1:temp*(node+1));  % LLR
                            a = Ln(1:temp/2); b = Ln(temp/2+1:end);
                            lnode = 2*node; ldepth = depth+1;
                            ltemp = temp/2;
                            ucapn = ucap(ldepth+1,ltemp*lnode+1:ltemp*(lnode+1));
                            node = node*2+1; depth = depth+1;
                            temp = temp/2;
                            L(depth+1,temp*node+1:temp*(node+1)) = g(a,b,ucapn);%g函数 更新LLR
                            ns(npos) = 2;
                        else %前往父节点
                            temp = 2^(n-depth);
                            lnode = 2*node; rnode = 2*node+1; cdepth = depth+1; %左右分支
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
        
        if ~RP(i) % 冻结比特
            visits(i) = visits(i) + 1;
            u0 = genparity_Rone_logical(0,state,polyL); %卷积编码
            m = calc_met_gaussian(ri,u0,E0(i));
            M(i+2) = M(i+1) + m;
            inpseq(i) = 0;
            parseq(i) = u0;
        else      %  信息比特
            visits(i) = visits(i) + 1;
            u0 = genparity_Rone_logical(0,state,polyL);
            m0 = calc_met_gaussian(ri,u0,E0(i)); % 分支度量
            
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
        
        if M(i+2) >= T % 适当增加阈值
            if M(i+1) < T + delta
                while T + delta <= M(i+2)
                    T = T + delta;
                end
            end
            reg = [inpseq(i) reg(1,1:end-1)];
            i = i + 1;
            VPC = VPC + 1;
            if i == N + 1 % 译码完成
                dobreakA = 1;
            end
            break         %前往循环A
        else  %  路径度量没超过阈值
            dobreak = 0;
            while(1) % 循环C
                if M(i) <  T
                    T = T - delta;
                    dobreak = 1;
                    break
                else
                    %  回溯开启新的译码器 判断
                    i = i - 1;
                    reg = [reg(2:end) 0];
                    t(i) = t(i) + 1;
                    if ~RP(i)    % 非分支 回溯
                        continue
                    end
                    if t(i) == 2
                    end
                    
                    if t(i) == 3 % 来自更差的节点 回溯
                        continue;       % 转跳循环C
                    end
                    break               % 前往循环B
                end
            end %循环C结束
            if(dobreak)
                break;                  %跳出循环B前往循环A
            end
        end
        
    end % 循环B结束
end
data = inpseq;%译码结果
metric = M;%路径度量
end


