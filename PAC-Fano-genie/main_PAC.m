clear; clc; close all%①可以加个crc 误帧率降低 早停策略 码率降低②△动态变化 ③不同polar构造方式比较和选取

%L个fano译码器 存储n*L个比特 异或单元个数 理论最长译码时间

s = rng(100);       %随机种子
delta = 2; %△
maxcycles = 1e8;%最大译码轮数

% EbNo_vec = 1:0.5:2; 
EbNo_vec=0:0.5:3;%信噪比
% EbNo_vec = 0;

K = 64; N = 128; R = K/N; Rc = R; SNR_Cons = 1;
% K = 4; N = 8; R = K/N; Rc = R; SNR_Cons = 1;

% rate profile
type = 4; % 1--> GA, 3--> Tse_RMpolar, 4 --> RM
RP = RM_Polar_Profile(N, K, SNR_Cons, type);
sum(RP);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

poly = 3211;
polyb = dec2bin(base2dec(num2str(poly), 8)) - '0'; % 1 011 011
constraint = length(polyb); % 7
C = constraint - 1;
tr = poly2trellis(constraint,poly); 
outputs = tr.outputs;
nextStates = tr.nextStates;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FER = zeros(1,length(EbNo_vec)); %误帧率
FE = zeros(1,length(EbNo_vec));  % #错误帧数
V = zeros(1,length(EbNo_vec));  % #每个SNR的访问次数
Nruns = zeros(1,length(EbNo_vec)); % #实际运行轮数
P = zeros(1,length(EbNo_vec));%暂存比特

maxRun = 1e7;%最大运行轮数

ClustNum = [500 500 1e3 5e3*(ones(1,length(EbNo_vec)-3))]; % 根据N和k进行调节
maxFE = 70; % 要求达到的误帧数

fprintf('-------------------------------------\n');
for EbNo_count = 1:length(EbNo_vec)
    tic;
    EbNo_dB = EbNo_vec(EbNo_count);
    EbNo = 10^(EbNo_dB/10);
    sigma = 1/sqrt(2*Rc*EbNo);
    
    ZW = Z_polarization_fast(N,sigma);
    E0 = log2(2./(1+ZW));
    
    Nblkerrs = 0;
    Visit = 0;
    Phi=0;
    
    fprintf('[%02d:%02d] Starting! SNR = %.2f\n',0,0, EbNo_dB);
    for i = 1:maxRun/ClustNum(EbNo_count)  
       for j = 1:ClustNum(EbNo_count)%仿真
            msg = randi([0 1],1,K);%信息比特
            code = zeros(1,N);
            code(RP) = msg; %添加信息比特
            %卷积
            u = convenc(code,tr);
            %极化
            x = polarencode(u);
            % BPSK调制
            modulated = 2*x-1 ;
            % AWGN过高斯信道
            r = modulated + randn(1,N)*sigma;
            % 解码
            [phi,VPC,decoded] = ...
                pac_decode_genie(r,delta,maxcycles,poly,RP,sigma,E0);
            Visit = Visit + VPC;
            Phi=Phi+phi;
            Nblkerrs = Nblkerrs + any(decoded(RP)~=msg);
        end
        if Nblkerrs >= maxFE
            break;
        end
        t = toc;
        elapsed_m = t/60;
        elapsed_h = floor(elapsed_m/60);
        elapsed_m = floor(elapsed_m - elapsed_h*60);
        fprintf(2,'[%02d:%02d] EbNo = %.2f, Frame = %d, FE = %d\n',elapsed_h,elapsed_m,EbNo_dB,i*ClustNum(EbNo_count),Nblkerrs);
    end
    t = toc;
    elapsed_m = t/60;
    elapsed_h = floor(elapsed_m/60);
    elapsed_m = floor(elapsed_m - elapsed_h*60);
    fprintf(2,'[%02d:%02d] SNR = %.1f, Frame = %d, FE = %d,FER=%f\n',elapsed_h,elapsed_m,EbNo_dB,i*ClustNum(EbNo_count),Nblkerrs,Nblkerrs/(i*ClustNum(EbNo_count)));
    
    temp = (i*ClustNum(EbNo_count));
    Nruns(EbNo_count) = temp;
    FER(EbNo_count) = Nblkerrs/temp;
    FE(EbNo_count) = Nblkerrs;
    V(EbNo_count) = Visit;
    P(EbNo_count) = Phi;
    fprintf('-------------------------------------\n');
end

rng(s);   

%  保存结果
res.trials = Nruns;
res.frame_errors = FE;%错误帧数
res.FER = FE./Nruns;%误帧率
res.time = V./Nruns;%访问次数

res.SNR = EbNo_vec;%信噪比
res.rate_profile = type;%构造方式
res.K = K;%码率
res.delta = delta;%△

res.max_cycles_allowed = maxcycles;
res.max_trials = maxRun;
res.max_frame_errors = maxFE;

filename = sprintf('outputs/pac_K%d_N%d_SNR%0.1f_delta%0.1f_poly%d.mat',K,N,EbNo_vec(1),delta,poly);
save(filename,'res');%保存文件