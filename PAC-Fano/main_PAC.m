clear; clc; close all
s = rng(100);       %�������
delta = 2; %��
maxcycles = 1e8;%�����������
EbNo_vec = 0:0.5:3;%�����
K = 64; N = 128; R = K/N; Rc = R; SNR_Cons = 1;
% K = 4; N = 8; R = K/N; Rc = R; SNR_Cons = 1;%��Ϣ���س���  �볤 ����

% rate profile
type = 4; % 1--> GA, 3--> Tse_RMpolar, 4 --> RM
RP = RM_Polar_Profile(N, K, SNR_Cons, type);%���췽ʽ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

poly = 3211;%�����
polyb = dec2bin(base2dec(num2str(poly), 8)) - '0'; % 1 011 011
constraint = length(polyb); % 7
C = constraint - 1;
tr = poly2trellis(constraint,poly); % #of states is 2^C
outputs = tr.outputs;
nextStates = tr.nextStates;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FER = zeros(1,length(EbNo_vec)); % ��֡��
FE = zeros(1,length(EbNo_vec));  %�����֡�ĸ���
V = zeros(1,length(EbNo_vec));  % #ÿ��SNR�ڵ�ķ��ʴ���
Nruns = zeros(1,length(EbNo_vec)); %ʵ����������

maxRun = 1e7;%������д���

ClustNum = [500 500 1e3 5e3*(ones(1,length(EbNo_vec)-3))]; % ���в���
maxFE = 80; % Ҫ��ﵽ�����֡�ĸ���

fprintf('-------------------------------------\n');
for EbNo_count = 1:length(EbNo_vec)
    tic;
    EbNo_dB = EbNo_vec(EbNo_count);
    EbNo = 10^(EbNo_dB/10);%�����ת��
    sigma = 1/sqrt(2*Rc*EbNo);%sigma
    
    ZW = Z_polarization_fast(N,sigma);%���ټ���
    E0 = log2(2./(1+ZW));
    load('E0.mat');
    E0=E0';%����fano·��������b
    Nblkerrs = 0;
    Visit = 0;
    
    fprintf('[%02d:%02d] Starting! SNR = %.2f\n',0,0, EbNo_dB);
    for i = 1:maxRun/ClustNum(EbNo_count)    
        for j = 1:ClustNum(EbNo_count)%����
            msg = randi([0 1],1,K);%��Ϣ����
            code = zeros(1,N);
            code(RP) = msg; %�����Ϣ����
            %���
            u = convenc(code,tr);
            %����
            x = polarencode(u);
            % BPSK����
            modulated = 2*x-1 ;
            % AWGN����˹�ŵ�
            r = modulated + randn(1,N)*sigma;
            % ����
            [~,~,decoded,VPC] = ...
                pac_decode(r,delta,maxcycles,poly,RP,sigma,E0);
            Visit = Visit + VPC;
            Nblkerrs = Nblkerrs + any(decoded(RP)~=msg');%����������
        end
        if Nblkerrs >= maxFE%�ﵽͳ������
            break;
        end
        t = toc;
        elapsed_m = t/60;
        elapsed_h = floor(elapsed_m/60);
        elapsed_m = floor(elapsed_m - elapsed_h*60);
    end
    t = toc;
    
    elapsed_m = t/60;
    elapsed_h = floor(elapsed_m/60);
    elapsed_m = floor(elapsed_m - elapsed_h*60);
    fprintf(2,'[%02f:%02f] SNR = %.1f, Frame = %d, FE = %d,FER=%f',elapsed_h,elapsed_m,EbNo_dB,i*ClustNum(EbNo_count),Nblkerrs,Nblkerrs/(i*ClustNum(EbNo_count)));
    
    temp = (i*ClustNum(EbNo_count));%ͳ��
    Nruns(EbNo_count) = temp;
    FER(EbNo_count) = Nblkerrs/temp;
    FE(EbNo_count) = Nblkerrs;
    V(EbNo_count) = Visit;
    fprintf('-------------------------------------\n');
    t/temp
    
end
rng(s);    
% ������
res.trials = Nruns;
res.frame_errors = FE;%����֡��
res.FER = FE./Nruns;%��֡��
res.avg_visit = V./Nruns;%���ʴ���

res.SNR = EbNo_vec;%�����
res.rate_profile = type;%���췽ʽ
res.K = K;%����
res.delta = delta;%��

res.max_cycles_allowed = maxcycles;
res.max_trials = maxRun;
res.max_frame_errors = maxFE;

filename = sprintf('outputs/pac_K%d_N%d_SNR%0.1f_delta%0.1f_poly%d.mat',K,N,EbNo_vec(1),delta,poly);
save(filename,'res');%�����ļ�

