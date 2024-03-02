clear; clc; close all

s = rng(100);       %�������
delta = 2; %��
maxcycles = 1e8;%�����������

% EbNo_vec = 1:0.5:2; % �����
EbNo_vec=0:0.5:3;
% EbNo_vec=0;

K = 64; N = 128; R = K/N; Rc = R; SNR_Cons = 1;
% K = 4; N = 8; R = K/N; Rc = R; SNR_Cons = 1;

% rate profile
type = 4; % 1--> GA, 3--> Tse_RMpolar, 4 --> RM
RP = RM_Polar_Profile(N, K, SNR_Cons, type);
sum(RP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u=GA(0.7449,128);%��˹���취
[~,sortIndex] = sort(u(:));%sortIndex �Ӹߵ������� ����ֵԽ���ŵ�Խ�� ǰN��λ��Ϊ�������
change_index=sortIndex(1:10);%flip ��ת����εļ������� PAC-Fano-modified-CRC-flip ֮ǰ�ı���Ҫ����ȥ��
crclen=16;%24 CRCУ���볤��
poly = 3211;
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
    EbNo = 10^(EbNo_dB/10);
    sigma = 1/sqrt(2*Rc*EbNo);
    ZW = Z_polarization_fast(N,sigma);
    E0 = log2(2./(1+ZW));
    load('E0.mat');
    E0=E0';
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
    
            for x=1:10%��ת����׵�10������ Fano-crc 
                final_crc_result=crc_check_result(decoded,K-crclen);
                if(~final_crc_result)%CRCУ��Ľ��
                    decoded(change_index(x))=1-decoded(change_index(x));
                    %  if(x~=1)
                    %   msg_crc(change_index(x-1))=1-msg_crc(change_index(x-1));
                    %   end
                else
                    break
                end
            end
            
            Nblkerrs = Nblkerrs + any(decoded(RP)~=msg');%ͳ�ƴ�����
        end
        if Nblkerrs >= maxFE
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
    fprintf(2,'[%02d:%02d] SNR = %.1f, Frame = %d, FE = %d,FER=%f��Z= %d',elapsed_h,elapsed_m,EbNo_dB,i*ClustNum(EbNo_count),Nblkerrs,Nblkerrs/(i*ClustNum(EbNo_count)),Visit);
    
    temp = (i*ClustNum(EbNo_count));
    Nruns(EbNo_count) = temp;
    FER(EbNo_count) = Nblkerrs/temp;
    FE(EbNo_count) = Nblkerrs;
    V(EbNo_count) = Visit;
    fprintf('-------------------------------------\n');
end

rng(s);    

%������
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
