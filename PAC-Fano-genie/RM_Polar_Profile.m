function RP = RM_Polar_Profile(N, K, SNR, type)
R = K/N;

EbNo_X = 10^(SNR/10); %�����ת��
sigma_X = 1/sqrt(2*R*EbNo_X);%sigma
ZW_X = Z_polarization_fast(N,sigma_X);%��˹����

if type == 1% type 1 -------> ��˹����
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PRP rate profile (K most reliable channels)%K�������ε��ŵ�
    RP = false(1,N);    
    [~,I] = mink(ZW_X, K);
    RP(I) = true;
%     sum(RP)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif type == 2% type 2 -------> RM-polar
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    Klow1 = 0;
    n = log2(N);
    i = n;
    while Klow1 < K    
        if Klow1 + nchoosek(n,i) < K
            Klow1 = Klow1 + nchoosek(n,i);
        else
            break;
        end
        i = i-1;
    end
    
    [RP, w] = RM_profile(Klow1,N);
    idx_w = find(w==i); % i ������Ȩ��
           
    E0_X = log2(2./(1+ZW_X));
    
    X = zeros(1,N);
    X(idx_w) = E0_X(idx_w);
    [~,I] = sort(X,'descend');  
    I(K-Klow1+1:end) = []; %��ɿ�����
    RP(I) = true;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif type == 3% type 3 ------->Tse RMpolar����
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ѡȡK�������ε�λ��
    Klow2 = 0;
    n = log2(N);
    i = n;
    while Klow2 <= K    
        Klow2 = Klow2 + nchoosek(n,i);
        i = i-1;
    end
    [RP, ~] = RM_profile(Klow2,N);
        
    E0_X = log2(2./(1+ZW_X));    
    
    E = 1e6*E0_X;
    RP_X = RP.*E;
    [~,I] = maxk(RP_X, K);
    
    RP = false(1,N);
    RP(I) = true;
%     sum(RP)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif type == 4% type 4 ------->RM����
    [RP, ~] = RM_profile(K,N);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, w] = RM_profile(K,N)
w = sum(dec2bin(0:N-1)-'0',2);%���е�Ȩ�ؽ��м���
[~, index] = sort(w,'descend');
kindex = index(1:K);
P = false(1,N);
P(kindex) = true;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









