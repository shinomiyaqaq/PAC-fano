function RP = RM_Polar_Profile(N, K, SNR, type)
R = K/N;

EbNo_X = 10^(SNR/10); %信噪比转换
sigma_X = 1/sqrt(2*R*EbNo_X);%sigma
ZW_X = Z_polarization_fast(N,sigma_X);%高斯构造

if type == 1% type 1 -------> 高斯构造
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PRP rate profile (K most reliable channels)%K个最信任的信道
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
    idx_w = find(w==i); % i 计算行权重
           
    E0_X = log2(2./(1+ZW_X));
    
    X = zeros(1,N);
    X(idx_w) = E0_X(idx_w);
    [~,I] = sort(X,'descend');  
    I(K-Klow1+1:end) = []; %最可靠的行
    RP(I) = true;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif type == 3% type 3 ------->Tse RMpolar构造
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %选取K个最信任的位置
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
elseif type == 4% type 4 ------->RM构造
    [RP, ~] = RM_profile(K,N);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, w] = RM_profile(K,N)
w = sum(dec2bin(0:N-1)-'0',2);%对行的权重进行计算
[~, index] = sort(w,'descend');
kindex = index(1:K);
P = false(1,N);
P(kindex) = true;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









