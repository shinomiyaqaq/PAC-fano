function ZW = Z_polarization_fast(N,sigma)%��˹���Ʒ� ���ټ���
% SNR_dB = SNR_vec(i); % 2.5dB
% ebn0 = 10^(SNR_dB/10); % R =0.5;
% sigma = 1/sqrt(2*R*ebn0); % sigma = 0.7499;

fc = @(t) phi_inv(1-(1-phi_fun(t))^2); %fc����
fv = @(t) 2*t; %fv����

n = log2(N);
m = zeros(1,2*N-1); %�������Ľڵ�
m(1) = 2/sigma^2; %��ʼ�����Ҳ�

for d = 1:n
    u = 0;
    for i = 0:2^d-1
        if u == 2^d
            break; % �ﵽ��ײ�
        end
        if mod(i+1,2) == 1 % ��֧����
            m(2^d+i) = fc(m(2^(d-1)+floor(i/2))) ;
        else %��֧����
            m(2^d+i) = fv(m(2^(d-1)+floor(i/2)));
        end
        u = u+1;
    end
end

ZW = exp(-m(N:2*N-1)/4); % Z(W) = e^(-1/2sigma^2), m = 2/sigma^2;

end

function phi = phi_fun(t)%phi����
load('phiF_vec.mat');%��ǰ����
if t == 0
    phi = 1;
elseif t < 50
    [~,I] = min(abs(A-t));
    phi = phi_vec(I) ;
    %     syms z
    %     F = tanh(z/2).*exp(-(z-t).^2./(4*t));
    %     phi = 1 - 1./(sqrt(4*pi*t)).*double(int(F,z,-100,100));
else
    phi = 1.1795e-06;
end

end

function t = phi_inv(y)%phi�溯��
load('phi_vector2.mat'); %���� x = .01:.01:45 ��ص�phiֵ
if y == 1
    t = 0;
else
    [~,I] = min(abs(phi_vector-y));
    t = x(I) ;
end
end
