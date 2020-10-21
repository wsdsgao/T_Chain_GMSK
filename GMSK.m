function [phi, I_sig_BB, Q_sig_BB, int_phi] = GMSK(bit_5, phi_last, g)
% input: 5��������              bit_5 
%        ��ǰ��Ƶ�ز�            fc 
%        ��һbitĩ��λ          phi_last(��bit_5�ĵ�����bitֵ)
%        g                    ��ӦBT=0.3��g����
% output����ǰbit��ĩ��λ        phi
%         GMSK����I·��������    I_sig_BB
%         GMSK����Q·��������    Q_sig_BB
%         ��Ӧ��ǰbit����λ·��   int_phi             

bit_rate = 16e6;  % ��������
T = 1/bit_rate;  % ����ʱ��
fs_IF = 64e6;  % ��Ƶ�źŲ���Ƶ��
oversamp_IF = T*fs_IF;  % ����������
T_s_IF = 1/fs_IF;  % ����ʱ��
BT = 0.3;  % ��˹�˲��� BT����

c_bit = bit_5;

% ��λ�ۼ�
intphi = zeros(1,T*fs_IF);

for i = 1:T*fs_IF
    if i == 1
        intphi(i) = phi_last + c_bit(5)*g(i)*T_s_IF + c_bit(4)*g(i + T*fs_IF)*T_s_IF + c_bit(3)*g(i + 2*T*fs_IF)*T_s_IF + c_bit(2)*g(i + 3*T*fs_IF)*T_s_IF + c_bit(1)*g(i + 4*T*fs_IF)*T_s_IF;
    else
        intphi(i) = intphi(i-1) + c_bit(5)*g(i)*T_s_IF + c_bit(4)*g(i + T*fs_IF)*T_s_IF + c_bit(3)*g(i + 2*T*fs_IF)*T_s_IF + c_bit(2)*g(i + 3*T*fs_IF)*T_s_IF + c_bit(1)*g(i + 4*T*fs_IF)*T_s_IF;
    end
end


t = -T/2:1/fs_IF:T/2-1/fs_IF;
t = t(1:end) + (T/oversamp_IF/2);  % ���ĶԳ�
int_phi = mod(pi/(2*T).*intphi, 2*pi);

phi = int_phi(T*fs_IF)/(pi/(2*T));

% I·�ź�
I_sig_BB = cos(int_phi);
% I_sig_IF = cos(2*pi*f_c.*t);


% Q·�ź�
Q_sig_BB = sin(int_phi);
% Q_sig_IF = sin(2*pi*f_c.*t);

