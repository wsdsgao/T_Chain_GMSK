clear all; 
close all;

rng(0);

bit_rate = 16e6;  % 符号速率
T = 1/bit_rate;  % 符号时间

f_IF = 20e6; %射频频率
fs_IF = 64e6;  % 射频、中频频信号采样速率
fs_BB = 128e6;  % 基带信号采样速率
num_bits_pulse = 300; 
oversamp_BB = T * fs_BB;  % 基带信号过采样速率
oversamp_IF = T * fs_IF;  % 射频、中频信号过采样速率
T_s_BB = 1/fs_BB;  % 基带采样间隔

load('lib/g_1024.mat');  % GMSK调制 g函数 

g = g(1:16:end);

% state = zeros(32,2);

% for k = 1:32
%     tmp = dec2bin(k-1, 5);
%     state(k,1) = bin2dec([tmp(2:end), '0']) + 1;
%     state(k,2) = bin2dec([tmp(2:end), '1']) + 1;
% end
error_rate = 0;

for ll = 1 : 1

% generate code

% I_single = randi(2,1,num_bits_pulse);
% I_single = I_single - 1;
load('test_I.mat');
I = 2*I_single - 1;
% I_single = [0,0,1,0,0,0,1,1,1,1];
% I = [1,-1,1,-1,-1,-1,1,1,1,1];

% coding
bit_5 = zeros(1,5);
phi_last = 0;

for i = 1:num_bits_pulse
    if i == 1
        bit_5 = [-1,-1,I(i:i+2)];
    elseif i == 2
        bit_5 = [-1,I(i-1:i+2)];
    elseif i == num_bits_pulse-1
        bit_5 = [I(i-2:i+1),-1];
    elseif i == num_bits_pulse
        bit_5 = [I(i-2:i),-1,-1];
    else
        bit_5 = I(i-2:i+2);
    end

    [phi_last, I_sig, Q_sig, phi_int] = GMSK(bit_5, phi_last, g);
    signal_trans_BB((i-1)*oversamp_IF+1:(i)*oversamp_IF) = complex(I_sig, Q_sig);
    phi_all((i-1)*oversamp_IF+1:(i)*oversamp_IF) = phi_int;
end

% figure;
% plot(mod(phi_all,2*pi))

% 射频信号

t = linspace(0, num_bits_pulse*T, oversamp_IF*num_bits_pulse);
signal_trans_IF = signal_trans_BB .* exp(1i*2*pi*0*t);

fre = t./(num_bits_pulse*T)*fs_IF;
% figure;
% plot(fre, 20*log10(abs(fft(signal_trans_IF))))

% 加频偏
% signal_trans_IF = signal_trans_IF .* exp(1i*2*pi*0.001*f_IF*t);


% 加噪声
Eb_N0 = 5;
SNRdB = Eb_N0 - 10*log10(oversamp_IF);
signal_recv_IF_noise = awgn(signal_trans_IF, SNRdB, 'measured');
% signal_recv_IF_noise = signal_trans_IF;
signal_recv_noise_IF_FFT = abs(fft(signal_recv_IF_noise));

% figure;
% plot(fre, 20*log10(signal_recv_noise_IF_FFT))

% 带通滤波

% 低通滤波

LPF = [-2.09996781279780e-06,1.48709839477551e-05,2.41084066936802e-05,7.67786559304659e-06,-4.09259577702220e-05,-8.00291950244231e-05,-4.28121775662192e-05,8.39217038174333e-05,0.000198996940063321,0.000143352961846203,-0.000126819855599412,-0.000410207352970788,-0.000369458463558583,0.000126671332563022,0.000733581280914653,0.000802565752592462,-6.20571789460999e-07,-0.00116213178237107,-0.00153703156923141,-0.000383445772287268,0.00164046974116608,0.00266364777074907,0.00121047725370158,-0.00204234618765731,-0.00424686114244369,-0.00271584085615062,0.00214956679334894,0.00629995763102076,0.00518195988408056,-0.00163029888540793,-0.00876460407608604,-0.00894949222251640,2.37850893950213e-06,0.0115012888139029,0.0144823392483473,0.00346706512785841,-0.0142957653923261,-0.0225989826221139,-0.0100919992380440,0.0168829707471250,0.0352960794540586,0.0230493817765340,-0.0189854455536604,-0.0595950632578754,-0.0546144621111168,0.0203591687539180,0.145950735187827,0.264171020676178,0.312496426755729,0.264171020676178,0.145950735187827,0.0203591687539180,-0.0546144621111168,-0.0595950632578754,-0.0189854455536604,0.0230493817765340,0.0352960794540586,0.0168829707471250,-0.0100919992380440,-0.0225989826221139,-0.0142957653923261,0.00346706512785841,0.0144823392483473,0.0115012888139029,2.37850893950213e-06,-0.00894949222251640,-0.00876460407608604,-0.00163029888540793,0.00518195988408056,0.00629995763102076,0.00214956679334894,-0.00271584085615062,-0.00424686114244369,-0.00204234618765731,0.00121047725370158,0.00266364777074907,0.00164046974116608,-0.000383445772287268,-0.00153703156923141,-0.00116213178237107,-6.20571789460999e-07,0.000802565752592462,0.000733581280914653,0.000126671332563022,-0.000369458463558583,-0.000410207352970788,-0.000126819855599412,0.000143352961846203,0.000198996940063321,8.39217038174333e-05,-4.28121775662192e-05,-8.00291950244231e-05,-4.09259577702220e-05,7.67786559304659e-06,2.41084066936802e-05,1.48709839477551e-05,-2.09996781279780e-06];
signal_recv_BB = filtfilt(LPF, 1, signal_recv_IF_noise);
% signal_recv_BB = signal_recv_IF_noise;
% figure;
% plot(fre, 20*log10(abs(fft(signal_recv_BB))))


% viterbi译码

de_out = zeros(size(I));
% signal_recv_BB = signal_trans_IF;
viterbi_deep = 10;
% signal_trans_IF(1 : oversamp_IF*viterbi_deep)
for n = 1:viterbi_deep:num_bits_pulse
    % de_out = GMSK_demod(signal_recv_BB(1+(i-1)*oversamp_IF*viterbi_deep:i*oversamp_IF*viterbi_deep), de_out, i, viterbi_deep, g);
    signal_recv = signal_recv_BB(1+(n-1)*oversamp_IF:(n - 1 + viterbi_deep)*oversamp_IF);
    init_path = zeros(8, 3);
    path_record = cell(32, viterbi_deep);
    % path_record_tmp = cell(32, 1);
    path_weight = zeros(32, viterbi_deep);
    
    init_path(1,:) = [0,0,0];
    init_path(2,:) = [0,0,1];
    init_path(3,:) = [0,1,0];
    init_path(4,:) = [0,1,1];
    init_path(5,:) = [1,0,0];
    init_path(6,:) = [1,0,1];
    init_path(7,:) = [1,1,0];
    init_path(8,:) = [1,1,1];

    if n == 1
        init_path = [repmat(zeros(1,2),8,1), init_path];
    else
        init_path = [repmat(de_out(n-2:n-1),8,1), init_path];
    end

    % 初始化放入path_record中

    phi_last = 0;
    for j = 1:8
        path_record{bin2dec(num2str(init_path(j,:)))+1, 1} = init_path(j,:);
        bit_5 = init_path(j,:) * 2 - 1;
        g_path_ang_dif = zeros(1,2);
        signal_recv_ang_dif = zeros(1,2);
        [~, I_sig, Q_sig, ~] = GMSK(bit_5, phi_last, g);
        signal_recv_ang = angle(signal_recv(1:4));
        g_path_ang = angle(complex(I_sig, Q_sig));
        
        for m = 1:2
            g_path_ang_dif(m) = g_path_ang(m) - g_path_ang(m+2);
            if  g_path_ang_dif(m) > pi
                g_path_ang_dif(m) =  g_path_ang_dif(m) - 2*pi;
            elseif  g_path_ang_dif(m) < -pi
                g_path_dif_dif(m) =  g_path_ang_dif(m) + 2*pi;
            end
            signal_recv_ang_dif(m) = signal_recv_ang(m) - signal_recv_ang(m+2);
            if signal_recv_ang_dif(m) > pi
                signal_recv_ang_dif(m) = signal_recv_ang_dif(m) - 2*pi;
            elseif signal_recv_ang_dif(m) < -pi
                signal_recv_ang_dif(m) = signal_recv_ang_dif(m) + 2*pi;
            end
        end

        path_weight(j,1) =  path_weight(j,1) + (sum(signal_recv_ang_dif) - sum(g_path_ang_dif)).^2;

    end

    % path_record_tmp = path_record;

    for i = 2:viterbi_deep
        signal_recv_ang = angle(signal_recv(1+(i-1)*4:i*4));
        for j = 1 : 32
            if ~isempty(path_record{j,i-1})
                bit_record = path_record{j,i-1};
                bit_5 = bit_record(end-3:end);
                
                bit_record_cur0 = [bit_record, 0];
                bit_record_cur1 = [bit_record, 1];
                [~, I_sig0, Q_sig0, ~] = GMSK([bit_5*2-1,-1], 0, g);
                [~, I_sig1, Q_sig1, ~] = GMSK([bit_5*2-1,1], 0, g);
                
                state_index0 = bin2dec(num2str([bit_record(end-3:end),0]))+1;
                state_index1 = bin2dec(num2str([bit_record(end-3:end),1]))+1;
                
                g_path_ang_dif = zeros(1,2);
                signal_recv_ang_dif = zeros(1,2);
                g_path_ang0 = angle(complex(I_sig0, Q_sig0));
                g_path_ang1 = angle(complex(I_sig1, Q_sig1));
                
                for m = 1:2
                    g_path_ang_dif0(m) = g_path_ang0(m) - g_path_ang0(m+2);
                    if  g_path_ang_dif0(m) > pi
                        g_path_ang_dif0(m) =  g_path_ang_dif0(m) - 2*pi;
                    elseif  g_path_ang_dif0(m) < -pi
                        g_path_dif_dif0(m) =  g_path_ang_dif0(m) + 2*pi;
                    end
                    g_path_ang_dif1(m) = g_path_ang1(m) - g_path_ang1(m+2);
                    if  g_path_ang_dif1(m) > pi
                        g_path_ang_dif1(m) =  g_path_ang_dif1(m) - 2*pi;
                    elseif  g_path_ang_dif1(m) < -pi
                        g_path_dif_dif1(m) =  g_path_ang_dif1(m) + 2*pi;
                    end
                    signal_recv_ang_dif(m) = signal_recv_ang(m) - signal_recv_ang(m+2);
                    if signal_recv_ang_dif(m) > pi
                        signal_recv_ang_dif(m) = signal_recv_ang_dif(m) - 2*pi;
                    elseif signal_recv_ang_dif(m) < -pi
                        signal_recv_ang_dif(m) = signal_recv_ang_dif(m) + 2*pi;
                    end
                end

                if path_weight(state_index0,i) == 0 || path_weight(j,i-1) + (sum(signal_recv_ang_dif) - sum(g_path_ang_dif0)).^2 < path_weight(state_index0,i)
                    path_weight(state_index0,i) =  path_weight(j,i-1) + (sum(signal_recv_ang_dif) - sum(g_path_ang_dif0)).^2;
                    path_record{state_index0, i} = bit_record_cur0;                   
                end
                
                if path_weight(state_index1,i) == 0 || path_weight(j,i-1) + (sum(signal_recv_ang_dif) - sum(g_path_ang_dif1)).^2 < path_weight(state_index1,i)
                    path_weight(state_index1,i) =  path_weight(j,i-1) + (sum(signal_recv_ang_dif) - sum(g_path_ang_dif1)).^2;
                    path_record{state_index1, i} = bit_record_cur1;                   
                end
            end
        end

        % path_record = path_record_tmp;
    end

    [~,index_min] = min(path_weight(:,end));
    out = path_record{index_min,end};
    de_out(n:n+viterbi_deep-1) = out(3:end-2);
end

error = I_single - de_out;

error(error~=0) = 1

error_rate = error_rate + sum(error);
end


error_rate/num_bits_pulse


