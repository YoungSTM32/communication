clc;
clear;
%原理参考了：https://zhuanlan.zhihu.com/p/57967971


%% 参数设置

N_sc=52;      %系统子载波数（不包括直流载波）、number of subcarrierA
N_fft=64;            % FFT 长度
N_cp=16;             % 循环前缀长度、Cyclic prefix
N_symbo=N_fft+N_cp;        % 1个完整OFDM符号长度
N_c=53;             % 包含直流载波的总的子载波数、number of carriers
M=4;               %4PSK调制
SNR=0:1:20;         %仿真信噪比
N_frm=10;            % 每种信噪比下的仿真帧数、frame
Nd=6;               % 每帧包含的OFDM符号数
P_f_inter=8;      %导频间隔 一共8个导频
data_station=[];    %导频位置
L=7;                %卷积码约束长度
tblen=6*L;          %Viterbi译码器回溯深度
stage = 3;          % m序列的阶数
ptap1 = [1 3];      % m序列的寄存器连接方式
regi1 = [1 1 1];    % m序列的寄存器初始值
alpha=0.3           %

%% 基带数据数据产生
P_data=randi([0 1],1,N_sc*Nd*N_frm);
%数据量=仿真帧数×每帧含OFDM符号数×子载波数

%% 信道编码（卷积码、或交织器）
%卷积码：前向纠错非线性码
%交织：使突发错误最大限度的分散化
trellis = poly2trellis(7,[133 171]);       %(2,1,7)卷积编码
code_data=convenc(P_data,trellis);


%% qpsk调制
data_temp1= reshape(code_data,log2(M),[])';             %以每组2比特进行分组，M=4
data_temp2= bi2de(data_temp1);                             %二进制转化为十进制
modu_data=pskmod(data_temp2,M,pi/M);              % 4PSK调制
scatterplot(modu_data);
title('QPSK调制星座图')             %星座图(也可以取实部用plot函数)

%% 扩频
%————————————————————————————————————————————————————————%
%扩频通信信号所占有的频带宽度远大于所传信息必需的最小带宽
%根据香农定理，扩频通信就是用宽带传输技术来换取信噪比上的好处，这就是扩频通信的基本思想和理论依据。
%扩频就是将一系列正交的码字与基带调制信号内积
%扩频后数字频率变成了原来的m倍。码片数量 = 2（符号数）* m（扩频系数）
%————————————————————————————————————————————————————————%

code = mseq(stage,ptap1,regi1,N_sc);     % 扩频码的生成
code = code * 2 - 1;         %将1、0变换为1、-1
modu_data=reshape(modu_data,N_sc,length(modu_data)/N_sc);
spread_data = spread(modu_data,code);        % 扩频
spread_data=reshape(spread_data,[],1);

%% 插入导频
P_f=3+3*1i;                       %Pilot frequency
P_f_station=[1:P_f_inter:N_fft];%导频位置
pilot_num=length(P_f_station);%导频数量

for img=1:N_fft                        %数据位置
    if mod(img,P_f_inter)~=1          %mod(a,b)就是求的是a除以b的余数
        data_station=[data_station,img];
    end
end
data_row=length(data_station);
data_col=ceil(length(spread_data)/data_row);

pilot_seq=ones(pilot_num,data_col)*P_f;%将导频放入矩阵
data=zeros(N_fft,data_col);%预设整个矩阵
data(P_f_station(1:end),:)=pilot_seq;%对pilot_seq按行取

if data_row*data_col>length(spread_data)
    data2=[spread_data;zeros(data_row*data_col-length(spread_data),1)];%将数据矩阵补齐，补0是虚载频~
else data2=spread_data;
end;
    

%% 串并转换
data_seq=reshape(data2,data_row,data_col);
data(data_station(1:end),:)=data_seq;%将导频与数据合并

%% IFFT
ifft_data=ifft(data);
figure(2)
subplot(3,1,1);
plot(abs(ifft_data(:,1)),'b');
title('原始单个OFDM符号');
xlabel('Time');
ylabel('Amplitude');

%% 插入保护间隔、循环前缀
Tx_cd=[ifft_data(N_fft-N_cp+1:end,:);ifft_data];%把ifft的末尾N_cp个数补充到最前面
subplot(3,1,2);
plot(abs(Tx_cd(:,1)));
xlabel('Time');
ylabel('Amplitude');
title('加CP的单个OFDM符号');

%% 加窗
Tx_window=zeros(size(Tx_cd));

% 通过矩阵点乘
Tx_window=Tx_cd.*repmat(rcoswindow(alpha,size(Tx_cd,1)),1,390);%390是矩阵Tx_cd的维度(80*390)
subplot(3,1,3);
plot(abs(Tx_window(:,1)))
title('加窗后的单个OFDM符号')
xlabel('Time');
ylabel('Amplitude');

%% 信号频谱
%未加窗信号频谱
%figure % 归一化
figure(3);
orgin_aver_power = 20*log10(mean(abs(fft(Tx_cd'))));
subplot(2,1,1)
plot((1:length(orgin_aver_power))/length(orgin_aver_power),orgin_aver_power)
hold on
plot(0:1/length(orgin_aver_power):1 ,-35,'rd')
hold on
axis([0 1 -40 max(orgin_aver_power)])
grid on
title('未加窗信号频谱')

%加窗信号频谱
orgin_aver_power = 20*log10(mean(abs(fft(Tx_window'))));
subplot(2,1,2)
plot((1:length(orgin_aver_power))/length(orgin_aver_power),orgin_aver_power)
hold on
plot(0:1/length(orgin_aver_power):1 ,-35,'rd')
hold on
axis([0 1 -40 max(orgin_aver_power)])
grid on
title('加窗信号频谱')

%% 并串转换
Tx_data=reshape(Tx_window,[],1);%由于传输需要

%% 信道（信道均衡1、信道均衡2）
 %预设好误码率的变量矩阵
 Ber_single=zeros(1,length(SNR)); %单径误码率矩阵设置 信道均衡1
 Ber2_single=zeros(1,length(SNR));
 
 Ber_mult=zeros(1,length(SNR));%多径误码率矩阵设置 信道均衡1
 Ber2_mult=zeros(1,length(SNR));
 
 Ber_spline=zeros(1,length(SNR));%单径误码率矩阵设置 信道均衡2
 Ber2_spline=zeros(1,length(SNR));
 
 Ber_offset=zeros(1,length(SNR));  %多径误码率矩阵设置 信道均衡2 有频率补偿
 Ber2_offset=zeros(1,length(SNR));

Tx_data = reshape(Tx_window,1,[]); % 变成时域一个完整信号，待传输
%signal_origin = reshape(Tx_window,1,[]); % 未加窗完整信号
mult_path_am = [1 0.2 0.1]; %  多径幅度
mult_path_time1=randi([1,50]); % 多径随机时延1-50
mult_path_time2=randi([1,50]);
path2 = 0.2*[zeros(1,mult_path_time1) Tx_data(1:end-mult_path_time1) ];
path3 = 0.1*[zeros(1,mult_path_time2) Tx_data(1:end-mult_path_time2) ];
Tx_mult = Tx_data + path2 + path3; % 多径信号

figure('menubar','none');
subplot(2,1,1);
plot(real(Tx_mult) )
title('多径下OFDM信号')
xlabel('Time/samples')
ylabel('Amplitude')
subplot(2,1,2)
plot(real(Tx_data) )
title('单径下OFDM信号')
xlabel('Time/samples')
ylabel('Amplitude')

%% 单径/多径信道后处理

for jj=1:length(SNR)
    channe_single=awgn(Tx_data,SNR(jj),'measured');
    channe_mult=awgn(Tx_mult ,SNR(jj),'measured');

    %添加高斯白噪声
    %y = awgn(x,SNR,SIGPOWER) 
    %如果SIGPOWER为'measured'，则函数将在加入噪声之前测定信号强度。
%% 串并转换
    Rx_data1_single=reshape(channe_single,N_fft+N_cp,[]);
    Rx_data1_mult=reshape(channe_mult,N_fft+N_cp,[]);
    
%% 去掉保护间隔、循环前缀
    Rx_data2_single=Rx_data1_single(N_cp+1:end,:);
    Rx_data2_mult=Rx_data1_mult(N_cp+1:end,:);

%% FFT
    fft_data_single=fft(Rx_data2_single);
    fft_data_mult=fft(Rx_data2_mult);
    
%% 信道估计与插值线性（信道均衡1）
    %单径信道
    data3_single=fft_data_single(1:N_fft,:); 
    Rx_pilot_single=data3_single(P_f_station(1:end),:); %接收到的导频
    h_single=Rx_pilot_single./pilot_seq; 
    H_single=interp1( P_f_station(1:end)',h_single,data_station(1:end)','linear','extrap');
    %分段线性插值：插值点处函数值由连接其最邻近的两侧点的线性函数预测。
    %对超出已知点集的插值点用指定插值方法(线性)计算函数值
   
    %多径信道
    data3_mult=fft_data_mult(1:N_fft,:); 
    Rx_pilota_mult=data3_mult(P_f_station(1:end),:); %接收到的导频
    h_mult=Rx_pilota_mult./pilot_seq; 
    H_mult=interp1( P_f_station(1:end)',h_mult,data_station(1:end)','linear','extrap');
    
   %% 多途径信道样条插值估计（信道均衡2）
     H_mult_spline=interp1( P_f_station(1:end)',h_mult,data_station(1:end)','spline','extrap');
     
%% 信道校正
    data_aftereq_single=data3_single(data_station(1:end),:)./H_single;
    data_aftereq_mult=data3_mult(data_station(1:end),:)./H_mult;
    
    data_aftereq_spline=data3_mult(data_station(1:end),:)./H_mult_spline;
    
%% 并串转换
    %单径信道(信道均衡方法1)
    data_aftereq_single=reshape(data_aftereq_single,[],1);
    data_aftereq_single=data_aftereq_single(1:length(spread_data));
    data_aftereq_single=reshape(data_aftereq_single,N_sc,length(data_aftereq_single)/N_sc);
    %多径信道(信道均衡方法1)
    data_aftereq_mult=reshape(data_aftereq_mult,[],1);
    data_aftereq_mult=data_aftereq_mult(1:length(spread_data));
    data_aftereq_mult=reshape(data_aftereq_mult,N_sc,length(data_aftereq_mult)/N_sc);
    %多径信道(信道均衡方法2)
    data_aftereq_spline=reshape(data_aftereq_spline,[],1);
    data_aftereq_spline=data_aftereq_spline(1:length(spread_data));
    data_aftereq_spline=reshape(data_aftereq_spline,N_sc,length(data_aftereq_spline)/N_sc);    
    
%% 解扩
    demspread_data_single = despread(data_aftereq_single,code);       % 数据解扩 
    demspread_data_mult = despread(data_aftereq_mult,code); 
   
    demspread_data_spline = despread(data_aftereq_spline,code); 

%% QPSK解调
    demodulation_data_single=pskdemod(demspread_data_single,M,pi/M);    
    De_data1_single = reshape(demodulation_data_single,[],1);
    De_data2_single = de2bi(De_data1_single);
    De_Bit_single = reshape(De_data2_single',1,[]);
    
    demodulation_data_mult=pskdemod(demspread_data_mult,M,pi/M);    
    De_data1_mult = reshape(demodulation_data_mult,[],1);
    De_data2_mult = de2bi(De_data1_mult);
    De_Bit_mult = reshape(De_data2_mult',1,[]);
    
    demodulation_data_spline=pskdemod(demspread_data_spline,M,pi/M);    
    De_data1_spline= reshape(demodulation_data_spline,[],1);
    De_data2_spline= de2bi(De_data1_spline);
    De_Bit_spline = reshape(De_data2_spline',1,[]);
%% （解交织）
%% 信道译码（维特比译码）
    trellis_single = poly2trellis(7,[133 171]);
    rx_c_de_single = vitdec(De_Bit_single,trellis_single,tblen,'trunc','hard');   %硬判决
    
    trellis_mult = poly2trellis(7,[133 171]);
    rx_c_de_mult = vitdec(De_Bit_mult,trellis_mult,tblen,'trunc','hard');   %硬判决
    
    trellis_spline = poly2trellis(7,[133 171]);
    rx_c_de_spline= vitdec(De_Bit_spline,trellis_spline,tblen,'trunc','hard');   %硬判决
%% 计算误码率
    [err,Ber2_single(jj)] = biterr(De_Bit_single(1:length(code_data)),code_data);%译码前的误码率
    [err, Ber_single(jj)] = biterr(rx_c_de_single(1:length(P_data)),P_data);%译码后的误码率
    
    [err,Ber2_mult(jj)] = biterr(De_Bit_mult(1:length(code_data)),code_data);%译码前的误码率
    [err, Ber_mult(jj)] = biterr(rx_c_de_mult(1:length(P_data)),P_data);%译码后的误码率
    
    [err,Ber2_spline(jj)] = biterr(De_Bit_spline(1:length(code_data)),code_data);%译码前的误码率
    [err, Ber_spline(jj)] = biterr(rx_c_de_spline(1:length(P_data)),P_data);%译码后的误码率   
    
    
end

%% 原始信号不通过信道均衡
%信道编码、QPSK调制后的信号为modu_data
%IFFT
ifft_data_original=ifft(modu_data);
%插入保护间隔，循环前缀
Tx_cd_original=[ifft_data_original(N_fft-N_cp+1:end,:);ifft_data_original];
%加窗 （通过矩阵点乘）
Tx_window_original=zeros(size(Tx_cd_original));
Tx_window_original=Tx_cd_original.*repmat(rcoswindow(alpha,size(Tx_cd_original,1)),1,60);%60是矩阵的维度
%并串转换
Tx_data_original=reshape(Tx_window_original,[],1);
%AWGN信道
 Ber_original=zeros(1,length(SNR));
 Ber2_original=zeros(1,length(SNR));
 
for ii=1:length(SNR)
    channe_original=awgn(Tx_data_original,SNR(ii),'measured');%添加高斯噪声
    % 串并转换
    Rx_data1_original=reshape(channe_original,N_fft+N_cp,[]);
    
    %去保护间隔、循环前缀
    Rx_data2_original=Rx_data1_original(N_cp+1:end,:);
    % FFT
    fft_data_original=fft(Rx_data2_original); 
    % 并串转换
    data_aftereq_original=reshape(fft_data_original,[],1);
    data_aftereq_original=data_aftereq_original(1:2652);
    data_aftereq_original=reshape(data_aftereq_original,N_sc,length(data_aftereq_original)/N_sc);
    
    %QPSK解调
    demodulation_data_original=pskdemod(data_aftereq_original,M,pi/M);    
    De_data1_original= reshape(demodulation_data_original,[],1);
    De_data2_original= de2bi(De_data1_original);
    De_Bit_original= reshape(De_data2_original',1,[]);
    
    % 信道译码（维特比译码）
    trellis_original = poly2trellis(7,[133 171]);
    rx_c_de_original = vitdec(De_Bit_original,trellis_original,tblen,'trunc','hard');   %硬判决
    
    %误码率
     [err,Ber2_original(jj)] = biterr(De_Bit_original,code_data(1:5304));%译码前的误码率
     [err, Ber_original(ii)] = biterr(rx_c_de_original,P_data(1:2652));%译码后的误码率
   
  
end

%% CFO估计(单径信道的频率估计和补偿)
FreqOffset=5*10^3/(5*10^9);%计算频偏比 载波频率为2GHz,频率偏移为5KHz

snr = 10^(-20/10);
snr1=2:2:20;
for jj=1:length(SNR)
channe_CFO=awgn(Tx_data,SNR(jj),'measured');
Rx = exp(i*2*pi*FreqOffset*(0:length(channe_CFO)-1)./N_fft).*channe_CFO ; %加频率偏移的信号

PHI_sum = zeros(1,Nd*(N_fft+N_cp)-N_fft);
GM_sum = zeros(1,Nd*(N_fft+N_cp)-N_fft);
%ML方法对于频偏估计,符号同步
for n = 64:Nd*(N_fft+N_cp)-(N_fft+N_cp)
    PHI=0;GM=0;
    for m = n:n+N_cp-1    
        PHI = PHI+ (Rx(m)*conj(Rx(m)) + Rx(m+N_fft)*conj(Rx(m+N_fft)));
        GM = GM+ Rx(m)*conj(Rx(m+N_fft));    
    end
    PHI_sum(n) = abs(GM)- (snr/(+1))*PHI;
    GM_sum(n) = -angle(GM)/(2*pi);
end

%确定频偏
[d_ml,deltaf_ml]=ml_estimate(Rx,1/4,N_fft,SNR);
max_peak=find_peak(d_ml,1/4,N_fft);
deltaf_ml1=FreqOffset;
for k=1:length(max_peak)
    if max_peak(k)~=0
        deltaf_ml1=deltaf_ml(k);
    end
    deltaf_ml(k)=deltaf_ml1;
end
deltaf_ml_add=zeros(1,(length(Rx)-length(deltaf_ml)));
deltaf_ml_add(1:end)=deltaf_ml(end);
deltaf_ml=cat(2,deltaf_ml,deltaf_ml_add);
%频率补偿
Rx_offset=Rx*exp(-i*2*pi*deltaf_ml(k)*(k-1)/N_fft);%频率补偿后的信号
%并串转换
Rx_data1_offset=reshape(Rx_offset,N_fft+N_cp,[]);
%去保护间隔
 Rx_data2_offset=Rx_data1_offset(N_cp+1:end,:);
%FFT
 fft_data_offset=fft(Rx_data2_offset);
%信道均衡
    data3_offset=fft_data_offset(1:N_fft,:); 
    Rx_pilot_offset=data3_offset(P_f_station(1:end),:); %接收到的导频
    h_offset=Rx_pilot_offset./pilot_seq; 
    H_offset=interp1( P_f_station(1:end)',h_offset,data_station(1:end)','linear','extrap');
%信道矫正
    data_aftereq_offset=data3_offset(data_station(1:end),:)./H_offset
%并串转换
    data_aftereq_offset=reshape(data_aftereq_offset,[],1);
    data_aftereq_offset=data_aftereq_offset(1:length(spread_data));
    data_aftereq_offset=reshape(data_aftereq_offset,N_sc,length(data_aftereq_offset)/N_sc);
%解扩
    demspread_data_offset = despread(data_aftereq_offset,code); 
%QPSK解调
    demodulation_data_offset=pskdemod(demspread_data_offset,M,pi/M);    
    De_data1_offset = reshape(demodulation_data_offset,[],1);
    De_data2_offset = de2bi(De_data1_offset);
    De_Bit_offset = reshape(De_data2_offset',1,[]);
%信道译码
    trellis_offset  = poly2trellis(7,[133 171]);
    rx_c_de_offset  = vitdec(De_Bit_offset ,trellis_offset ,tblen,'trunc','hard');   %硬判决
%计算误码率
     [err,Ber2_offset(jj)] = biterr(De_Bit_offset(1:length(code_data)),code_data);%译码前的误码率   
    [err, Ber_offset(jj)] = biterr(rx_c_de_offset(1:length(P_data)),P_data);%译码后的误码率
     
end 


%% 误比特率曲线
 %单径信道信道均衡1
 figure(5);
 subplot(3,1,1);
 semilogy(SNR,Ber2_single,'b-s');
 hold on;
 semilogy(SNR,Ber_single,'r-o');
 hold on;
 legend('4PSK调制、卷积码译码前（有扩频）','4PSK调制、卷积码译码后（信道均衡1）');
 hold on;
 xlabel('SNR');
 ylabel('BER');
 title('AWGN单径信道（无频偏补偿、信道均衡1）下误比特率曲线');

 %多径信道信道均衡1
 subplot(3,1,2);
 semilogy(SNR,Ber2_mult,'b-s');
 hold on;
 semilogy(SNR,Ber_mult,'r-o');
 hold on;
 legend('4PSK调制、卷积码译码前（有扩频）','4PSK调制、卷积码译码后（道信道均衡1）');
 hold on;
 xlabel('SNR');
 ylabel('BER');
 title('AWGN多径信道（无频偏补偿、信道均衡1）下误比特率曲线');
 
 %多径信道信道均衡2
 subplot(3,1,3);
 semilogy(SNR,Ber2_spline,'b-s');
 hold on;
 semilogy(SNR,Ber_spline,'r-o');
 hold on;
 legend('4PSK调制、卷积码译码前（扩频）','4PSK调制、卷积码译码后（信道均衡2）');
 hold on;
 xlabel('SNR');
 ylabel('BER');
 title('AWGN多径信道（无频偏补偿、信道均衡2）下误比特率曲线');
 
 figure(6);
 %画出频偏的图像
 subplot(3,1,1);plot(PHI_sum);title('定时偏移估计');grid on;
 subplot(3,1,2);plot(GM_sum);title('频率偏移估计');grid on;
 %频偏补偿
 subplot(3,1,3);
 semilogy(SNR,Ber2_offset,'b-s');
 hold on;
 semilogy(SNR,Ber_offset,'r-o');
 legend('4PSK调制、卷积码译码前（无行信道均衡）','4PSK调制、卷积码译码后（无行信道估计）');
 hold on;
 xlabel('SNR');
 ylabel('BER');
 title('AWGN单径信道下（频偏补偿、信道均衡1）误比特率曲线');
   



 %% 数据结果曲线
 figure(7)
 subplot(5,1,1);
 x=0:1:30;
 stem(x,P_data(1:31));
 ylabel('amplitude');
 title('发送数据（以前30个数据为例)');
 legend('4PSK调制、卷积译码、有扩频');

 subplot(5,1,2);
 x=0:1:30;
 stem(x,rx_c_de_single(1:31));
 ylabel('amplitude');
 title('接收数据单径信道(信道均衡1)');
 legend('4PSK调制、卷积译码、扩频、单径');
 
 subplot(5,1,3);
 x=0:1:30;
 stem(x,rx_c_de_mult(1:31));
 ylabel('amplitude');
 title('接收数据多径信道(信道均衡1)');
 legend('4PSK调制、卷积译码、扩频、多径');
 
 subplot(5,1,4);
 x=0:1:30;
 stem(x,rx_c_de_spline(1:31));
 ylabel('amplitude');
 title('接收数据多径信道(信道均衡2)');
 legend('4PSK调制、卷积译码、扩频、多径');
 
subplot(5,1,5);
 x=0:1:30;
 stem(x,rx_c_de_offset(1:31));
 ylabel('amplitude');
 title('接收数据单信道(信道均衡1、频率补偿)');
 legend('4PSK调制、卷积译码、有扩频、多径、有频率补偿');

 %% 星座图

  %单径 无均衡
 scatterplot(data_aftereq_original(:));
 title('AWGN单径无信道均衡接收信号的星座图'); 
 
 %单径 信道均衡1
 scatterplot(demspread_data_single(:));
 title('AWGN单径信道（信道均衡1）接收信号的星座图');
 
 %多径 信道均衡1
 scatterplot(demspread_data_mult(:));
 title('AWGN多径信道（信道均衡1）接收信号的星座图');

 %多径 信道均衡2
 scatterplot(demspread_data_spline(:));
 title('AWGN多径信道（信道均衡2）接收信号的星座图');
 
 %单径 信道均衡1 有频率补偿
 scatterplot(demspread_data_offset(:));
 title('AWGN单径信道（信道均衡1、有频率补偿）接收信号的星座图');
 

 
