clc;
clear;
%ԭ��ο��ˣ�https://zhuanlan.zhihu.com/p/57967971


%% ��������

N_sc=52;      %ϵͳ���ز�����������ֱ���ز�����number of subcarrierA
N_fft=64;            % FFT ����
N_cp=16;             % ѭ��ǰ׺���ȡ�Cyclic prefix
N_symbo=N_fft+N_cp;        % 1������OFDM���ų���
N_c=53;             % ����ֱ���ز����ܵ����ز�����number of carriers
M=4;               %4PSK����
SNR=0:1:20;         %���������
N_frm=10;            % ÿ��������µķ���֡����frame
Nd=6;               % ÿ֡������OFDM������
P_f_inter=8;      %��Ƶ��� һ��8����Ƶ
data_station=[];    %��Ƶλ��
L=7;                %�����Լ������
tblen=6*L;          %Viterbi�������������
stage = 3;          % m���еĽ���
ptap1 = [1 3];      % m���еļĴ������ӷ�ʽ
regi1 = [1 1 1];    % m���еļĴ�����ʼֵ
alpha=0.3           %

%% �����������ݲ���
P_data=randi([0 1],1,N_sc*Nd*N_frm);
%������=����֡����ÿ֡��OFDM�����������ز���

%% �ŵ����루����롢��֯����
%����룺ǰ������������
%��֯��ʹͻ����������޶ȵķ�ɢ��
trellis = poly2trellis(7,[133 171]);       %(2,1,7)�������
code_data=convenc(P_data,trellis);


%% qpsk����
data_temp1= reshape(code_data,log2(M),[])';             %��ÿ��2���ؽ��з��飬M=4
data_temp2= bi2de(data_temp1);                             %������ת��Ϊʮ����
modu_data=pskmod(data_temp2,M,pi/M);              % 4PSK����
scatterplot(modu_data);
title('QPSK��������ͼ')             %����ͼ(Ҳ����ȡʵ����plot����)

%% ��Ƶ
%����������������������������������������������������������������������������������������������������������������%
%��Ƶͨ���ź���ռ�е�Ƶ�����Զ����������Ϣ�������С����
%������ũ������Ƶͨ�ž����ÿ�����似������ȡ������ϵĺô����������Ƶͨ�ŵĻ���˼����������ݡ�
%��Ƶ���ǽ�һϵ����������������������ź��ڻ�
%��Ƶ������Ƶ�ʱ����ԭ����m������Ƭ���� = 2����������* m����Ƶϵ����
%����������������������������������������������������������������������������������������������������������������%

code = mseq(stage,ptap1,regi1,N_sc);     % ��Ƶ�������
code = code * 2 - 1;         %��1��0�任Ϊ1��-1
modu_data=reshape(modu_data,N_sc,length(modu_data)/N_sc);
spread_data = spread(modu_data,code);        % ��Ƶ
spread_data=reshape(spread_data,[],1);

%% ���뵼Ƶ
P_f=3+3*1i;                       %Pilot frequency
P_f_station=[1:P_f_inter:N_fft];%��Ƶλ��
pilot_num=length(P_f_station);%��Ƶ����

for img=1:N_fft                        %����λ��
    if mod(img,P_f_inter)~=1          %mod(a,b)���������a����b������
        data_station=[data_station,img];
    end
end
data_row=length(data_station);
data_col=ceil(length(spread_data)/data_row);

pilot_seq=ones(pilot_num,data_col)*P_f;%����Ƶ�������
data=zeros(N_fft,data_col);%Ԥ����������
data(P_f_station(1:end),:)=pilot_seq;%��pilot_seq����ȡ

if data_row*data_col>length(spread_data)
    data2=[spread_data;zeros(data_row*data_col-length(spread_data),1)];%�����ݾ����룬��0������Ƶ~
else data2=spread_data;
end;
    

%% ����ת��
data_seq=reshape(data2,data_row,data_col);
data(data_station(1:end),:)=data_seq;%����Ƶ�����ݺϲ�

%% IFFT
ifft_data=ifft(data);
figure(2)
subplot(3,1,1);
plot(abs(ifft_data(:,1)),'b');
title('ԭʼ����OFDM����');
xlabel('Time');
ylabel('Amplitude');

%% ���뱣�������ѭ��ǰ׺
Tx_cd=[ifft_data(N_fft-N_cp+1:end,:);ifft_data];%��ifft��ĩβN_cp�������䵽��ǰ��
subplot(3,1,2);
plot(abs(Tx_cd(:,1)));
xlabel('Time');
ylabel('Amplitude');
title('��CP�ĵ���OFDM����');

%% �Ӵ�
Tx_window=zeros(size(Tx_cd));

% ͨ��������
Tx_window=Tx_cd.*repmat(rcoswindow(alpha,size(Tx_cd,1)),1,390);%390�Ǿ���Tx_cd��ά��(80*390)
subplot(3,1,3);
plot(abs(Tx_window(:,1)))
title('�Ӵ���ĵ���OFDM����')
xlabel('Time');
ylabel('Amplitude');

%% �ź�Ƶ��
%δ�Ӵ��ź�Ƶ��
%figure % ��һ��
figure(3);
orgin_aver_power = 20*log10(mean(abs(fft(Tx_cd'))));
subplot(2,1,1)
plot((1:length(orgin_aver_power))/length(orgin_aver_power),orgin_aver_power)
hold on
plot(0:1/length(orgin_aver_power):1 ,-35,'rd')
hold on
axis([0 1 -40 max(orgin_aver_power)])
grid on
title('δ�Ӵ��ź�Ƶ��')

%�Ӵ��ź�Ƶ��
orgin_aver_power = 20*log10(mean(abs(fft(Tx_window'))));
subplot(2,1,2)
plot((1:length(orgin_aver_power))/length(orgin_aver_power),orgin_aver_power)
hold on
plot(0:1/length(orgin_aver_power):1 ,-35,'rd')
hold on
axis([0 1 -40 max(orgin_aver_power)])
grid on
title('�Ӵ��ź�Ƶ��')

%% ����ת��
Tx_data=reshape(Tx_window,[],1);%���ڴ�����Ҫ

%% �ŵ����ŵ�����1���ŵ�����2��
 %Ԥ��������ʵı�������
 Ber_single=zeros(1,length(SNR)); %���������ʾ������� �ŵ�����1
 Ber2_single=zeros(1,length(SNR));
 
 Ber_mult=zeros(1,length(SNR));%�ྶ�����ʾ������� �ŵ�����1
 Ber2_mult=zeros(1,length(SNR));
 
 Ber_spline=zeros(1,length(SNR));%���������ʾ������� �ŵ�����2
 Ber2_spline=zeros(1,length(SNR));
 
 Ber_offset=zeros(1,length(SNR));  %�ྶ�����ʾ������� �ŵ�����2 ��Ƶ�ʲ���
 Ber2_offset=zeros(1,length(SNR));

Tx_data = reshape(Tx_window,1,[]); % ���ʱ��һ�������źţ�������
%signal_origin = reshape(Tx_window,1,[]); % δ�Ӵ������ź�
mult_path_am = [1 0.2 0.1]; %  �ྶ����
mult_path_time1=randi([1,50]); % �ྶ���ʱ��1-50
mult_path_time2=randi([1,50]);
path2 = 0.2*[zeros(1,mult_path_time1) Tx_data(1:end-mult_path_time1) ];
path3 = 0.1*[zeros(1,mult_path_time2) Tx_data(1:end-mult_path_time2) ];
Tx_mult = Tx_data + path2 + path3; % �ྶ�ź�

figure('menubar','none');
subplot(2,1,1);
plot(real(Tx_mult) )
title('�ྶ��OFDM�ź�')
xlabel('Time/samples')
ylabel('Amplitude')
subplot(2,1,2)
plot(real(Tx_data) )
title('������OFDM�ź�')
xlabel('Time/samples')
ylabel('Amplitude')

%% ����/�ྶ�ŵ�����

for jj=1:length(SNR)
    channe_single=awgn(Tx_data,SNR(jj),'measured');
    channe_mult=awgn(Tx_mult ,SNR(jj),'measured');

    %��Ӹ�˹������
    %y = awgn(x,SNR,SIGPOWER) 
    %���SIGPOWERΪ'measured'���������ڼ�������֮ǰ�ⶨ�ź�ǿ�ȡ�
%% ����ת��
    Rx_data1_single=reshape(channe_single,N_fft+N_cp,[]);
    Rx_data1_mult=reshape(channe_mult,N_fft+N_cp,[]);
    
%% ȥ�����������ѭ��ǰ׺
    Rx_data2_single=Rx_data1_single(N_cp+1:end,:);
    Rx_data2_mult=Rx_data1_mult(N_cp+1:end,:);

%% FFT
    fft_data_single=fft(Rx_data2_single);
    fft_data_mult=fft(Rx_data2_mult);
    
%% �ŵ��������ֵ���ԣ��ŵ�����1��
    %�����ŵ�
    data3_single=fft_data_single(1:N_fft,:); 
    Rx_pilot_single=data3_single(P_f_station(1:end),:); %���յ��ĵ�Ƶ
    h_single=Rx_pilot_single./pilot_seq; 
    H_single=interp1( P_f_station(1:end)',h_single,data_station(1:end)','linear','extrap');
    %�ֶ����Բ�ֵ����ֵ�㴦����ֵ�����������ڽ������������Ժ���Ԥ�⡣
    %�Գ�����֪�㼯�Ĳ�ֵ����ָ����ֵ����(����)���㺯��ֵ
   
    %�ྶ�ŵ�
    data3_mult=fft_data_mult(1:N_fft,:); 
    Rx_pilota_mult=data3_mult(P_f_station(1:end),:); %���յ��ĵ�Ƶ
    h_mult=Rx_pilota_mult./pilot_seq; 
    H_mult=interp1( P_f_station(1:end)',h_mult,data_station(1:end)','linear','extrap');
    
   %% ��;���ŵ�������ֵ���ƣ��ŵ�����2��
     H_mult_spline=interp1( P_f_station(1:end)',h_mult,data_station(1:end)','spline','extrap');
     
%% �ŵ�У��
    data_aftereq_single=data3_single(data_station(1:end),:)./H_single;
    data_aftereq_mult=data3_mult(data_station(1:end),:)./H_mult;
    
    data_aftereq_spline=data3_mult(data_station(1:end),:)./H_mult_spline;
    
%% ����ת��
    %�����ŵ�(�ŵ����ⷽ��1)
    data_aftereq_single=reshape(data_aftereq_single,[],1);
    data_aftereq_single=data_aftereq_single(1:length(spread_data));
    data_aftereq_single=reshape(data_aftereq_single,N_sc,length(data_aftereq_single)/N_sc);
    %�ྶ�ŵ�(�ŵ����ⷽ��1)
    data_aftereq_mult=reshape(data_aftereq_mult,[],1);
    data_aftereq_mult=data_aftereq_mult(1:length(spread_data));
    data_aftereq_mult=reshape(data_aftereq_mult,N_sc,length(data_aftereq_mult)/N_sc);
    %�ྶ�ŵ�(�ŵ����ⷽ��2)
    data_aftereq_spline=reshape(data_aftereq_spline,[],1);
    data_aftereq_spline=data_aftereq_spline(1:length(spread_data));
    data_aftereq_spline=reshape(data_aftereq_spline,N_sc,length(data_aftereq_spline)/N_sc);    
    
%% ����
    demspread_data_single = despread(data_aftereq_single,code);       % ���ݽ��� 
    demspread_data_mult = despread(data_aftereq_mult,code); 
   
    demspread_data_spline = despread(data_aftereq_spline,code); 

%% QPSK���
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
%% ���⽻֯��
%% �ŵ����루ά�ر����룩
    trellis_single = poly2trellis(7,[133 171]);
    rx_c_de_single = vitdec(De_Bit_single,trellis_single,tblen,'trunc','hard');   %Ӳ�о�
    
    trellis_mult = poly2trellis(7,[133 171]);
    rx_c_de_mult = vitdec(De_Bit_mult,trellis_mult,tblen,'trunc','hard');   %Ӳ�о�
    
    trellis_spline = poly2trellis(7,[133 171]);
    rx_c_de_spline= vitdec(De_Bit_spline,trellis_spline,tblen,'trunc','hard');   %Ӳ�о�
%% ����������
    [err,Ber2_single(jj)] = biterr(De_Bit_single(1:length(code_data)),code_data);%����ǰ��������
    [err, Ber_single(jj)] = biterr(rx_c_de_single(1:length(P_data)),P_data);%������������
    
    [err,Ber2_mult(jj)] = biterr(De_Bit_mult(1:length(code_data)),code_data);%����ǰ��������
    [err, Ber_mult(jj)] = biterr(rx_c_de_mult(1:length(P_data)),P_data);%������������
    
    [err,Ber2_spline(jj)] = biterr(De_Bit_spline(1:length(code_data)),code_data);%����ǰ��������
    [err, Ber_spline(jj)] = biterr(rx_c_de_spline(1:length(P_data)),P_data);%������������   
    
    
end

%% ԭʼ�źŲ�ͨ���ŵ�����
%�ŵ����롢QPSK���ƺ���ź�Ϊmodu_data
%IFFT
ifft_data_original=ifft(modu_data);
%���뱣�������ѭ��ǰ׺
Tx_cd_original=[ifft_data_original(N_fft-N_cp+1:end,:);ifft_data_original];
%�Ӵ� ��ͨ�������ˣ�
Tx_window_original=zeros(size(Tx_cd_original));
Tx_window_original=Tx_cd_original.*repmat(rcoswindow(alpha,size(Tx_cd_original,1)),1,60);%60�Ǿ����ά��
%����ת��
Tx_data_original=reshape(Tx_window_original,[],1);
%AWGN�ŵ�
 Ber_original=zeros(1,length(SNR));
 Ber2_original=zeros(1,length(SNR));
 
for ii=1:length(SNR)
    channe_original=awgn(Tx_data_original,SNR(ii),'measured');%��Ӹ�˹����
    % ����ת��
    Rx_data1_original=reshape(channe_original,N_fft+N_cp,[]);
    
    %ȥ���������ѭ��ǰ׺
    Rx_data2_original=Rx_data1_original(N_cp+1:end,:);
    % FFT
    fft_data_original=fft(Rx_data2_original); 
    % ����ת��
    data_aftereq_original=reshape(fft_data_original,[],1);
    data_aftereq_original=data_aftereq_original(1:2652);
    data_aftereq_original=reshape(data_aftereq_original,N_sc,length(data_aftereq_original)/N_sc);
    
    %QPSK���
    demodulation_data_original=pskdemod(data_aftereq_original,M,pi/M);    
    De_data1_original= reshape(demodulation_data_original,[],1);
    De_data2_original= de2bi(De_data1_original);
    De_Bit_original= reshape(De_data2_original',1,[]);
    
    % �ŵ����루ά�ر����룩
    trellis_original = poly2trellis(7,[133 171]);
    rx_c_de_original = vitdec(De_Bit_original,trellis_original,tblen,'trunc','hard');   %Ӳ�о�
    
    %������
     [err,Ber2_original(jj)] = biterr(De_Bit_original,code_data(1:5304));%����ǰ��������
     [err, Ber_original(ii)] = biterr(rx_c_de_original,P_data(1:2652));%������������
   
  
end

%% CFO����(�����ŵ���Ƶ�ʹ��ƺͲ���)
FreqOffset=5*10^3/(5*10^9);%����Ƶƫ�� �ز�Ƶ��Ϊ2GHz,Ƶ��ƫ��Ϊ5KHz

snr = 10^(-20/10);
snr1=2:2:20;
for jj=1:length(SNR)
channe_CFO=awgn(Tx_data,SNR(jj),'measured');
Rx = exp(i*2*pi*FreqOffset*(0:length(channe_CFO)-1)./N_fft).*channe_CFO ; %��Ƶ��ƫ�Ƶ��ź�

PHI_sum = zeros(1,Nd*(N_fft+N_cp)-N_fft);
GM_sum = zeros(1,Nd*(N_fft+N_cp)-N_fft);
%ML��������Ƶƫ����,����ͬ��
for n = 64:Nd*(N_fft+N_cp)-(N_fft+N_cp)
    PHI=0;GM=0;
    for m = n:n+N_cp-1    
        PHI = PHI+ (Rx(m)*conj(Rx(m)) + Rx(m+N_fft)*conj(Rx(m+N_fft)));
        GM = GM+ Rx(m)*conj(Rx(m+N_fft));    
    end
    PHI_sum(n) = abs(GM)- (snr/(+1))*PHI;
    GM_sum(n) = -angle(GM)/(2*pi);
end

%ȷ��Ƶƫ
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
%Ƶ�ʲ���
Rx_offset=Rx*exp(-i*2*pi*deltaf_ml(k)*(k-1)/N_fft);%Ƶ�ʲ�������ź�
%����ת��
Rx_data1_offset=reshape(Rx_offset,N_fft+N_cp,[]);
%ȥ�������
 Rx_data2_offset=Rx_data1_offset(N_cp+1:end,:);
%FFT
 fft_data_offset=fft(Rx_data2_offset);
%�ŵ�����
    data3_offset=fft_data_offset(1:N_fft,:); 
    Rx_pilot_offset=data3_offset(P_f_station(1:end),:); %���յ��ĵ�Ƶ
    h_offset=Rx_pilot_offset./pilot_seq; 
    H_offset=interp1( P_f_station(1:end)',h_offset,data_station(1:end)','linear','extrap');
%�ŵ�����
    data_aftereq_offset=data3_offset(data_station(1:end),:)./H_offset
%����ת��
    data_aftereq_offset=reshape(data_aftereq_offset,[],1);
    data_aftereq_offset=data_aftereq_offset(1:length(spread_data));
    data_aftereq_offset=reshape(data_aftereq_offset,N_sc,length(data_aftereq_offset)/N_sc);
%����
    demspread_data_offset = despread(data_aftereq_offset,code); 
%QPSK���
    demodulation_data_offset=pskdemod(demspread_data_offset,M,pi/M);    
    De_data1_offset = reshape(demodulation_data_offset,[],1);
    De_data2_offset = de2bi(De_data1_offset);
    De_Bit_offset = reshape(De_data2_offset',1,[]);
%�ŵ�����
    trellis_offset  = poly2trellis(7,[133 171]);
    rx_c_de_offset  = vitdec(De_Bit_offset ,trellis_offset ,tblen,'trunc','hard');   %Ӳ�о�
%����������
     [err,Ber2_offset(jj)] = biterr(De_Bit_offset(1:length(code_data)),code_data);%����ǰ��������   
    [err, Ber_offset(jj)] = biterr(rx_c_de_offset(1:length(P_data)),P_data);%������������
     
end 


%% �����������
 %�����ŵ��ŵ�����1
 figure(5);
 subplot(3,1,1);
 semilogy(SNR,Ber2_single,'b-s');
 hold on;
 semilogy(SNR,Ber_single,'r-o');
 hold on;
 legend('4PSK���ơ����������ǰ������Ƶ��','4PSK���ơ������������ŵ�����1��');
 hold on;
 xlabel('SNR');
 ylabel('BER');
 title('AWGN�����ŵ�����Ƶƫ�������ŵ�����1���������������');

 %�ྶ�ŵ��ŵ�����1
 subplot(3,1,2);
 semilogy(SNR,Ber2_mult,'b-s');
 hold on;
 semilogy(SNR,Ber_mult,'r-o');
 hold on;
 legend('4PSK���ơ����������ǰ������Ƶ��','4PSK���ơ����������󣨵��ŵ�����1��');
 hold on;
 xlabel('SNR');
 ylabel('BER');
 title('AWGN�ྶ�ŵ�����Ƶƫ�������ŵ�����1���������������');
 
 %�ྶ�ŵ��ŵ�����2
 subplot(3,1,3);
 semilogy(SNR,Ber2_spline,'b-s');
 hold on;
 semilogy(SNR,Ber_spline,'r-o');
 hold on;
 legend('4PSK���ơ����������ǰ����Ƶ��','4PSK���ơ������������ŵ�����2��');
 hold on;
 xlabel('SNR');
 ylabel('BER');
 title('AWGN�ྶ�ŵ�����Ƶƫ�������ŵ�����2���������������');
 
 figure(6);
 %����Ƶƫ��ͼ��
 subplot(3,1,1);plot(PHI_sum);title('��ʱƫ�ƹ���');grid on;
 subplot(3,1,2);plot(GM_sum);title('Ƶ��ƫ�ƹ���');grid on;
 %Ƶƫ����
 subplot(3,1,3);
 semilogy(SNR,Ber2_offset,'b-s');
 hold on;
 semilogy(SNR,Ber_offset,'r-o');
 legend('4PSK���ơ����������ǰ�������ŵ����⣩','4PSK���ơ����������������ŵ����ƣ�');
 hold on;
 xlabel('SNR');
 ylabel('BER');
 title('AWGN�����ŵ��£�Ƶƫ�������ŵ�����1�������������');
   



 %% ���ݽ������
 figure(7)
 subplot(5,1,1);
 x=0:1:30;
 stem(x,P_data(1:31));
 ylabel('amplitude');
 title('�������ݣ���ǰ30������Ϊ��)');
 legend('4PSK���ơ�������롢����Ƶ');

 subplot(5,1,2);
 x=0:1:30;
 stem(x,rx_c_de_single(1:31));
 ylabel('amplitude');
 title('�������ݵ����ŵ�(�ŵ�����1)');
 legend('4PSK���ơ�������롢��Ƶ������');
 
 subplot(5,1,3);
 x=0:1:30;
 stem(x,rx_c_de_mult(1:31));
 ylabel('amplitude');
 title('�������ݶྶ�ŵ�(�ŵ�����1)');
 legend('4PSK���ơ�������롢��Ƶ���ྶ');
 
 subplot(5,1,4);
 x=0:1:30;
 stem(x,rx_c_de_spline(1:31));
 ylabel('amplitude');
 title('�������ݶྶ�ŵ�(�ŵ�����2)');
 legend('4PSK���ơ�������롢��Ƶ���ྶ');
 
subplot(5,1,5);
 x=0:1:30;
 stem(x,rx_c_de_offset(1:31));
 ylabel('amplitude');
 title('�������ݵ��ŵ�(�ŵ�����1��Ƶ�ʲ���)');
 legend('4PSK���ơ�������롢����Ƶ���ྶ����Ƶ�ʲ���');

 %% ����ͼ

  %���� �޾���
 scatterplot(data_aftereq_original(:));
 title('AWGN�������ŵ���������źŵ�����ͼ'); 
 
 %���� �ŵ�����1
 scatterplot(demspread_data_single(:));
 title('AWGN�����ŵ����ŵ�����1�������źŵ�����ͼ');
 
 %�ྶ �ŵ�����1
 scatterplot(demspread_data_mult(:));
 title('AWGN�ྶ�ŵ����ŵ�����1�������źŵ�����ͼ');

 %�ྶ �ŵ�����2
 scatterplot(demspread_data_spline(:));
 title('AWGN�ྶ�ŵ����ŵ�����2�������źŵ�����ͼ');
 
 %���� �ŵ�����1 ��Ƶ�ʲ���
 scatterplot(demspread_data_offset(:));
 title('AWGN�����ŵ����ŵ�����1����Ƶ�ʲ����������źŵ�����ͼ');
 

 
