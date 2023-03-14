%�����Ȼ���ƺ���ml_estimate
function [d_ml,deltaf_ml]=ml_estimate(prefix_data,g,fft_length,SNR)%���������루���ݣ�ѭ��ǰ׺������fft���ȣ�����ȣ��������ʱ���ƣ�Ƶƫ���ƣ�
len3=length(prefix_data);
for nc=1:len3-g*fft_length-fft_length
    for m=1:g*fft_length
        syn1(m)=(prefix_data(nc+m-1))*(conj(prefix_data(nc+m-1+fft_length)));
       sig_fi(m)=abs(prefix_data(nc+m-1))^2+abs(prefix_data(nc+m-1+fft_length))^2;
    end
    r(nc)=sum(syn1);
    fi(nc)=0.5*sum(sig_fi);
end
ro=SNR/(1+SNR);
d_ml=abs(r)-ro*fi;
deltaf_ml=-angle(r)/2/pi;