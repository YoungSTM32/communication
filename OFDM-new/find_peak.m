%Ѱ�Ҷ�ʱ��ֵλ�ú���
function max_peak=find_peak(syn2,g,fft_length)%���������루ML��ʱ��ѭ��ǰ׺������fft���ȣ��������ֵλ�ü���ȷ��λ��
max_peak=imregionalmax(syn2);%�ҵ������ط�
a_max=1:length(syn2);
position_max=a_max(max_peak);
data_max=syn2(position_max);
%figure(8)
%stem(position_max,data_max)

max_peak=zeros(1,length(syn2));
for i=1:length(position_max)
    a1_max=position_max(i)+(g*fft_length+fft_length)/2;
    a1_min=position_max(i)-(g*fft_length+fft_length)/2;
    a2=find(position_max<a1_max&position_max>a1_min);%λ��
    a31=a2(1);
    a32=a2(end);
    a4=data_max(a31:a32);
    a5=max(a4);
    if data_max(i)==a5
        max_peak(position_max(i))=data_max(i);
    end
end
