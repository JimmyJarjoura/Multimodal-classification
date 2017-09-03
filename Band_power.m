function Bp=Band_power(data,ws,hs)
%w in secs
%h in secs
len=length(data);%lenght of EEG sample
sr=128;%eeg sample rate hz
w=floor(ws*sr);%window size
h=floor(hs*sr);%hop size
pin=1;%point in 
pout=w;%point out
FDlen=floor((len-w)/h)+1;%calculate leght of the FD vecctor
Bp=zeros(1,FDlen);%create empty vector of arousal values, pin and pout and pmid

    for i=1:FDlen;%for each sliding window do
        a=data(pin:pout);%get band data
%         a2=log(1+mean(a.^2));%calculate logaritmic power 
        a2=log(mean(a.^2));%calculate logaritmic power
        Bp(1,i)=a2
        pin=pin+h;%advance window by hop size(not overlap!!!)
        pout=pout+h;
    end
    
end
 