%% Extract Epochs of different classes from main GDF file
addpath('/Users/jimmyjarjoura/Downloads/eeglab14_0_0b/')
run eeglab
fprintf('press any hey to continue\n');
pause();
%% Check mapping of events to numbers in eeglab: EDIT -> Event values or ALLEEG.event variable 
test_out= pop_epoch( ALLEEG(1), {'775','1041','1042','1043','1044','773'}, [0 12],'newname','test');

%% EDiting event types


for i=1:size(test_out.event,2)
            switch test_out.event(i).type
                case 775
           test_out.event(i).type='baseline';
                 case 1042
           test_out.event(i).type='fear';
                 case 1041
           test_out.event(i).type='angry';
                case 1043
           test_out.event(i).type='happy';
                case 773
           test_out.event(i).type='tender';
                case 1044
           test_out.event(i).type='sad';
            end
end
       
    % mood_class = num2str(mood_class);
    % mood_class(mood_class=='1042') = 'fear';
    % mood_class(mood_class=='1041') = 'angry';
    % mood_class(mood_class=='1043') = 'happy';
    % mood_class(mood_class=='773') = 'tender'
    % mood_class(mood_class=='1044') = 'sad';
        
 
testing=eeg_epoch2continuous(test_out);

%% Specifiying band pass filters


d_delta = fdesign.bandpass('N,F3dB1,F3dB2',4,1,4,128);
Delta = design(d_delta,'butter');

d_theta = fdesign.bandpass('N,F3dB1,F3dB2',4,5,8,128);
Theta = design(d_theta,'butter');

d_alpha = fdesign.bandpass('N,F3dB1,F3dB2',4,8,12,128);
Alpha = design(d_alpha,'butter');

d_beta = fdesign.bandpass('N,F3dB1,F3dB2',4,12,30,128);
Beta = design(d_beta,'butter');

d_gamma = fdesign.bandpass('N,F3dB1,F3dB2',4,31,39,128);
Gamma = design(d_gamma,'butter');



sr=128

%% Epoch settings
ws_in = input('input window size (secs) of the epochs (default 4secs):');
hs_in = input('input hop size (secs) of the epochs (default 1secs):');

if isempty(ws_in)
    ws=4;%window in seconds
else
    ws=ws_in;
end
if isempty(hs_in)
    hs=1;%window in seconds
else
    hs=hs_in;
end
%        ws2=20;%historic window
w=floor(ws*sr);%in samples...
%w2=floor(ws2*sr);
h=floor(hs*sr);
len=size(test_out.data,2);
FDlen=floor((len-w)/h)+1;

%% Filtering data through 5 bands, Alpha, delta,theta,gamma and beta
for j=1:size(test_out.data,3)
         for i=1:size(test_out.data,1)
         alpha_data(i,:,j)= filter(Alpha,test_out.data(i,:,j));
         delta_data(i,:,j)= filter(Delta,test_out.data(i,:,j));
         theta_data(i,:,j)= filter(Theta,test_out.data(i,:,j));
         gamma_data(i,:,j)= filter(Gamma,test_out.data(i,:,j));
         beta_data(i,:,j)= filter(Beta,test_out.data(i,:,j));
         end
end

%% Computing band power using the Band_power function
for j=1:size(test_out.data,3)
    for i=1:size(test_out.data,1)   
        
     Bp=Band_power(delta_data(i,:,j),ws,hs);
     Bp_delta(1:numel(Bp),i,j)=Bp;
     
     Bp=Band_power(theta_data(i,:,j),ws,hs);
     Bp_theta(1:numel(Bp),i,j)=Bp;
     
     Bp=Band_power(alpha_data(i,:,j),ws,hs);
     Bp_alpha(1:numel(Bp),i,j)=Bp;
     
     Bp=Band_power(beta_data(i,:,j),ws,hs);
     Bp_beta(1:numel(Bp),i,j)=Bp;
     
     Bp=Band_power(gamma_data(i,:,j),ws,hs);
     Bp_gamma(1:numel(Bp),i,j)=Bp;
    end
end

Bp_delta_concat=reshape(permute(Bp_delta,[2,1,3]),size(Bp_delta,2),[])';
Bp_theta_concat=reshape(permute(Bp_theta,[2,1,3]),size(Bp_theta,2),[])';
Bp_alpha_concat=reshape(permute(Bp_alpha,[2,1,3]),size(Bp_alpha,2),[])';
Bp_beta_concat=reshape(permute(Bp_beta,[2,1,3]),size(Bp_beta,2),[])';
Bp_gamma_concat=reshape(permute(Bp_gamma,[2,1,3]),size(Bp_gamma,2),[])';


Bp_global= [Bp_delta_concat Bp_theta_concat Bp_alpha_concat Bp_beta_concat Bp_gamma_concat];


%% Prefrontal assymetry difference(F3,F4) of all bands
delta_prefront_diff=[];
theta_prefront_diff=[];
alpha_prefront_diff=[];
beta_prefront_diff=[];
gamma_prefront_diff=[];

delta_prefront_diff_AF=[];
theta_prefront_diff_AF=[];
alpha_prefront_diff_AF=[];
beta_prefront_diff_AF=[];
gamma_prefront_diff_AF=[];

for i=1:size(Bp_alpha,3)
    delta_prefront_diff=cat(1,delta_prefront_diff,Bp_delta(:,12,i)-Bp_delta(:,3,i));
    theta_prefront_diff=cat(1,theta_prefront_diff,Bp_theta(:,12,i)-Bp_theta(:,3,i));
    alpha_prefront_diff=cat(1,alpha_prefront_diff,Bp_alpha(:,12,i)-Bp_alpha(:,3,i));
    beta_prefront_diff=cat(1,beta_prefront_diff,Bp_beta(:,12,i)-Bp_beta(:,3,i));
    gamma_prefront_diff=cat(1,gamma_prefront_diff,Bp_gamma(:,12,i)-Bp_gamma(:,3,i));
    

    % Prefrontal assymetry difference(AF3,AF4) of all bands
    delta_prefront_diff_AF=cat(1,delta_prefront_diff_AF,Bp_delta(:,14,i)-Bp_delta(:,1,i));
    theta_prefront_diff_AF=cat(1,theta_prefront_diff_AF,Bp_theta(:,14,i)-Bp_theta(:,1,i));
    alpha_prefront_diff_AF=cat(1,alpha_prefront_diff_AF,Bp_alpha(:,14,i)-Bp_alpha(:,1,i));
    beta_prefront_diff_AF=cat(1,beta_prefront_diff_AF,Bp_beta(:,14,i)-Bp_beta(:,1,i));
    gamma_prefront_diff_AF=cat(1,gamma_prefront_diff_AF,Bp_gamma(:,14,i)-Bp_gamma(:,1,i));

end

Bp_prefront_global= [delta_prefront_diff theta_prefront_diff alpha_prefront_diff beta_prefront_diff gamma_prefront_diff delta_prefront_diff_AF theta_prefront_diff_AF alpha_prefront_diff_AF beta_prefront_diff_AF gamma_prefront_diff_AF ];




%% Generate Column headers to be included in excel csv
l=1;
    for j=1:14
    myfilename = sprintf('Delta_Channel%d',j);
    mydata{l} = myfilename;
    myfilename = sprintf('Theta_Channel%d',j);
    mydata{l+14} = myfilename;
    myfilename = sprintf('Alpha_Channel%d',j);
    mydata{l+2*14} = myfilename;
    myfilename = sprintf('Beta_Channel%d',j);
    mydata{l+3*14} = myfilename;
    myfilename = sprintf('Gamma_Channel%d',j);
    mydata{l+4*14} = myfilename;
    l=l+1;
    end
   
l=1;
    myfilename = sprintf('Delta_pref');
    mydata_pref{l} = myfilename;
    myfilename = sprintf('Theta_pref');
    mydata_pref{l+1} = myfilename;
    myfilename = sprintf('Alpha_pref');
    mydata_pref{l+2} = myfilename;
    myfilename = sprintf('Beta_pref');
    mydata_pref{l+3} = myfilename;
    myfilename = sprintf('Gamma_pref');
    mydata_pref{l+4} = myfilename;
    myfilename = sprintf('Delta_pref_AF');
    mydata_pref{l+5} = myfilename;
    myfilename = sprintf('Theta_pref_AF');
    mydata_pref{l+6} = myfilename;
    myfilename = sprintf('Alpha_pref_AF');
    mydata_pref{l+7} = myfilename;
    myfilename = sprintf('Beta_pref_AF');
    mydata_pref{l+8} = myfilename;
    myfilename = sprintf('Gamma_pref_AF');
    mydata_pref{l+9} = myfilename;
    
 
    
   %% Generating mood tags depending on window and hop size 
   mood_class={};
    for i=1:size(test_out.event,2)
            switch test_out.event(i).type
                case 'baseline'
           names = repmat({'baseline'}, FDlen, 1);
           mood_class=cat(1,mood_class,names);
                 case 'fear'
           names = repmat({'fear'}, FDlen, 1);
           mood_class=cat(1,mood_class,names);
                 case 'angry'
            names = repmat({'angry'}, FDlen, 1);
           mood_class=cat(1,mood_class,names);
                case 'happy'
           names = repmat({'happy'}, FDlen, 1);
           mood_class=cat(1,mood_class,names);
                case 'tender'
           names = repmat({'tender'}, FDlen, 1);
           mood_class=cat(1,mood_class,names);
                case 'sad'
           names = repmat({'sad'}, FDlen, 1);
           mood_class=cat(1,mood_class,names);
            end
    end
   
    a='mood';
    mood=cat(1,a,mood_class);
    
    
% mood_class = num2str(mood_class);
% mood_class(mood_class=='1042') = 'fear';
% mood_class(mood_class=='1041') = 'angry';
% mood_class(mood_class=='1043') = 'happy';
% mood_class(mood_class=='773') = 'tender'
% mood_class(mood_class=='1044') = 'sad';
%     
    

%% Extract chanels AF3 AF4 F3 and F4


  for i=1:size(test_out.data,3)
    Af3(i,1:numel(test_out.data(1,:,i)))=test_out.data(1,:,i);
    Af4(i,1:numel(test_out.data(14,:,i)))=test_out.data(14,:,i);
    F3(i,1:numel(test_out.data(3,:,i)))=test_out.data(3,:,i);
    F4(i,1:numel(test_out.data(12,:,i)))=test_out.data(12,:,i); 
  end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Band waves filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Spatial filter
for i=1:size(Af3,1)
    data(i,1:numel(Af3(i,:)))=0.25*Af3(i,:)+0.25*Af4(i,:)+0.25*F3(i,:)+0.25*F4(i,:);
end

% data(1:numel(Af3))=0.25*Af3+0.25*Af4+0.25*F3+0.25*F4;
% len=length(data);%length of EEG sample
% sr=128;%eeg sample rate hz

%% Band pass filter (Alpha): 
% IIR Butterworth filter of order 4 with 3?dB frequencies of 8 and 12 Hz. The sampling frequency is 128 Hz
% d_alpha = fdesign.bandpass('N,F3dB1,F3dB2',4,8,12,128);
% Alpha = design(d_alpha,'butter');

for i=1:size(data,1)
    alpha_data_Assym(i,:)= filter(Alpha,data(i,:));
    alpha_F3(i,:)= filter(Alpha,F3(i,:));
    alpha_F4(i,:)= filter(Alpha,F4(i,:));
end

%% Band pass filter (Beta): 
% IIR Butterworth filter of order 4 with 3?dB frequencies of 12 and 30 Hz. The sampling frequency is 128 Hz
% d_beta = fdesign.bandpass('N,F3dB1,F3dB2',4,12,30,128);
% Beta = design(d_beta,'butter');

for i=1:size(data,1)
    beta_data_Assym(i,:)= filter(Beta,data(i,:));
    beta_F3(i,:)= filter(Beta,F3(i,:));
    beta_F4(i,:)= filter(Beta,F4(i,:));
end

%% Arousal and Valence calculations

%changed Ar=zeros(4,FDlen) to Ar=zeros(1,FDlen) in arousal_sig
% and Ar_F3=zeros(1,FDlen)and Ar_F4=zeros(1,FDlen) as second,third and
% fourth row were not used. and to be able to include the information into
% a single array.
for i=1:size(alpha_data_Assym,1)
    %calculate arousal over all the signal w=seg, h=1seg...
    Arr(i,:)=Arousal_sig(alpha_data_Assym(i,:),beta_data_Assym(i,:),ws,hs);
    %scale Arousal
    Arrs(i,:)=scale1_1(max(Arr(i,:)),min(Arr(i,:)),Arr(i,:));
    %% Valence calculations
    %calculate valence over all the signal w=seg, h=1seg...
    Val(i,:)=Valence_sig(alpha_F3(i,:),alpha_F4(i,:),beta_F3(i,:),beta_F4(i,:),ws,hs);
    %scale Valence
    Vals(i,:)=scale1_1(max(Val(i,:)),min(Val(i,:)),Val(i,:));
end

Arrs_reshaped=reshape(Arrs.',1,[]);
Vals_reshaped=reshape(Vals.',1,[]);



%% Preparing Finale file    
%% Alfeatures concatenated
All_features=[Bp_global Bp_prefront_global transpose(Arrs_reshaped) transpose(Vals_reshaped)]
header = cat(2, mydata , mydata_pref);
header= cat(2,header,'Arousal');
header= cat(2,header,'Valence');
Finale = [header;num2cell(All_features)];
OUT=cat(2,Finale,mood);
%% Writing to csv
cell2csv('/Users/jimmyjarjoura/Desktop/Final_project/experiment3/EEG_experiment3.csv',OUT)


% %% PLOTS
% %% Arousal over all the signal and mean of each time slot
% figure(2);
% 
%     %     title('EEG signal');
%     %     xlabel('Seconds');
%     
% %     plot([1:length(Arrs(i,:))]*(h/sr),Arrs(i,:));
%     plot([1:length(Arrs_reshaped)]*(h/sr),Arrs_reshaped);
% 
%     title('Normalized Arousal and Interval Average');
%     xlabel('Seconds');
%     ylabel('Arousal Normalized');
%     FDlen=floor((len-w)/h)+1;
% for j=1:FDlen:(length(Arrs_reshaped)-FDlen)
%     meanArrs(1,j)=mean(Arrs_reshaped(j:(j-1)+FDlen));
%     hline(meanArrs(1,j),j,j+FDlen,'k');
%     vline(j,'r');
% end
% figure(3); 
% 
% plot([1:length(Vals_reshaped)]*(h/sr),Vals_reshaped);
% 
%     title('Normalized Valence and Interval Average');
%     xlabel('Seconds');
%     ylabel('Valence Normalized');
%     FDlen=floor((len-w)/h)+1;
% for j=1:FDlen:(length(Vals_reshaped)-FDlen)
%     meanVals(1,j)=mean(Vals_reshaped(j:(j-1)+FDlen));
%     hline(meanVals(1,j),j,j+FDlen,'k');
%     vline(j,'r');
% end
% for i=1:size(Vals,1)
%     subplot(size(Vals,1),1,i),
%     plot([1:length(Vals(i,:))]*(h/sr),Vals(i,:));
%     hline(mean(Vals(i,:)),0,length(Vals(i,:)),'r');
%     title('Normalized Valence and Interval Average');
%     xlabel('Seconds');
%     ylabel('Valence Normalized');
% end
