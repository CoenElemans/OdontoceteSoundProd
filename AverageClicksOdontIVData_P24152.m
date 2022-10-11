
% Odontocete In Vitro sound production
% This script computes several features of filmed porpoise phonic lip opening and closing and 
% signatures in acoustic recordings (hydrophone/micrphone) in individual
% porpoises.

% The script computes the average upsampled Digital Kymogram (DKG) of DKG snippets
% alligned to hydophone click onset detections. 

% Author: Coen P.H. Elemans
% Date: 22/10/10

clc
clear all
close all

% Dropbox
rootdir=pwd;

HI1_sens    =   1e6*10^(-211/20); % Hydrophone sensitivity [V/Pa];

% Microphone parameters
Data.time_Op_SI1= [];
Data.SI1_peak_mag= [];
Data.SI1_peak_time= [];
Data.SI1_diff_min= [];
Data.SI1_diff_min_time=[];
Data.SI1_diff_max= [];
Data.SI1_diff_max_time=[];

% Hydrophone parameters
Data.NumClicks= [];
Data.time_Click_HI1= [];
Data.HI1_RMS_peak_mag= [];
Data.HI1_RMS_peak_time= [];

% Bursa time and kinematics
Data.Bursa_Time_Op= [];
Data.Bursa_Time_Cl= [];
Data.Bursa_MaxOp= [];

% Pressure/Flow
Data.Psub= [];
Data.Fsub= [];
Data.Numclicks=[];

% Plot style
FontName = 'Arial';
FontSize = 18;
FontWeight = 'Normal';
FontColor = 'k';
FontName = 'Arial';
LabelFontSizeMultiplier=1.5;
colmap = lines(7);


%% Cycle over individuals
plotfig=false; 

for k = [1]
        
    % Hydrophone parameters
    time_Click_HI1= [];
    HI1_RMS_peak_mag= [];
    HI1_RMS_peak_time= [];

    % Bursa time and kinematics
    Bursa_Time_Op= [];
    Bursa_Time_Cl= [];
    Bursa_MaxOp= [];

    % Pressure/Flow
    Psub= [];
    Fsub= [];

    switch k
        
        % PHOTRON FILE
        case 1 % P24152 015             
            Data_dir = [ rootdir '/P24152_R_2_015'];
            WAV_filename = 'P24152_R_2_64_1.0V(peak)__30-01-2020_00_07_56.wav';
            DAQ_filename = 'Data_P24152_R_2_015#001.mat'; 
            DKG_filename = 'P24152_R_C001H001S0010_S0001_DKG_PX348.333_PY50.233.mat';
            HI1_sens    =   1e6*10^(-211/20); % Hydrophone sensitivity [V/Pa];
              
            H1Ch=2; CTICh=1; % Channel for hydrophone and trigger
            
            NThr    =   0.0005/HI1_sens;   % Threshold for noise removal [a.u.]
            F_LP   =   [100e3]    ;   % Low pass filter to detect switch noise
            SegDel  =   600/2; 	% Remove segment [points] 
                
            SI1_Thr     =   1;            % Threshold [Pa] to detect openings in microphone signal
            SI1_holdoff =   1.2e-3;           % Hold off after detection [s]
            SI1_segmentsize = [2e-3 2e-3];  % Segment before and after threshold crossing
           
            DKG_segmentsize_NoiseRem = [-1 2];  % Number of frames before and after threshold crossing to check for power switch noise
            DKG_segmentsize = [-20 40];          % Number of frames before and after threshold crossing
               
            DetTresh_SI1_snippet=-5e-1; % [Pa] Detection of opening in microphone snippet 
            DetTresh_HI1_snippet= 200;  % [Pa] Detection of opening in hydrophone snippet 
            Opthr=30;   %38                % [intensity] Detection threshold of opening in DKG
                      
            d=.08;  % [m] distance to microphone       
            H1dist=.1;
            St_Click=1;                
            SplitInd=1;         
    end

     
    %% Load Data
    cd (Data_dir);
    'Reading files and conditioning'
    % Load WAV from file
    [y, Fs_wav]= audioread(WAV_filename);
    HI1     =	y(:,H1Ch)/HI1_sens;    % Hydrophone 1 [Pa]
    CTI2    =   y(:,CTICh);            % Cam Trig
    time_wav=   [0:length(HI1)-1]/Fs_wav;
    trig_wav=   time_wav(find(CTI2>0.7,1,'first')); % Find trigger
    time_wav=   time_wav'-trig_wav;
    Hlag    =   round((H1dist/1500)*Fs_wav);  % number of samples to shift due to TOF microphone signal

    % Load DAQ data
    load (DAQ_filename)
    Fs_daq  =   1/(time(2)-time(1));
    trig_daq=   time(find(CTI>0.7,1,'first'))
    time_daq=   time-trig_daq;
    lag      =   round((d/340)*Fs_daq);  % number of samples to shift due to TOF microphone signal

    % band pass filter hydrophone 
    Wn       =   [2e3 180e3]/(Fs_wav/2);
    [B,A]    =   butter(3,Wn,'bandpass');
    HI1      =   filtfilt(B, A, HI1);
    
    % band pass filter hydrophone  to isolate clicks
    Wn       =   [100e3 180e3]/(Fs_wav/2);
    [B,A]    =   butter(3,Wn,'bandpass');
    HI1_fil  =   filtfilt(B, A, HI1);
    HI1(1:Hlag-1)=[]; HI1=[HI1 ; zeros(Hlag-1,1)]; % Offset time of flight lag sound in melon
    HI1_fil(1:Hlag-1)=[]; HI1_fil=[HI1_fil ; zeros(Hlag-1,1)]; % Offset time of flight lag sound in melon

    HI1_RMSbackgroundnoiselevel=rms(HI1_fil(1:round(.1*Fs_wav)));

    % band pass filter microphone
    Wn       =   [500 10000]/(Fs_daq/2);
    [B,A]    =   butter(3,Wn,'bandpass');
    SI1      =   filtfilt(B, A, SI1);
    SI1(1:lag-1)=[]; SI1=[SI1 ; zeros(lag-1,1)]; % Offset time of flight lag sound in air
    % low pass filter pressure and flow
    Wn       =   [100]/(Fs_daq/2);
    [B,A]    =   butter(3,Wn,'low');
    FI1      =   filtfilt(B, A, FI1);    
    FI2      =   filtfilt(B, A, FI2);    
    PI1      =   filtfilt(B, A, PI1);    
    
    % Load DKG (single one)
    load(DKG_filename);
    NumFrames=length(DKGMetaData.time_vid);
        
    % Find indices of WAV and DAQ files per frame. 
    DAQ_ind = zeros(1, length(DKGMetaData.time_vid));
    WAV_ind = zeros(1, length(DKGMetaData.time_vid));

    DAQ_ind(1)=find(time_daq>=DKGMetaData.time_vid(1), 1, 'first');
    WAV_ind(1)=find(time_wav>=DKGMetaData.time_vid(1), 1, 'first');
        
    for j = 2:NumFrames
    	DAQ_ind(j)=find(time_daq(DAQ_ind(j-1):DAQ_ind(j-1)+1000)>=DKGMetaData.time_vid(j), 1, 'first')+DAQ_ind(j-1)-1;
    	WAV_ind(j)=find(time_wav(WAV_ind(j-1):WAV_ind(j-1)+1000)>=DKGMetaData.time_vid(j), 1, 'first')+WAV_ind(j-1)-1;
    end
    DAQ_NSampleFr   =  round(mean(diff(DAQ_ind)))-1; % Average number of DAQ samples per Video frame
    WAV_NSampleFr   =  round(mean(diff(WAV_ind)))-1; % Average number of WAV samples per Video frame
    DAQ_Plot_Ind    =  DAQ_ind(1):DAQ_ind(NumFrames)+DAQ_NSampleFr; % Indices for plotting signals over entire DKG
    WAV_Plot_Ind    =  WAV_ind(1):WAV_ind(NumFrames)+WAV_NSampleFr;
    
%     Compute RMS per frame for all pressure and flow signals
    PI1_rms = zeros(1, length(DKGMetaData.time_vid));    
    PIM1_rms= zeros(1, length(DKGMetaData.time_vid));
    PIM2_rms= zeros(1, length(DKGMetaData.time_vid));
    FI1_rms = zeros(1, length(DKGMetaData.time_vid));    
    FI2_rms = zeros(1, length(DKGMetaData.time_vid));    
    HI1_rms = zeros(1, length(DKGMetaData.time_vid)); 
    t_rms   = zeros(1, length(DKGMetaData.time_vid)); 
    
    for j = 1:NumFrames-1
    	t_rms(j)    =   mean(time_wav(DAQ_ind(j):DAQ_ind(j+1)-1));
    	PI1_rms(j)  =   mean(PI1(DAQ_ind(j):DAQ_ind(j+1)-1));
    	PIM1_rms(j) =   mean(PIM1(DAQ_ind(j):DAQ_ind(j+1)-1));
    	PIM2_rms(j) =   mean(PIM2(DAQ_ind(j):DAQ_ind(j+1)-1));
    	FI1_rms(j)  =   mean(FI1(DAQ_ind(j):DAQ_ind(j+1)-1));
    	FI2_rms(j)  =   mean(FI2(DAQ_ind(j):DAQ_ind(j+1)-1));
        
        HI1_rms(j)  =   rms(HI1_fil(WAV_ind(j):WAV_ind(j+1)-1));
    end    
    
    SI1_diff=[diff(SI1); 0];
    
    %% Plot overview figure with data
  
    figure(1) 
    cla, hold on
    plot (time_wav(WAV_Plot_Ind), HI1_fil(WAV_Plot_Ind), 'Color', colmap(1,:))
    plot (DKGMetaData.time_vid, HI1_rms, 'k-')
    ylabel ('Hydrophone [Pa]')
    grid on
    title (WAV_filename)
    xlim([.0 .25])

    figure(2)
    cla
    imagesc([DKGMetaData.time_vid(1) DKGMetaData.time_vid(NumFrames)],[0 size(DKG_Data,1)/DKGMetaData.scale_len], DKG_Data)
    ylabel ('DKG line [mm]')
    grid on
    colormap gray
    xlim([.0 .25])

    figure(3)
    hold on
    plot (time_daq(DAQ_Plot_Ind), SI1(DAQ_Plot_Ind),'Color', colmap(2,:))
    ylabel ('Microphone [Pa]')
    title (DAQ_filename)
    grid on
    xlim([.0 .25])

    figure(4)
    yyaxis left
    plot (DKGMetaData.time_vid, PI1_rms, 'Color', colmap(4,:))
    ylabel ('Pressure [Pa]')
    ylim ([6 8])         
    yyaxis right
    plot (DKGMetaData.time_vid, FI1_rms+FI2_rms, 'Color', colmap(6,:))
    ylabel ('Flow [LPM]')
    ylim ([5000 6000]) 
    title (DAQ_filename)
    grid on
    xlabel ('Time [s]')
    xlim([.0 .25])
      
    %% Detect clicks in SI1 and HI1

% Detect openings on microphone signal SI1
    SI1_holdoff_ind =   round((SI1_holdoff)*Fs_daq);

    above_thr       =   find((SI1(DAQ_Plot_Ind)) > SI1_Thr);
    diffvector      =   diff(above_thr);
    SI_pulses_ind   =   DAQ_Plot_Ind([above_thr(1); above_thr(find(diffvector > SI1_holdoff_ind)+1)]);
    N_SIpulses      =   length(SI_pulses_ind);
    IPI             =   diff(time_daq(SI_pulses_ind));     % Find Inter pulse interval [s]
    SI1_reprate     =   1./IPI;            % Repetition rate equals 1/ICI 
      
% Detect clicks on hydrophone signal HI1
    HI1_holdoff_ind =   round((SI1_holdoff)*Fs_wav);
    [yupper,ylower] =   envelope (HI1_fil, round(.1e-3*Fs_wav), 'rms');

    above_thr       =   find((yupper(WAV_Plot_Ind)) > 2*HI1_RMSbackgroundnoiselevel);
    diffvector      =   diff(above_thr);
    HI1_pulses_ind   =   WAV_Plot_Ind([above_thr(1); above_thr(find(diffvector > HI1_holdoff_ind)+1)]);
    N_HI1pulses      =   length(HI1_pulses_ind);

    IPI             =   diff(time_wav(HI1_pulses_ind));     % Find Inter pulse interval [s]
    HI1_reprate         =   1./IPI;            % Repetition rate equals 1/ICI 
    if plotfig
        figure (1)
        plot(time_wav(HI1_pulses_ind),0 , 'k.') % detections
    end

        
%% Remove clicks with SI1 switch noise or NaNs

    Rem_Nan=zeros(1,N_SIpulses);
    
    for j=1:N_SIpulses-1;

        % clear HI1_snippet Video_ind_snippet SI1_snippet time_daq_snippet time_wav_snippet
        % find video frame boundaries for snippet
        Video_ind_snippet   =   DKG_segmentsize_NoiseRem+find(DKGMetaData.time_vid>=time_daq(SI_pulses_ind(j)), 1, 'first');
        % Isolate snippet 
        HI1_snippet         =   HI1(WAV_ind(Video_ind_snippet(1)):WAV_ind(Video_ind_snippet(2)));
           
        % if segment contains NANs then skip
        if sum(isnan(HI1_snippet))>0
            continue 
        end
        Rem_Nan(j)=1;
    end
    
    SI_pulses_ind_NoiseRem=SI_pulses_ind(find(Rem_Nan));  

    %% Pre allocate mean DKG  
    Fs_interp            =  100000;
    DKG_mean_segmentsize =  [3e-3 5e-3];
    DKG_click_allignInd  =  round(DKG_mean_segmentsize(1)*Fs_interp);
    DKG_click_mean       =  NaN*zeros(length(SI_pulses_ind_NoiseRem)-1, size(DKG_Data,1), round(sum(DKG_mean_segmentsize)*Fs_interp));
    SI1_click_mean       =  NaN*zeros(length(SI_pulses_ind_NoiseRem)-1, round(sum(DKG_mean_segmentsize)*Fs_interp));
    HI1_click_mean       =  NaN*zeros(length(SI_pulses_ind_NoiseRem)-1, round(sum(DKG_mean_segmentsize)*Fs_interp));


    %% Isolate remaining clicks
    clc
    ['Of #' num2str(N_SIpulses) ' pulses in SI1, #' num2str(length(SI_pulses_ind_NoiseRem)) ' pulses contain no noise switch in HI1']
    
    % Signal Conditioning; band pass filter microphone
    Wn       =   [4e3]/(Fs_daq/2);
    [B,A]    =   butter(3,Wn,'low');
    SI1_fil  =   filtfilt(B, A, SI1);
    FI1      =   filtfilt(B, A, FI1);
    FI2      =   filtfilt(B, A, FI2);
    PIM1     =   filtfilt(B, A, PIM1);

    t1=1e-3; % segment to determine SI1 offset

    clear time_Op_SI1 time_Click_HI1 Psub Fsub HI1_RMS_peak_mag HI1_RMS_peak_time 
    clear Bursa_MaxOp Bursa_Time_Cl Bursa_Time_Op SI1_diff_max SI1_diff_min
    for j=1:30% length(SI_pulses_ind_NoiseRem)-1;    % when plotting only averaged click
%     for j=St_Click: length(SI_pulses_ind_NoiseRem)-1;   % for all clicks  

        j
        % find video frame boundaries for snippet
        Video_ind_snippet   =   round(DKG_segmentsize+find(DKGMetaData.time_vid>=time_daq(SI_pulses_ind_NoiseRem(j)), 1, 'first'));
        
        % Isolate snippet 
        HI1_snippet         =   HI1_fil(WAV_ind(Video_ind_snippet(1)):WAV_ind(Video_ind_snippet(2))-1);  
        Spec_snippet         =   HI1(WAV_ind(Video_ind_snippet(1)):WAV_ind(Video_ind_snippet(2))-1);  
        SI1_snippet         =   SI1(DAQ_ind(Video_ind_snippet(1)):DAQ_ind(Video_ind_snippet(2)));
        SI1_fil_snippet     =   SI1_fil(DAQ_ind(Video_ind_snippet(1)):DAQ_ind(Video_ind_snippet(2)));
        time_daq_snippet    =   time_daq(DAQ_ind(Video_ind_snippet(1)):DAQ_ind(Video_ind_snippet(2)));
        time_wav_snippet    =   time_wav(WAV_ind(Video_ind_snippet(1)):WAV_ind(Video_ind_snippet(2))-1);
        time_DKG_snippet    =   [DKGMetaData.time_vid(Video_ind_snippet(1)): 1/DKGMetaData.Fs_vid: DKGMetaData.time_vid(Video_ind_snippet(2))-1/DKGMetaData.Fs_vid];
        
        % Pressure parameters
        Psub(j)= mean(PIM1_rms(Video_ind_snippet(1):Video_ind_snippet(2)-1));
        Fsub(j)= mean(FI1_rms(Video_ind_snippet(1):Video_ind_snippet(2)-1))+ mean(FI2_rms(Video_ind_snippet(1):Video_ind_snippet(2)-1));
        
        % Microphone signal parameters
        SI1_offset=mean(SI1(1:round(t1*Fs_daq)));
       
        % Detect opening with threshold
        DKG_snip = DKG_Data(:,Video_ind_snippet(1):Video_ind_snippet(2)-1);
        BW=(DKG_snip<Opthr);
        Opening= (sum(BW,1))/DKGMetaData.scale_len;
        [Opening_vel,Opening_acc,Lpath]   =   pos2velacc(Opening,zeros(1,length(Opening)),DKGMetaData.Fs_vid);
            
        % Interpolate DKG image
        F = griddedInterpolant(DKG_snip);
        [sx,sy] = size(DKG_snip);
        xq = (1:.25:sx)';
        yq = (1:sy)';
        F.Method = 'spline';
        vq = (F({xq,yq}));
        BW=(vq<Opthr);
        Opening_interp= .25*(sum(BW,1))/DKGMetaData.scale_len;
        [Opening_vel_interp,Opening_acc_interp,Lpath]   =   pos2velacc(Opening,zeros(1,length(Opening)),DKGMetaData.Fs_vid); % Kinematic data

        % Find opening from microphone signal
        if isempty(time_daq_snippet(find(SI1_fil_snippet<DetTresh_SI1_snippet, 1, 'first')));
            time_Op_SI1(j)=NaN;
            SI1_peak_mag(j)=NaN;
            SI1_peak_time(j)=NaN;
            SI1_diff_max(j)=NaN;
            SI1_diff_max_time(j)=NaN;
            SI1_diff_min(j)=NaN;
            SI1_diff_min_time(j)=NaN;
            time_Click_HI1(j)= NaN;
            HI1_RMS_peak_mag(j)=NaN;
            HI1_RMS_peak_time(j)=NaN;
            HI1_peakind=[];
            Bursa_MaxOp(j)=NaN; Bursa_Time_Cl(j)=NaN;  Bursa_Time_Op(j)=NaN;            
            'no SI1 threshold crossed'
            continue
        else
            time_Op_SI1(j)=time_daq_snippet(find(SI1_fil_snippet<DetTresh_SI1_snippet, 1, 'first')); % time when signal crosses negative threshold
            [SI1_peak_mag(j) ind]=max(SI1_fil_snippet);
            SI1_peak_time(j)=time_daq_snippet(ind);
            [SI1_diff_max(j) ind] =max(diff(SI1_snippet)); % max of derivative of sound signal
            SI1_diff_max_time(j) = time_daq_snippet(ind); 
            [SI1_diff_min(j) ind]=min(diff(SI1_snippet)); % max of derivative of sound signal
            SI1_diff_min_time(j) = time_daq_snippet(ind); 
        end
   
        % Hydrophone parameters
        [yupper,ylower]=envelope (HI1_snippet, round(.1e-3*Fs_wav), 'rms');
        
        if isempty(find(yupper>2*HI1_RMSbackgroundnoiselevel, 1, 'first'))% check if click is above 2x HI1_RMSbackgroundnoiselevel
            time_Click_HI1(j)= NaN;
            [HI1_RMS_peak_mag(j) ]=max(yupper);
            HI1_RMS_peak_time(j)=NaN;
            HI1_peakind=[];
            Bursa_MaxOp(j)=NaN; Bursa_Time_Cl(j)=NaN;  Bursa_Time_Op(j)=NaN;
            'no HI1 threshold crossed'

        else
            time_Click_HI1(j)=time_wav_snippet(find(yupper>2*HI1_RMSbackgroundnoiselevel, 1, 'first')); % time when rsm HI1 signal crosses 2x noise level 
            [HI1_RMS_peak_mag(j) HI1_peakind]=max(yupper); % max RMS HI1
            HI1_RMS_peak_time(j)=time_wav_snippet(HI1_peakind);
            [Bursa_MaxOp(j) ind] =max(Opening_interp);            
              
            % Find opening and closing of bursa
            switch 'MaxSI1'
                case 'MaxSI1' % Search for first closure before and after maximal SI1 aplitude
                    ind=find(time_DKG_snippet>=SI1_peak_time(j), 1, 'first')-SplitInd;
                case 'MaxOp' % Search for first closure before and after maximal opening
            end
            
            if ~isempty(find(Opening_interp(ind:end)==0, 1, 'first' ))
                Bursa_Time_Cl(j) =  time_DKG_snippet(ind-1+find(Opening_interp(ind:end)==0, 1, 'first' ));                
            else
                Bursa_Time_Cl(j) =  time_DKG_snippet(end);
            end
            if ~isempty(find(fliplr(Opening_interp(1:ind))==0, 1, 'first'))
                Bursa_Time_Op(j) =            time_DKG_snippet(ind+1-find(fliplr(Opening_interp(1:ind))==0, 1, 'first') );
            else
                Bursa_Time_Op(j) =            time_DKG_snippet(1);
            end   
            
        end
        

        % Plot individual click
        if plotfig
            figure(5)
            cla  
            nfft=200;
            noverlap=floor(.95*nfft);
            w   =   window(@flattopwin,nfft);
            [Sp,F,T,P] = spectrogram(Spec_snippet,w, noverlap,nfft, Fs_wav);  % calc spectrogram data
            surf(T,F,10*log10(abs(P)),'EdgeColor','none');
            axis xy; colormap((jet)); view(0,90);
            ylim([20e3 200e3])
            colorbar off
            ylabel('Frequency (Hz)');
            caxis([-10 30])


            figure(6)
            yyaxis left
            cla
            hold on
            plot(time_wav_snippet, abs(HI1_snippet), 'Color', colmap(1,:))          
            plot(time_Click_HI1(j), 2*HI1_RMSbackgroundnoiselevel, 'ro' ) % HI1 theshold crossing      
            plot(time_wav_snippet, yupper, 'k-')
            ylabel ('Hydrophone [Pa]')
            ylim ([0 1400])
            
            yyaxis right
            cla
            plot(time_DKG_snippet, Opening_interp, 'Color',colmap(3,:))
            hold on
            plot(Bursa_Time_Op(j), 0, 'kv') 
            plot(Bursa_Time_Cl(j), 0, 'kv') 
            ylim ([0 .30])
            ylabel ('Opening [mm]')
            grid on
            if isnan(time_Click_HI1(j))
                 title ([ DAQ_filename '\color{black}  SI1 click nr: ' num2str(j) ', NO HI1 click detected'] ) 
            else
                 title ([DAQ_filename '\color{red} SI1 click nr: ' num2str(j) ', HI1 click detected'] ) 
            end

            figure(7)
            cla
            imagesc([DKGMetaData.time_vid(Video_ind_snippet(1)) DKGMetaData.time_vid(Video_ind_snippet(2))],[0 size(DKG_Data,1)/DKGMetaData.scale_len], DKG_snip)
            hold on
%             Make a truecolor all-green image.
            red = cat(3, colmap(3,1)*ones(size(BW)), colmap(3,2)*ones(size(BW)), colmap(3,3)*ones(size(BW)));
            h = imagesc([DKGMetaData.time_vid(Video_ind_snippet(1)) DKGMetaData.time_vid(Video_ind_snippet(2))],[0 size(DKG_Data,1)/DKGMetaData.scale_len],red);
            hold off
            % Use influence map image as the AlphaData for the solid red image.
            set(h, 'AlphaData', BW)
            colormap((gray))
            ylabel ('DKG line [mm]')
            xlabel('Time (s)');
            title ([DAQ_filename '\color{red} SI1 click nr: ' num2str(j) ', HI1 click detected'] ) 
            drawnow
            pause

        end
    
        
    % Compute average DKG over many clicks, allign to peak in microphone signal        
    [SI1_snippet_resample] = resample(SI1_snippet,time_daq_snippet, Fs_interp);
    [HI1_snippet_rms_resample] = resample(yupper,time_wav_snippet, Fs_interp);
    [ MinVal MaxInd]=    min(diff(SI1_snippet_resample)); % allign to minimal velocity of SI1 signal
    SI1_click_mean(j, DKG_click_allignInd-MaxInd:DKG_click_allignInd-MaxInd+length(SI1_snippet_resample)-1) = SI1_snippet_resample;
    HI1_click_mean(j, DKG_click_allignInd-MaxInd:DKG_click_allignInd-MaxInd+length(HI1_snippet_rms_resample)-1) = HI1_snippet_rms_resample;
    
    % Reshape DKG_snip so that its repeated 5 or 10x to Fs_interp 
    DKG_stretch=[];
    for k=1:Fs_interp/DKGMetaData.Fs_vid
        DKG_stretch=[ DKG_stretch; DKG_snip ];
    end
    AB = reshape([DKG_stretch], size(DKG_snip,1), []);

    DKG_click_mean(j, :, DKG_click_allignInd-MaxInd: DKG_click_allignInd-MaxInd+size(AB,2)-1) = AB; 

    sum(~isnan(time_Click_HI1)&~isnan(time_Op_SI1));  
    sum(~isnan(HI1_RMS_peak_mag)); % all clicks detected on hydrophone
    
    end
    
    
%% Average performance at higher sampling rate
DKG_snip_mean=squeeze(mean(DKG_click_mean, 1,'omitnan'));
DKG_snip_mean(isnan(DKG_snip_mean))=255;

BW=(DKG_snip_mean<42);
Opening_mean= (sum(BW,1))/DKGMetaData.scale_len;

t_interp=[1:size(HI1_click_mean,2)]/Fs_interp*1000;

fig=figure(8);
ax3(1)=subplot(3,1,1);
plot(t_interp, mean(HI1_click_mean,1, 'omitnan'), 'Color', colmap(1,:))
hold on
plot(t_interp, mean(HI1_click_mean,1, 'omitnan')+std(HI1_click_mean,0,1, 'omitnan'), ':', 'Color', colmap(1,:))
plot(t_interp, mean(HI1_click_mean,1, 'omitnan')-std(HI1_click_mean,0,1, 'omitnan'), ':', 'Color', colmap(1,:))
grid on
% xlabel('Time (index)', 'FontName', FontName, 'FontSize', FontSize,  'FontWeight', FontWeight, 'Color', FontColor);
ylabel('Hydrophone (Pa)', 'FontName', FontName, 'FontSize', FontSize,  'FontWeight', FontWeight, 'Color', FontColor);
ylim([0 800])

% Interpolate DKG image
F = griddedInterpolant(DKG_snip_mean);
[sx,sy] = size(DKG_snip_mean);
xq = (0:.25:sx)';
yq = (0:sy)';
F.Method = 'spline';
vq = (F({xq,yq}));
title('Higher Resolution')

ax3(2)=subplot(3,1,2:3);
imagesc(t_interp,[], vq)
colormap gray
caxis([0 100])
grid on
xlabel('Time (ms)', 'FontName', FontName, 'FontSize', FontSize,  'FontWeight', FontWeight, 'Color', FontColor);
ylabel('DKG mean', 'FontName', FontName, 'FontSize', FontSize,  'FontWeight', FontWeight, 'Color', FontColor);
linkaxes(ax3,'x')
xlim([1 4.5])
ylim([300 500])

% fig.Renderer='painters';
% saveas(gcf,'P24152_R_2_015_averagedclick1_35_221003.svg');



%% Append data  
    Ind = [];
    Ind2=zeros(size(Bursa_Time_Op)); Ind2(St_Click:end)=1; % Include clicks from start index on
    Ind= ~isnan(time_Op_SI1) & Ind2; 
        
    % Microphone parameters
    Data.SI1_diff_max= [Data.SI1_diff_max , SI1_diff_max(Ind)];  
    Data.SI1_diff_max_time= [Data.SI1_diff_max_time , SI1_diff_max_time(Ind)];  
    Data.SI1_diff_min= [Data.SI1_diff_min , SI1_diff_min(Ind)];    
    Data.SI1_diff_min_time= [Data.SI1_diff_min_time , SI1_diff_min_time(Ind)];  
    Data.time_Op_SI1= [Data.time_Op_SI1 , time_Op_SI1(Ind)];    
    Data.SI1_peak_mag= [Data.SI1_peak_mag , SI1_peak_mag(Ind)];    
    Data.SI1_peak_time= [Data.SI1_peak_time , SI1_peak_time(Ind)];    

    % Hydrophone parameters
    Data.time_Click_HI1= [Data.time_Click_HI1 , time_Click_HI1(Ind)];
    Data.HI1_RMS_peak_mag= [Data.HI1_RMS_peak_mag , HI1_RMS_peak_mag(Ind)];
    Data.HI1_RMS_peak_time= [Data.HI1_RMS_peak_time , HI1_RMS_peak_time(Ind)];
            
    % Bursa time and kinematics
    Data.Bursa_Time_Op= [Data.Bursa_Time_Op , Bursa_Time_Op(Ind)];
    Data.Bursa_Time_Cl= [Data.Bursa_Time_Cl , Bursa_Time_Cl(Ind)];
    Data.Bursa_MaxOp= [Data.Bursa_MaxOp , Bursa_MaxOp(Ind)];

    % Pressure/Flow
    Data.Psub= [Data.Psub , Psub(Ind)];
    Data.Fsub= [Data.Fsub , Fsub(Ind)];
    Data.Numclicks = [Data.Numclicks sum(Ind)];
  
    
end

%% 
IPI = Data.time_Click_HI1(2:end)-Data.time_Click_HI1(1:end-1);
ClickRate=1./IPI;

Ind= ~isnan(Data.time_Click_HI1); % Find indices with HI1 clicks

figure % histogram
cla
histogram(1000*(Data.Bursa_Time_Op-Data.time_Click_HI1),'BinWidth',.05);
pd = fitdist(1000*(Data.Bursa_Time_Op-Data.time_Click_HI1)','Normal')
hold on
histogram(1000*(Data.Bursa_Time_Cl-Data.time_Click_HI1),'BinWidth',.05);   
pd = fitdist(1000*(Data.Bursa_Time_Cl-Data.time_Click_HI1)','Normal')
    
xlabel('Delay (ms)', 'FontName', FontName, 'FontSize', FontSize,  'FontWeight', FontWeight, 'Color', FontColor);
ylabel('Number in bin', 'FontName', FontName, 'FontSize', FontSize,  'FontWeight', FontWeight, 'Color', FontColor);
legend ('Bursa Opening - Click', ' Bursa Closing - Click')
title ([ DAQ_filename   '  Nclicks: ' num2str(numel(Data.time_Click_HI1)) ])
xlim ([-2.5 .5])
ax = gca;ax.FontName = FontName; ax.FontSize; ax.LabelFontSizeMultiplier=LabelFontSizeMultiplier;

%%
figure 
subplot 121  % Timing SI1 parameters to bursa state
histogram(1000*(Data.time_Op_SI1-Data.Bursa_Time_Op),'BinWidth',.05);
xlim([-.5 .5])
xlabel('Delay (ms)', 'FontName', FontName, 'FontSize', FontSize,  'FontWeight', FontWeight, 'Color', FontColor);
ylabel('Number in bin', 'FontName', FontName, 'FontSize', FontSize,  'FontWeight', FontWeight, 'Color', FontColor);
legend ('SI1 Opening - Bursa Opening')
title ([ DAQ_filename   '  Nclicks: ' num2str(numel(Data.time_Click_HI1)) ])
ax = gca;ax.FontName = FontName; ax.FontSize; ax.LabelFontSizeMultiplier=LabelFontSizeMultiplier;

subplot 122
histogram(1000*(Data.SI1_diff_min_time-Data.Bursa_Time_Cl),'BinWidth',.05);
xlim([-.5 .5])
xlabel('Delay (ms)', 'FontName', FontName, 'FontSize', FontSize,  'FontWeight', FontWeight, 'Color', FontColor);
ylabel('Number in bin', 'FontName', FontName, 'FontSize', FontSize,  'FontWeight', FontWeight, 'Color', FontColor);
legend ('min diff(SI1) - Bursa Closing')
title ([ DAQ_filename   '  Nclicks: ' num2str(numel(Data.time_Click_HI1)) ])
ax = gca;ax.FontName = FontName; ax.FontSize; ax.LabelFontSizeMultiplier=LabelFontSizeMultiplier;
