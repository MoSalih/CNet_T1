clc; clear; close all

nfiles = 16;

for c = 0:15
%     c = 0;
    lction1 = '\\ntapprdfs01n02.rmit.internal\eh4\e29974\Configuration\Desktop\Randoms\DNN\Databases\IEMOCAP\IEMOCAP_full_release\Session1\sentences\wav\Ses01F_impro01\';
    if c < 10
        lction2 = 'Ses01F_impro01_F00';
    else
        lction2 = 'Ses01F_impro01_F0';
    end
    lction3 = num2str(c);
    lction4 = '.wav';
    lction = [lction1 lction2 lction3 lction4];
    
    [wav,fs] = audioread(lction);
    fftw = fft(wav);
    
    N = length(wav);
    ts = 1/fs;
    fstep = fs/N;
%     nf = fs/2;
    t = 0:ts:(N*ts)-ts;
    f = 0:fstep:fs-fstep;
    
    Tw = 25;                % analysis frame duration (ms)
    Ts = 10;                % analysis frame shift (ms)
    alpha = 0.97;           % preemphasis coefficient
    M = 20;                 % number of filterbank channels 
    C = 39;                 % number of cepstral coefficients
    L = 22;                 % cepstral sine lifter parameter
    LF = 200;               % lower frequency limit (Hz)
    HF = 3600;              % upper frequency limit (Hz)

%     figure;
%     plot(t,wav);
%     xlabel('Time(s)'); 
%     ylabel('Amplitude');
%     
%     figure (c + 1 + (nfiles*1))
%     plot(f,real(fftshift(fftw)));
%     xlabel('Frequency(Hz)'); 
%     ylabel('Amplitude');
    
    i = 0;
    j = 0;
    while j < 0.025
        
        i = i + 1;
        j = t(i);
        if j == 0.01
            inc = i;
        end
        
    end
    
    win_lng = i;
    nfft = 2^nextpow2(i);
    
    [spec_full,freq,time] = spectrogram(wav,win_lng,inc,nfft,fs,'yaxis');
    [ wav_mfcc, FBEs, frames ] = mfcc( wav, fs, Tw, Ts, alpha, @hamming, [LF HF], M, C+1, L );

    [ Nw, NF ] = size( frames );                % frame length and number of frames
    time_frames = [0:NF-1]*Ts*0.001+0.5*Nw/fs;  % time vector (s) for frames 
    time = [ 0:length(wav)-1 ]/fs;           % time vector (s) for signal samples 
    logFBEs = 20*log10( FBEs );                 % compute log FBEs for plotting
    logFBEs_floor = max(logFBEs(:))-50;         % get logFBE floor 50 dB below max
    logFBEs( logFBEs<logFBEs_floor ) = logFBEs_floor; % limit logFBE dynamic range
    
    
    figure;
    
    subplot(2,1,1);
%     spectrogram(wav,win_lng,inc,nfft,fs,'yaxis'); colorbar('off');
%     colormap('jet');
    imagesc( time_frames, [1:C], wav_mfcc(2:end,:) ); % HTK's TARGETKIND: MFCC
    %imagesc( time_frames, [1:C+1], wav_mfcc );       % HTK's TARGETKIND: MFCC_0
    axis( 'xy' );
    xlim( [ min(time_frames) max(time_frames) ] );
    xlabel( 'Time (s)' ); 
    ylabel( 'Cepstrum index' );
    title( 'Mel frequency cepstrum' );
    colorbar;

    subplot(2,1,2);
    plot(t,wav); xlabel('Time(s)'); ylabel('Amplitude'); title('Speech Wave');
    set(gca,'xlim',[0 t(length(t))]);    
    
end