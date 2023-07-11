%First message modulation
[r1,fs1] = audioread('Short_FM9090.wav');  %Read the first message
channel1_1 = r1(:,1);
channel2_1 = r1(:,2);
m1 = channel1_1+channel2_1;         %adding the two channels
message1 = interp(m1,16);           %changing the sampling rate
Fs = fs1*16;                              
T = 1/Fs;         
L = length(message1);           
t = (0:L-1)*T;       
carrier1 = cos(2*pi*1e5*t);
carrier1=carrier1.';                %Transpose the carrier
s1 = message1.*carrier1;            %Multipling the message and the carrier

%Second message modulation
[r2,fs2] = audioread('Short_SkyNewsArabia.wav');  %Read the second message
channel1_2 = r2(:,1);
channel2_2 = r2(:,2);
m2 = channel1_2+channel2_2;         %Adding the two channels
message2 = interp(m2,16);           %Changing the sampling rate                                 
L2 = length(message2);        
t2 = (0:L2-1)*T;       
carrier2 = cos(2*pi*1.5e5*t2);
carrier2=carrier2.';                %Transpose the carrier
s2 = message2.*carrier2;            %Multipling the message and the carrier

s1(numel(s2)) = 0;                  %pad the short signal with zeros
sent = s1 + s2;                     %Adding the two messages

%Plotting the transmitted message
Lsent = length(sent);
Ysent = fftshift(fft(sent));
fsent = Fs*(-(Lsent/2):(Lsent/2)-1)/Lsent;
figure('Name','Transmitted message') 
plot(fsent,abs(Ysent))
xlabel("F (Hz)")
ylabel("Amplitude")


%RF Stage
n=input("Choose the message (0 or 1): ");     %To choose the first or the second message
A_stop1 = 60;		% Attenuation in the first stopband
A_pass = 1;     	% Amount of ripple allowed in the passband
F_stop1 = 0.6e5+50e3*n;	% Edge of the stopband 
F_pass1 = 0.8e5+50e3*n;	% Edge of the passband 
F_pass2 = 1.15e5+50e3*n;	% Closing edge of the passband 
F_stop2 = 1.2e5+50e3*n;	% Edge of the second stopband
A_stop2 =60;		% Attenuation in the second stopband


RFBandPassSpecObj = ...
   fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
		F_stop1, F_pass1, F_pass2, F_stop2, A_stop1, A_pass, ...
		A_stop2, Fs);

RFFilter = design(RFBandPassSpecObj, 'equiripple');
RFoutput = filter(RFFilter,sent);

%Plotting RF stage output
LRFoutput =length(RFoutput);
YRFoutput = fftshift(fft(RFoutput));
fRFoutput = Fs*(-(LRFoutput/2):(LRFoutput/2)-1)/LRFoutput;
figure('Name','RF Stage output') 
plot(fRFoutput,abs(YRFoutput))
xlabel("F (Hz)")
ylabel("Amplitude")


%The Oscillator 
t = (0:LRFoutput-1)*T;  
mixer1 = cos(2*pi*(1.25+n*0.5)*1e5*t);
mixer1=mixer1.';
mixer1out = RFoutput.*mixer1;

%Plotting the Oscillator output
Ymixer1 = fftshift(fft(mixer1out));
Lmixer1 =length(mixer1out);
fmixer1 = Fs*(-(Lmixer1/2):(Lmixer1/2)-1)/Lmixer1;
figure('Name','The Oscillator output') 
plot(fmixer1,abs(Ymixer1))
xlabel("F (Hz)")
ylabel("Amplitude")

%IF Stage
A_stop1 = 60;		% Attenuation in the first stopband
A_pass = 1;     	% Amount of ripple allowed in the passband
F_stop1 = 0.5e3;	% Edge of the stopband
F_pass1 = 1e3;	% Edge of the passband
F_pass2 = 50e3;	% Closing edge of the passband
F_stop2 = 60e3;	% Edge of the second stopband
A_stop2 =60;		% Attenuation in the second stopband


IFBandPassSpecObj = ...
   fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
		F_stop1, F_pass1, F_pass2, F_stop2, A_stop1, A_pass, ...
		A_stop2, Fs);

IFFilter = design(IFBandPassSpecObj, 'equiripple');
IFoutput = filter(IFFilter,mixer1out);

%Plotting IF stage output
LIFoutput =length(IFoutput);
YIFoutput = fftshift(fft(IFoutput));
fIFoutput = Fs*(-(LIFoutput/2):(LIFoutput/2)-1)/LIFoutput;
figure('Name','IF stage output') 
plot(fIFoutput,abs(YIFoutput))
xlabel("F (Hz)")
ylabel("Amplitude")


%The Baseband detection
%The second Mixer
t = (0:LIFoutput-1)*T;  
mixer2 = cos(2*pi*25e3*t);
mixer2=mixer2.';
mixer2out = IFoutput.*mixer2;

%Plotting The Baseband detection mixer output
Lmixer2 =length(mixer2out);
Ymixer2 = fftshift(fft(mixer2out));
fmixer2 = Fs*(-(Lmixer2/2):(Lmixer2/2)-1)/Lmixer2;
figure('Name','The Baseband detection mixer output') 
plot(fmixer2,abs(Ymixer2))
xlabel("F (Hz)")
ylabel("Amplitude")

%Lowpass Filter
LowPassSpecObj=fdesign.lowpass('Fp,Fst,Ap,Ast',25e3,30e3,1,60,Fs);
lpFilter = design(LowPassSpecObj,'equiripple');
Receivedsignal = filter(lpFilter,mixer2out);

%Plotting The Received signal
LReceived =length(Receivedsignal);
YReceived = fftshift(fft(Receivedsignal));
fReceived = Fs*(-(LReceived/2):(LReceived/2)-1)/LReceived;
figure('Name','The Received signal') 
plot(fReceived,abs(YReceived))
xlabel("F (Hz)")
ylabel("Amplitude")