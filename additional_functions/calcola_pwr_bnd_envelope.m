% calcola PSD  del segnale audio

function   [pwr_bnd] = calcola_pwr_bnd_envelope(S,Fs,bande,TR)

tc = 1 /Fs;
time = 0 : tc : (length(S)-1)*tc;
n_bin = ceil((length(time)-TR*Fs)/(TR*Fs));
n_bnd=size(bande,2);
pwr_bnd = zeros(n_bnd, n_bin);

%% ------ spectral power density of the signal

% ------ definizione finestra -----------------
N_win = Fs;    
SD = Fs/6; 
alpha = (N_win - 1)/(2*SD);
window = gausswin(N_win,alpha);
% ---------------------------------------------

% concateniamo il segnale in blocchi per evitare un for col pwelch
S_final=zeros(TR*Fs,n_bin);
for tt = 1 : n_bin
    start = 1+(tt-1)*TR*Fs;
    stop = (start-1) + TR*Fs;
    S_final(:,tt)=S(start:stop);
end

disp(sprintf('Starting computations...it may take a while...'));
 tic
 [Pxx,freq] = pwelch(S_final,window,[],int32(2*Fs),int32(Fs)); % overlap 75%
 toc
 all_Pxx=Pxx';

% tic
% all_Pxx=zeros(n_bin,Fs+1);
% for tt = 1 : n_bin
%     start = 1+(tt-1)*TR*Fs;
%     stop = (start-1) + TR*Fs;
% 
%     [Pxx,freq] = pwelch(S(start:stop),window,int32(N_win-0.001*Fs),int32(2*Fs),int32(Fs)); % overlap 75%
%     all_Pxx(tt,:)=Pxx;
% 
%     if mod(tt,100)==0
%         disp(sprintf('processing TR %d', tt));
%         toc
%         tic
%     end
% end
% clear start stop Pxx
% 

ind_bande= zeros(size(bande));
for k = 1 : size(bande,2)
    ind_bande(1,k) = find(freq >= bande(1,k),1,'first');
    ind_bande(2,k) = find(freq <= bande(2,k),1,'last');
end

for i = 1 : size(bande,2)
    f1 = ind_bande(1,i);
    f2 = ind_bande(2,i);
    for tt = 1 :  n_bin
        pwr_bnd(length(bande)-i+1,tt) = trapz(freq(f1:f2),abs(all_Pxx(tt,f1:f2)));
    end
    clear f1 f2
end

% pwr_bnd = pow2db(pwr_bnd);



