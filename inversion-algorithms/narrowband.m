function I = narrowband(s,t,I1,fb, fa,n)

%INVERSÃO RECUSIRVA

% s ... traços sísmicos
% t ... vetor do tempo do traço sísmico
% znot ... primeiro valor da impedância
% fb ... frequência mais baixa no traço
% fa ... frequência mais alta no traço


%Dados sinteitcos
%[w, tw] = ricker(25, 0.003); %WAVELET TIPO RICKER
%plot(tw,w);
%aux=s_spectrum(w) %Espectro

%[s,r,t] = make_traces(200,200,0.00003,w,'gauss_mixture',[0.001, 0.1, 0.5],0.85);

%fa=120; 
%[refl,dp] = sparse_decon(s,w,0.85,20); %Desconvolução
%S = conv2(refl,w); % Modelo Convolucional 

%%

S=padpow2(s);
T=xcoord(t(1),t(2)-t(1),S);

%computar o spectro
[G,f] = fftrl(S,T); %transformada de Fourier
                    % G espectro 
                    % f frequencia

%Prever as baixas frequencias com um filtro
%T = predlowf(G,f,fb,fa,n);

rc = ifftrl(G,f);
if(length(rc)~=length(s))
    %rc=rc(1:length(s));
end

%irc=cumprod(rc); 
%I = I1*irc;
%ou
ircs=cumsum(rc);
I=I1*exp(2*ircs);


%%

