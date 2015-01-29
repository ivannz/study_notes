clear all; format compact;

% Data

N=2^15; % sample size
PROC=[3]; % select process3

    switch PROC;
        case 1 % FBM
            H = 0.6 ;
            PSYNTH='[data,fgn] =synthfbmcircul(N+1,H) ; data=data(2:end);'; eval(PSYNTH);
            t= [0:N-1]./(N-1);
            ynamep = ['fBm ',num2str(H)]         
        case 21 % Rosenblatt (Hermite ordre 2)
            %Hfgn = 0.8 ;    % H=0.6 
            %Hfgn = 0.85 ;   % H=0.7 
            Hfgn = 0.9 ;    % H=0.8 
            %Hfgn = 0.95 ;    % H=0.9
            H=1+2*(Hfgn-1)
            sigma2=1;
            % increasing downsample = better approximation
            %PSYNTH='downsample = 128 ;[B,xn] = synthfbmcircul2(N*downsample+1,Hfgn,sigma2) ; data = cumsum(xn.^2 - 1) ; data = data(2:downsample:end) ; data=data/max(data);';
            PSYNTH='downsample = 16 ;[B,xn] = synthfbmcircul2(N*downsample+1,Hfgn,sigma2) ; data = cumsum(xn.^2-1) ; data = data(2:downsample:end) ; data=data/max(data);';  eval(PSYNTH);
            t= [0:N-1]./(N-1);
            ynamep = ['Rosenblatt ',num2str(H)]
        case 22 % Hermite ordre 3
            %Hfgn = 0.9-1/30 ;    % H=0.6 
            Hfgn = 0.9 ;         % H=0.7 
            %Hfgn = 0.9+1/30 ;    % H=0.8 
            %Hfgn = 0.9+2/30 ;    % H=0.9
            H=1+3*(Hfgn-1)
            sigma2=1;
            % increasing downsample = better approximation
            %PSYNTH='downsample = 128 ;[B,xn] = synthfbmcircul2(N*downsample+1,Hfgn,sigma2) ; data = cumsum(xn.^3 - 3*xn) ; data = data(2:downsample:end) ; data=data/max(data);';
            PSYNTH='downsample = 16 ;[B,xn] = synthfbmcircul2(N*downsample+1,Hfgn,sigma2) ; data = cumsum(xn.^3 - 3*xn) ; data = data(2:downsample:end) ; data=data/max(data);';    eval(PSYNTH);
            t= [0:N-1]./(N-1);
            ynamep = ['Hermite 3 ',num2str(H)]
            
        case 23 % Hermite ordre 4
            %Hfgn = 0.90 ;    % H=0.6 
            Hfgn = 0.925 ;   % H=0.7 
            %Hfgn = 0.95 ;    % H=0.8 
            %Hfgn = 0.975 ;    % H=0.9
            H=1+4*(Hfgn-1)
            sigma2=1;
            % increasing downsample = better approximation
            %PSYNTH='downsample = 128 ;[B,xn] = synthfbmcircul2(N*downsample+1,Hfgn,sigma2) ; data = cumsum(xn.^4 - 6*xn.^2 + 3) ; data = data(2:downsample:end) ; data=data/max(data);';
            PSYNTH='downsample = 16 ;[B,xn] = synthfbmcircul2(N*downsample+1,Hfgn,sigma2) ; data = cumsum(xn.^4 - 6*xn.^2 + 3) ; data = data(2:downsample:end) ; data=data/max(data);';    eval(PSYNTH);
            t= [0:N-1]./(N-1);
            ynamep = ['Hermite 4 ',num2str(H)]
            
        case 3 % Weierstrass
            H=0.7;
            choix='a';
            PSYNTH='[data,t] =weir(N,1.2,1000,H,choix);';eval(PSYNTH);
            t= [0:N-1]./(N-1);
            ynamep = ['Weierstrass H=',num2str(H)]    

    end
    
    plot(t,data)