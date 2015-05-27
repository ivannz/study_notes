clear all ; format compact ;

basedir = tempname
mkdir( basedir ) ;

%%%%% HRM - 2
Hfgn = 0.8 ;    % H=0.6
N = 2^20 ; % sample size
downsample = 16 ;
sigma2 = 1 ;
H = 1+2*(Hfgn-1) ;
for i = 1:1000
      [B,xn] = synthfbmcircul2( N * downsample + 1, Hfgn, sigma2 ) ;
      data = cumsum( xn.^2 - 1 ) ; data = data( 2 : downsample : end ) ;
      data = data / max( data );
      filename = [basedir, '\', 'HRM-2_20-16_', num2str( H, 4 ), '_', num2str( i ), '.mat'] ;
      save( filename, 'data' );
end

Hfgn = 0.85 ;   % H=0.7 
N = 2^20 ; % sample size
downsample = 16 ;
sigma2 = 1 ;
H = 1+2*(Hfgn-1) ;
for i = 1:1000
      [B,xn] = synthfbmcircul2( N * downsample + 1, Hfgn, sigma2 ) ;
      data = cumsum( xn.^2 - 1 ) ; data = data( 2 : downsample : end ) ;
      data = data / max( data );
      filename = [basedir, '\', 'HRM-2_20-16_', num2str( H, 4 ), '_', num2str( i ), '.mat'] ;
      save( filename, 'data' );
end

Hfgn = 0.9 ;    % H=0.8 
N = 2^20 ; % sample size
downsample = 16 ;
sigma2 = 1 ;
H = 1+2*(Hfgn-1) ;
for i = 1:1000
      [B,xn] = synthfbmcircul2( N * downsample + 1, Hfgn, sigma2 ) ;
      data = cumsum( xn.^2 - 1 ) ; data = data( 2 : downsample : end ) ;
      data = data / max( data );
      filename = [basedir, '\', 'HRM-2_20-16_', num2str( H, 4 ), '_', num2str( i ), '.mat'] ;
      save( filename, 'data' );
end

Hfgn = 0.95 ;    % H=0.9
N = 2^20 ; % sample size
downsample = 16 ;
sigma2 = 1 ;
H = 1+2*(Hfgn-1) ;
for i = 1:1000
      [B,xn] = synthfbmcircul2( N * downsample + 1, Hfgn, sigma2 ) ;
      data = cumsum( xn.^2 - 1 ) ; data = data( 2 : downsample : end ) ;
      data = data / max( data );
      filename = [basedir, '\', 'HRM-2_20-16_', num2str( H, 4 ), '_', num2str( i ), '.mat'] ;
      save( filename, 'data' );
end

%%%%% HRM - 3
Hfgn = 0.9-1/30 ;    % H=0.6
N = 2^20 ; % sample size
downsample = 16 ;
sigma2 = 1 ;
H = 1+3*(Hfgn-1) ;
for i = 1:1000
      [B,xn] = synthfbmcircul2( N * downsample + 1, Hfgn, sigma2 ) ;
      data = cumsum(xn.^3 - 3*xn) ; data = data( 2 : downsample : end ) ;
      data = data / max( data );
      filename = [basedir, '\', 'HRM-3_20-16_', num2str( H, 4 ), '_', num2str( i ), '.mat'] ;
      save( filename, 'data' );
end

Hfgn = 0.9 ;   % H=0.7 
N = 2^20 ; % sample size
downsample = 16 ;
sigma2 = 1 ;
H = 1+3*(Hfgn-1) ;
for i = 1:1000
      [B,xn] = synthfbmcircul2( N * downsample + 1, Hfgn, sigma2 ) ;
      data = cumsum(xn.^3 - 3*xn) ; data = data( 2 : downsample : end ) ;
      data = data / max( data );
      filename = [basedir, '\', 'HRM-3_20-16_', num2str( H, 4 ), '_', num2str( i ), '.mat'] ;
      save( filename, 'data' );
end

Hfgn = 0.9+1/30 ;    % H=0.8 
N = 2^20 ; % sample size
downsample = 16 ;
sigma2 = 1 ;
H = 1+3*(Hfgn-1) ;
for i = 1:1000
      [B,xn] = synthfbmcircul2( N * downsample + 1, Hfgn, sigma2 ) ;
      data = cumsum(xn.^3 - 3*xn) ; data = data( 2 : downsample : end ) ;
      data = data / max( data );
      filename = [basedir, '\', 'HRM-3_20-16_', num2str( H, 4 ), '_', num2str( i ), '.mat'] ;
      save( filename, 'data' );
end

Hfgn = 0.9+2/30 ;    % H=0.9
N = 2^20 ; % sample size
downsample = 16 ;
sigma2 = 1 ;
H = 1+3*(Hfgn-1) ;
for i = 1:1000
      [B,xn] = synthfbmcircul2( N * downsample + 1, Hfgn, sigma2 ) ;
      data = cumsum(xn.^3 - 3*xn) ; data = data( 2 : downsample : end ) ;
      data = data / max( data );
      filename = [basedir, '\', 'HRM-3_20-16_', num2str( H, 4 ), '_', num2str( i ), '.mat'] ;
      save( filename, 'data' );
end

%%%%% HRM - 4
Hfgn = 0.90 ;    % H=0.6
N = 2^20 ; % sample size
downsample = 16 ;
sigma2 = 1 ;
H = 1+4*(Hfgn-1) ;
for i = 1:1000
      [B,xn] = synthfbmcircul2( N * downsample + 1, Hfgn, sigma2 ) ;
      data = cumsum(xn.^4 - 6*xn.^2 + 3) ; data = data( 2 : downsample : end ) ;
      data = data / max( data );
      filename = [basedir, '\', 'HRM-4_20-16_', num2str( H, 4 ), '_', num2str( i ), '.mat'] ;
      save( filename, 'data' );
end

Hfgn = 0.925 ;   % H=0.7 
N = 2^20 ; % sample size
downsample = 16 ;
sigma2 = 1 ;
H = 1+4*(Hfgn-1) ;
for i = 1:1000
      [B,xn] = synthfbmcircul2( N * downsample + 1, Hfgn, sigma2 ) ;
      data = cumsum(xn.^4 - 6*xn.^2 + 3) ; data = data( 2 : downsample : end ) ;
      data = data / max( data );
      filename = [basedir, '\', 'HRM-4_20-16_', num2str( H, 4 ), '_', num2str( i ), '.mat'] ;
      save( filename, 'data' );
end

Hfgn = 0.95 ;    % H=0.8 
N = 2^20 ; % sample size
downsample = 16 ;
sigma2 = 1 ;
H = 1+4*(Hfgn-1) ;
for i = 1:1000
      [B,xn] = synthfbmcircul2( N * downsample + 1, Hfgn, sigma2 ) ;
      data = cumsum(xn.^4 - 6*xn.^2 + 3) ; data = data( 2 : downsample : end ) ;
      data = data / max( data );
      filename = [basedir, '\', 'HRM-4_20-16_', num2str( H, 4 ), '_', num2str( i ), '.mat'] ;
      save( filename, 'data' );
end

Hfgn = 0.975 ;    % H=0.9
N = 2^20 ; % sample size
downsample = 16 ;
sigma2 = 1 ;
H = 1+4*(Hfgn-1) ;
for i = 1:1000
      [B,xn] = synthfbmcircul2( N * downsample + 1, Hfgn, sigma2 ) ;
      data = cumsum(xn.^4 - 6*xn.^2 + 3) ; data = data( 2 : downsample : end ) ;
      data = data / max( data );
      filename = [basedir, '\', 'HRM-4_20-16_', num2str( H, 4 ), '_', num2str( i ), '.mat'] ;
      save( filename, 'data' );
end




basedir
