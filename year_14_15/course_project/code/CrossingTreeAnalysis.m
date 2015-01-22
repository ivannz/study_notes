clear all; 

% Data
N = 2 ^ 18 ; % sample size
t = [ 0 : N - 1 ] ./ ( N - 1 ) ; % Time scale
MC = 100 ; % # of monte carlo realizations
H = 0.95 ; %0.5 for Brownian motion

PSYNTH = '[ data, fgn ] = synthfbmcircul( N + 1, H ) ; data = data( 2 : end ) ;'; 
eval( PSYNTH ) ; 

% Offspring distribution: we estimate P(Z = 2), P(Z = 4),..., P(Z = 12)
Z1 = 2 ; Z2 = 40 ; Zr = Z1 : 2 : Z2 ;
PrZ = zeros( MC, ( Z2 - Z1 ) / 2 + 1 ) ;
for imc = 1 : MC
    if rem( imc, round( MC / 10 ) ) == 0 ; fprintf( '%2g ', imc ) ; end
    eval( PSYNTH ) ;
    dataR = max( data ) - min( data ) ;
    j1 = floor( log2( dataR ) ) - 6 ;
    j2 = floor( log2( dataR ) ) ;
    % Crossing tree analysis: w = crossing durations
    [ w, subx, hp, ht ] = f_get_w( data, t, [ j1 : j2 ], 1, 0 ) ;
    Z = [ subx{ 3 } subx{ 4 } ] ;
    for z = Zr
        PrZ( imc, z / 2 ) = sum( Z == z ) / length( Z ) ;
    end
end

% Estimated offspring distribution and std error
mPrZ = mean( PrZ )
sPrZ = std( PrZ )

% hist( PrZ, 100 )

figure
plot( PrZ )

figure
plot( log( mPrZ( mPrZ > 0 ) ) )

%% Results for H=.55
%% 
%% Looks exponential for H=.65, but has a downward kink
%%  for larger-band crossings, hinting, possibly, at a
%%  distribution with compact support.
%% Resutls for H=.75
%% 
%% Let's see the MC results for H=.85
%% 
%% And, at last, for H=.95
%% 

%% Strange, there are cases when the probability of a
%%  narrower crossing is lower than that of a wider one.



figure
hold on
    plot( t, data, '-k', 'LineWidth', 1.5 )
    for i = 1:length( subx ) ; stairs( ht{i}, hp{i}, 'r' ) ; end
hold off











