function [w, subx, hit_point, hit_time] = f_get_w( x, t, levels, delta, deleteFirst )

% times, points and subcrossings for continuous process
%
% x: signal to be analysed
% t: time ( mandatory )
% levels: for each level in levels crossings will be calculated at scale delta * 2^level
% delta: if not 0 then is used as base scale, o/w use stddev( diff( x ) )
% deleteFirst: if 1 then delete first crossing at each level
% 
% w: crossing duration ( cell structure, size of the number of levels )
% subx : number of subcrossings ( cell structure, size of the number of levels, first empty )
% hit_point: cell struct., hitting levels 
% hit_time: cell struct., hitting times
% Y.Shen 12.2002, O.D.Jones 7.2003, 6.2007
% P - O. Amblard & G. Decrouez 2011
%
%  [w, subx, hit_point, hit_time] = f_get_w( fname, levels, delta, deleteFirst )


% hitting times and points
lx = length( x );   
if delta == 0
%% why is std( diff ) better than min( abs( diff ) )
    % delta = std( diff( x ) );
    delta = min( abs( diff ) ) * 2 ;
end
%% Define lists of consequitive hitting times and levels
hit_point = cell( length( levels ), 1 );
hit_time = cell( length( levels ), 1 );

w = hit_time;
comptscale = 1;

%% Precomputations
x_0 = ( x - x( 1 ) ) / delta;
resolution = max( diff( x_0 ) );

for n = levels
    scale = 2^n;

%% Define hitting times and points
    % h_t = zeros( 1, lx * ceil( max( diff( y ) ) ) ); %% at this scale, upper bound on the total # of crossings
    h_t = zeros( 1, lx * ceil( resolution / scale ) ); %% at this scale, upper bound on the total # of crossings
    h_p = h_t;

    if deleteFirst
        last_hit = 0;
    else
        last_hit = 1;
    end

%% Might be slow due to O(N) rescaling
    y = x_0 / scale;

    compt = 1;
    for i = 2:lx
%% Won't this be true alomst every time due to numerical issues?
        if y( i - 1 ) ~= y( i )
%% NOTICE! Conditions do not depend on the current scale
            if y( i - 1 ) < y( i )
                x_init = ceil( y( i - 1 ) );
                x_final = floor( y( i ) );
                step = 1;
            else
                x_init = floor( y( i - 1 ) );
                x_final = ceil( y( i ) );
                step = - 1;
            end
%% x_init and x_final represent the first and the last normalised levels, which
%%   the sample path has crossed. Obviously, if the 

            for j = x_init : step : x_final
%% It is not a crossing if the process returned to the level of the previous crossing
%%  T^n_{k+1} \defn \inf\{  . t > T^n_k : \, X_t\in \partial B_n,\, X_t \neq X_{T^n_k} \}
                if j ~= last_hit
%% Due to discrete data the the hitting times are approximated using linear interpolation
                    h_t( compt ) = t( i - 1 ) + ( j - y( i - 1 ) ) / ( y( i ) - y( i - 1 ) ) * ( t( i ) - t( i - 1 ) );

                    h_p( compt ) = j * delta * scale + x( 1 );
                    compt = compt + 1;             
                    last_hit = j;
                end
            end
        end
    end
    hit_point{ comptscale } = h_p( 1:compt - 1 );
    hit_time{ comptscale } = h_t( 1:compt - 1 );
    w{ comptscale } = diff( h_t( 1:compt - 1 ) );
    comptscale = comptscale + 1;   
end   

% subcrossing numbers     
hit0 = [hit_time{ 1 }' hit_point{ 1 }'];
subx = cell( length( levels ), 1 );

for level = 2:length( levels )   
    hit1 = [hit_time{ level }' hit_point{ level }'];  
    if ~isempty( hit1 )
        j0 = 1;
        
        sx = zeros( 1, size( hit1, 1 ) - 1 );
        compt = 1;
        while hit0( j0, 1 ) ~= hit1( 1, 1 ), j0 = j0 + 1; end
        for i = 2:size( hit1, 1 )
            j1 = j0 + 1;
            while hit0( j1, 1 ) ~= hit1( i, 1 ), j1 = j1 + 1; end           
            sx( compt ) = j1 - j0;
            compt = compt + 1;
            j0 = j1;
        end
	end
    subx{ level } = sx;
    hit0 = hit1;
end