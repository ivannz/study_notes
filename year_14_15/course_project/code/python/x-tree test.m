clear all; format compact;

N = 2^15; % sample size

H = 0.7;

choix='b';

PSYNTH='[data,t] = weir(N,1.2,1000,H,choix);';eval(PSYNTH);

t = [0:N-1]./(N-1);

delta = std( diff( data ));
Z = (data - data(1))/delta;
% plot( t, Z )

% [ w, subx, hp, ht ] = f_get_w( Z, t, [ 0:0 ], 1, 0 ) ;

lx = length(Z);
h_t=zeros(1,lx*ceil(max(diff(Z)))); %% at this scale, upper bound on the total # of crossings
h_p=h_t;
compt=1;
last_hit = 0;
fZ = floor( Z );
for i = 1:(lx-1)
    if Z(i) ~= Z(i+1)
        if Z(i) < Z(i+1)
            x_init = fZ(i)+1 ;%ceil(Z(i));
            x_final = fZ(i+1) ;%floor(Z(i+1));
            step = 1;
        else
            x_init = fZ(i) ;%floor(Z(i));
            x_final = fZ(i+1)+1 ;%ceil(Z(i+1));
            step = -1;
        end
        for j = x_init:step:x_final
            if j ~= last_hit
                h_t(compt) = t(i) + (j - Z(i))/(Z(i+1) - Z(i))*(t(i+1) - t(i));
                h_p(compt) = j;
                compt=compt+1;             
                last_hit = j;
            end
        end
    end
end
h_p=h_p(1:compt-1);
h_t=h_t(1:compt-1);
