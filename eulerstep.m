function anew = eulerstep(aold, f, C, K, N, tf)
    h = tf/N;           % Stepsize
    aOLD = aold;
    anew = zeros([size(aold) 1]);
    for i=1:N
        anew = (K+C./h)\(f+C*aOLD./h);
        aOLD = anew;
    end
end