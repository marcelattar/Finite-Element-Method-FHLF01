function b =isMemberL1(i,e,p)

    b = e(6,i)==3 && e(7,i)==0 && p(2,e(1,i)) > 0.599*10^(-3) && p(2,e(2,i)) > 0.599*10^(-3) && p(1,e(1,i)) < 0.21*10^(-3) && p(1,e(2,i)) < 0.21*10^(-3);
    % multiplies with 10^(-3) to convert to meter.
end