function b = isMemberL5(i,e,p)

    b = (e(6,i) == 2 && e(7,i) == 0 && p(2,e(1,i)) < 0.001*10^(-3) && p(2,e(2,i)) < 0.001*10^(-3));

end