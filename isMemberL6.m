function b = isMemberL6(i,e,p)

    b = p(1,e(1,i)) < 0.001*10^(-3) && p(1,e(2,i)) < 0.001*10^(-3);

end