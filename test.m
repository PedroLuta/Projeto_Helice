syms Cr Cm rm Ct a b c d

eqn1 = (d)^0.5 == Cr;
eqn2 = (a + b + c + d)^0.5 == Ct;
eqn3 = ((a*(rm^3)) + (b*(rm^2)) + (c*rm) + d)^0.5 == Cm;
eqn4 = ((3*a*(rm^2)) + (2*b*rm) + c)/(2*(((a*(rm^3)) + (b*(rm^2)) + (c*rm) + d)^0.5)) == 0;
solve(eqn1, eqn2, eqn3, eqn4, [a b c d])