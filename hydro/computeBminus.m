function Bm = computeBminus(n,m,omega,k,theta,h)
    % compute the “–” quadratic transfer function B⁻_{nm} (Longuet–Higgins, Eq.4.38)
    g      = 9.81;
    wdiff  = omega(n) - omega(m);
    kn     = k(n);    km    = k(m);
    thn    = theta(n); thm   = theta(m);
    kn_km  = kn*km*cos(thn-thm);
    num    = wdiff * ( kn_km + kn*km );
    denom  = wdiff^2  - g*abs(kn - km);
    Bm     = (g^2/(omega(n)*omega(m))) * (num/denom);
end
