function Bp = computeBplus(n,m,omega,k,theta,h)
    % compute the “+” quadratic transfer function B⁺_{nm} (Longuet–Higgins, Eq.4.38)
    g     = 9.81;
    wsum  = omega(n) + omega(m);
    kn    = k(n);    km   = k(m);
    thn   = theta(n); thm  = theta(m);
    % dot product kn·km = kn*km*cos(θn−θm)
    kn_km = kn*km*cos(thn-thm);
    num   = wsum * ( kn_km - kn*km );
    denom = wsum^2  - g*(kn + km);
    Bp    = (g^2/(omega(n)*omega(m))) * (num/denom);
end
