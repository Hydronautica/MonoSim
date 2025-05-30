function [a2x,a2y,a2z] = secondOrderAcceleration(A,phi,omega,k,theta,h,x,zVec,t)
    % second-order acceleration over all depths zVec and times t
    N     = numel(omega);
    Nz    = numel(zVec);
    Nt    = numel(t);
    a2x   = zeros(Nz, Nt);
    a2y   = zeros(Nz, Nt);
    a2z   = zeros(Nz, Nt);
    epsw  = 1e-6;

    for iz = 1:Nz
        z = zVec(iz);
        tmpX = zeros(1, Nt);
        tmpY = zeros(1, Nt);
        tmpZ = zeros(1, Nt);

        for n = 1:N
            for m = 1:N
                % — sum-frequency term
                wsum = omega(n) + omega(m);
                if abs(wsum) > epsw
                    Bp   = computeBplus(n,m,omega,k,theta,h);
                    xUp  = Bp*(k(n)*cos(theta(n)) + k(m)*cos(theta(m)));
                    yUp  = Bp*(k(n)*sin(theta(n)) + k(m)*sin(theta(m)));
                    ksum = abs(k(n)*exp(1i*theta(n)) + k(m)*exp(1i*theta(m)));
                    zUp  = 1i*Bp*ksum.*tanh(ksum*(h+z));
                    phaseP = sin(wsum*t + phi(n) + phi(m));
                    tmpX = tmpX - A(n)*A(m) * xUp  .* (wsum.*phaseP);
                    tmpY = tmpY - A(n)*A(m) * yUp  .* (wsum.*phaseP);
                    tmpZ = tmpZ - A(n)*A(m) * real(zUp .* (wsum.*phaseP));
                end

                % — diff-frequency term
                wdiff = omega(n) - omega(m);
                if abs(wdiff) > epsw
                    Bm    = computeBminus(n,m,omega,k,theta,h);
                    xUm   = Bm*(k(n)*cos(theta(n)) - k(m)*cos(theta(m)));
                    yUm   = Bm*(k(n)*sin(theta(n)) - k(m)*sin(theta(m)));
                    kdiff = abs(k(n)*exp(1i*theta(n)) - k(m)*exp(1i*theta(m)));
                    zUm   = 1i*Bm*kdiff.*tanh(kdiff*(h+z));
                    phaseM = sin(wdiff*t + phi(n) - phi(m));
                    tmpX = tmpX - A(n)*A(m) * xUm  .* (wdiff.*phaseM);
                    tmpY = tmpY - A(n)*A(m) * yUm  .* (wdiff.*phaseM);
                    tmpZ = tmpZ - A(n)*A(m) * real(zUm .* (wdiff.*phaseM));
                end
            end
        end

        a2x(iz,:) = tmpX;
        a2y(iz,:) = tmpY;
        a2z(iz,:) = tmpZ;
    end
end
