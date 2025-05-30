function [u2x,u2y,u2z] = secondOrderVelocity(A,phi,omega,k,theta,h,x,zVec,t)
    % second-order velocity over all depths zVec and times t
    N     = numel(omega);
    Nz    = numel(zVec);
    Nt    = numel(t);
    u2x   = zeros(Nz, Nt);
    u2y   = zeros(Nz, Nt);
    u2z   = zeros(Nz, Nt);
    epsw  = 1e-6;

    for iz = 1:Nz
        z = zVec(iz);
        tmpX = zeros(1, Nt);
        tmpY = zeros(1, Nt);
        tmpZ = zeros(1, Nt);

        for n = 1:N
            for m = 1:N
                % — sum-frequency term
                Bp   = computeBplus(n,m,omega,k,theta,h);
                if abs(Bp) > epsw
                    wsum = omega(n) + omega(m);
                    xUp  = Bp*(k(n)*cos(theta(n)) + k(m)*cos(theta(m)));
                    yUp  = Bp*(k(n)*sin(theta(n)) + k(m)*sin(theta(m)));
                    ksum = abs(k(n)*exp(1i*theta(n)) + k(m)*exp(1i*theta(m)));
                    zUp  = 1i*Bp*ksum.*tanh(ksum*(h+z));
                    phaseP = cos(wsum*t + phi(n) + phi(m));
                    tmpX = tmpX + A(n)*A(m) * xUp  .* phaseP;
                    tmpY = tmpY + A(n)*A(m) * yUp  .* phaseP;
                    tmpZ = tmpZ + A(n)*A(m) * real(zUp .* phaseP);
                end

                % — diff-frequency term
                Bm   = computeBminus(n,m,omega,k,theta,h);
                if abs(Bm) > epsw
                    wdiff = omega(n) - omega(m);
                    xUm   = Bm*(k(n)*cos(theta(n)) - k(m)*cos(theta(m)));
                    yUm   = Bm*(k(n)*sin(theta(n)) - k(m)*sin(theta(m)));
                    kdiff = abs(k(n)*exp(1i*theta(n)) - k(m)*exp(1i*theta(m)));
                    zUm   = 1i*Bm*kdiff.*tanh(kdiff*(h+z));
                    phaseM = cos(wdiff*t + phi(n) - phi(m));
                    tmpX = tmpX + A(n)*A(m) * xUm  .* phaseM;
                    tmpY = tmpY + A(n)*A(m) * yUm  .* phaseM;
                    tmpZ = tmpZ + A(n)*A(m) * real(zUm .* phaseM);
                end
            end
        end

        u2x(iz,:) = tmpX;
        u2y(iz,:) = tmpY;
        u2z(iz,:) = tmpZ;
    end
end
