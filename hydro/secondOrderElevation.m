% Second–Order Surface Elevation η^{(2)}
function eta2 = secondOrderElevation(varargin)
    % Function can be called in two ways:
    % 1. eta2 = secondOrderElevation(A, phi, omega, k, h, t)
    % 2. eta2 = secondOrderElevation(t, params)
    %
    % eta2(t) = sum_{n,m} A(n)*A(m)*[ L+_{nm} cos((ωn+ωm)t+φn+φm)
    %                               + L-_{nm} cos((ωn-ωm)t+φn-φm) ]
    % where L±_{nm} is given by Eq.3.26 in Kim et al. (2014).
    
    % Check input arguments
    if nargin == 2 && isstruct(varargin{2})
        % Called as secondOrderElevation(t, params)
        t = varargin{1};
        params = varargin{2};
        
        % Extract parameters from params struct
        A = params.A_m;
        phi = params.phiRand;
        omega = params.omegaVec;
        k = params.k_m_vec;
        h = params.h;
    else
        % Called as secondOrderElevation(A, phi, omega, k, h, t)
        A = varargin{1};
        phi = varargin{2};
        omega = varargin{3};
        k = varargin{4};
        h = varargin{5};
        t = varargin{6};
    end

    g    = 9.81;
    N    = numel(omega);
    eta2 = zeros(size(t));
    epsk = 1e-6;

    for n = 1:N
        Rn = k(n)*tanh(k(n)*h);               % Eq.3.28 :contentReference[oaicite:1]{index=1}
        for m = 1:N
            Rm = k(m)*tanh(k(m)*h);
            ksum  = k(n) + k(m);
            kdiff = abs(k(n) - k(m));

            % ——— sum-frequency transfer function L+_{nm} ———
            if ksum > epsk
                % numerator D+_{nm}, Eq.3.29 with “+”
                Dp = ( (sqrt(Rn)+sqrt(Rm)) * ...
                       ( sqrt(Rm)*(k(n)^2 - Rn^2) + sqrt(Rn)*(k(m)^2 - Rm^2) ) / ...
                       ( (sqrt(Rn)+sqrt(Rm))^2 - ksum*tanh(ksum*h) ) ) ...
                   + ( 2*(sqrt(Rn)+sqrt(Rm))^2 * ( k(n)*k(m) - Rn*Rm ) / ...
                       ( (sqrt(Rn)+sqrt(Rm))^2 - ksum*tanh(ksum*h) ) );
                % L+ per Eq.3.26
                Lp = 0.25 * ( Dp ...
                          - (k(n)*k(m) - Rn*Rm)/sqrt(Rn*Rm) ...
                          + (Rn + Rm) );
                wsum = omega(n) + omega(m);
                eta2 = eta2 + A(n)*A(m) * Lp .* cos(wsum*t + phi(n) + phi(m));
            end

            % ——— difference-frequency transfer function L-_{nm} ———
            if kdiff > epsk
                Dm = ( (sqrt(Rn)-sqrt(Rm)) * ...
                       ( sqrt(Rm)*(k(n)^2 - Rn^2) - sqrt(Rn)*(k(m)^2 - Rm^2) ) / ...
                       ( (sqrt(Rn)-sqrt(Rm))^2 - kdiff*tanh(kdiff*h) ) ) ...
                   + ( 2*(sqrt(Rn)-sqrt(Rm))^2 * ( k(n)*k(m) + Rn*Rm ) / ...
                       ( (sqrt(Rn)-sqrt(Rm))^2 - kdiff*tanh(kdiff*h) ) );
                Lm = 0.25 * ( Dm ...
                          - (k(n)*k(m) + Rn*Rm)/sqrt(Rn*Rm) ...
                          + (Rn + Rm) );
                wdiff = omega(n) - omega(m);
                eta2 = eta2 + A(n)*A(m) * Lm .* cos(wdiff*t + phi(n) - phi(m));
            end
        end
    end
end