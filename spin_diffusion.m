function signals = spin_diffusion(m, length_target, rstep, tmax, tstep, ...
    T1out, T1in, epsilon0, Din, Dout, length_source, slope, theta, ...
    dep, concjuice, concX, Cp, radii)
    
    r = linspace(0, length_target+length_source, rstep);
    t = linspace(0, tmax, tstep);
    
    
    % Extract the first solution of the pdepe equation, as p and pref.
    sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,r,t);
    p = sol(:,:,1);
    sol2 = pdepe(m,@pdex1pdenoDNP,@pdex1icnoDNP,@pdex1bcnoDNP,r,t); %solution of the diffusion equation without microwaves
    pref = sol2(:,:,1); %pref is the polarization without microwaves
    
    %% --- Visible polarization ---
    p_cut = zeros(length(t),length(r));
    pref_cut = zeros(length(t),length(r));
    for l = 1:length(r)
        p_cut(:,l) = p(:,l)*theta0(r(l)); %Visible polarization ON
        pref_cut(:,l) = pref(:,l)*theta0(r(l)); %Visible polarization
    end

    n_shells = length(radii) - 1;
    % extracted signals for each shell
    signals = struct('S_off', {{}}, 'S_on', {{}}, 't', {{}});
 
    %% --- Loop over shells ---
    for i = 1:n_shells
        r_min = radii(i);
        r_max = radii(i + 1);

        i_min = find(r >= r_min, 1, 'first');
        i_max = find(r <= r_max, 1, 'last');

        %Initialization of matrixes
        Presign = zeros(length(t), length(r));
        Presignref = zeros(length(t), length(r));
        SignalPol = zeros(1, length(t));
        SignalRef = zeros(1, length(t));
        NormSignalPol = zeros(1, length(t));
        NormSignalRef = zeros(1, length(t));

        %% -- Pre-signal --
        for j = 1:length(t)
            for k = i_min:i_max % integration only on the right portion
                Presign(j,k) = fc(r(k)) * p_cut(j,k) * r(k)^m; 
                Presignref(j,k) = fc(r(k)) * pref_cut(j,k) * r(k)^m;
            end
        end

        
        %% Calculation as a function of time
        % Calcul of the Signals by integrating the presignals
 
        for w = 1:length(t)
            SignalPol(w) = trapz(Presign(w,:));
            SignalRef(w) = trapz(Presignref(w,:));
        end
        
     
        % The normalized values of Signal Pol and Ref
        for c = 1:length(t)
            NormSignalPol(c) = SignalPol(c)./SignalPol(end);  %Range
            NormSignalRef(c) = SignalRef(c)./SignalRef(end);
        end

        signals.S_off{i} = SignalRef;
        signals.S_on{i} = SignalPol;
        signals.t{i} = t;
    end

%% ============================================================
%  PDE DEFINITIONS
%  c(x,t,u,ux)*u_t = x^{-m} * d/dx [x^m * f(x,t,u,ux)] + s
%% ============================================================

    % Source term WITH DNP.
    function [c,f,s] = pdex1pde(x,~,u,DuDx)
        c = Cp*fc(x);
        f = fD(x)*Cp*fc(x)*DuDx;
        s = -(u-fe0(x))*(1/fT1(x))*Cp*fc(x); % DNP ON: target = epsilon*dep
    end

    % Source term WITHOUT DNP.
    function [c,f,s] = pdex1pdenoDNP(x,~,u,DuDx)
        c = Cp*fc(x);
        f = fD(x)*Cp*fc(x)*DuDx;
        s = -(u-depo(x))*(1/fT1(x))*Cp*fc(x); % DNP OFF: target = dep
    end

    % Initial polarization is 0 all over the crystal WITH DNP
    function u0 = pdex1ic(~)
        u0 = 0;
    end

    % Initial polarization is 0 all over the crystal WITHOUT DNP
    function u0 = pdex1icnoDNP(~)
        u0 = 0;
    end

    % Boundary conditions WITH DNP
    function [pl,ql,pr,qr] = pdex1bc(~,~,~,~,~)
        pl = 0; ql = 1; pr = 0; qr = 1;
    end

    % Boundary condition WITHOUT DNP
    function [pl,ql,pr,qr] = pdex1bcnoDNP(~,~,~,~,~)
        pl = 0; ql = 1; pr = 0; qr = 1;
    end

%% ============================================================
%  SPATIAL PROFILE FUNCTIONS  (smooth tanh transitions)
%% ============================================================

% T1 function that will determine the relaxation term
    function T1 = fT1(y)
        T1 = (T1in+T1out)/2 + (T1in-T1out)/2 * tanh(slope*(length_target-y));
    end

% Build up epsilon spatial function MW ON
    function [epsilonzero] = fe0(y)
        epsilonzero = (epsilon0*dep+1)/2 + (epsilon0*dep-1)/2 * tanh(-slope*(length_target-y));
    end

% Contribution factor (quenching) spatial function
    function [theta3] = theta0(y) 
        theta3 = (1+theta)/2 + (1-theta)/2 * tanh(slope*(length_target-y));
    end

% Spin diffusion spatial function
    function [D] = fD(y)
        D = (Din+Dout)/2 + (Din-Dout)/2 * tanh(slope*(length_target-y));
    end

% Build−up spatial function MW OFF (depolarization)
    function [depol] = depo(y)
        depol = (1+dep)/2 + (1-dep)/2 * tanh(slope*(length_target-y));
    end

% Concentration spatial function
    function [c] = fc(y)
        c = (concX+concjuice)/2 + (concX-concjuice)/2 * tanh(slope*(length_target-y));
    end

end