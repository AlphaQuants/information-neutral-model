function [res,InsiderOption,BSOption] =PDESolver4()
global T dt t kgrid Xgrid kappaStar;
global tau ;
Upsilon = TerminalCost(tau);

nx= size(Upsilon,1);
nt = T/dt;

InsiderOption =  InsiderOptionPayoff();
BSOption =   InsiderOptionPayoff();

t=T;
nk=size(kgrid,2);
nx =size(Xgrid,2);
UpsilonStar =zeros(nx,nk);
shifted_x_post_execution =zeros(nx,nk);
while t>0
    t= t-dt;
    %for i = 1:5
    [dV_indexchange , UpsilonStar, shifted_x_post_execution ] = OptimizeKappa( Upsilon);
    % kappaStar= make_smooth(kappaStar ,false);
    %end
    
    Upsilon= solve(UpsilonStar,dV_indexchange ,true); %kappaStar
    [InsiderOption] = ReOrganizeGrid(InsiderOption, dV_indexchange, shifted_x_post_execution);
    InsiderOption =  solve(InsiderOption, dV_indexchange,false);
    BSOption = solve(BSOption,zeros(nx,nk),false);
    %applyBarrier(InsiderOption);
    %applyBarrier(BSOption);
    %applyBarrier2(Upsilon);
    %Upsilon=applyTemporaryCost(Upsilon , kappaStar);
    %meshc(kgrid,exp(Xgrid),Upsilon);
    %xlabel('Current Shares (V_0)')
    %ylabel('Spot (S_0)')
    %zlabel('Cost Function J(0,.,.)')
end
res=Upsilon;
end


function Upsilon = applyBoundary(Upsilon)
    Upsilon(1,  :) = 2*Upsilon(2,  :) - Upsilon(3,  :);
    Upsilon(end,:) = 2*Upsilon(end-1,:) - Upsilon(end-2,:);
    
    %Upsilon(:,1) = 2*Upsilon(:,2) - Upsilon(:,3);
    %Upsilon(end,:) = 2*Upsilon(:,end-1) - Upsilon(:,end-2);
end


function [UpsilonNew] = solve(  Upsilon ,dV_indexchange,Smeasure)
    global dt Xgrid kgrid
    nk=size(kgrid,2);
    nx =size(Xgrid,2);
    Upsilon0 =zeros(nx,nk);

    for k =1:nk
        tridiag_matrix=tridiagMat( dV_indexchange(:,k),Smeasure);
        UpsilonVec=tridiag_matrix\(Upsilon(:,k));
        % micro-optim: copy whole column (same numbers)
        Upsilon0(:,k) = UpsilonVec;
        % original equivalent:
        % for l=1:nx
        %     Upsilon0(l,k) = UpsilonVec(l);
        % end
    end
  
    UpsilonNew=applyBoundary(Upsilon0);
end


function tridiag_matrix = tridiagMat(dV_indexchange, Smeasure)
    % TRIDIAGMAT Build tridiagonal matrix for PDE solver
    %   tridiag_matrix = tridiagMat(dV_indexchange, Smeasure)
    %   Smeasure (optional, logical): if true, use stock measure drift tweak.
    %   Default is false.

    global r sigma Xgrid dt dX dk

    nx = numel(Xgrid);
    tridiag_matrix = zeros(nx, nx);

    % Argument to Omega / Lambda
    kappa  = dV_indexchange * dk / dt;
    omega  = Omega(kappa);
    lambda = Lambda(kappa);

    % local aliases / constants (no math change)
    r_loc    = r;
    sigma_loc= sigma;
    sig2     = sigma_loc^2;
    dt_loc   = dt;
    dX_loc   = dX;

    if Smeasure
        driftSign = -1;
    else
        driftSign = 1;
    end

    mid_interior = 1 + (~Smeasure) * r_loc * dt_loc + sig2 * dt_loc / dX_loc^2;
    mid_boundary = 1 + (~Smeasure) * r_loc * dt_loc;

    for i = 1:nx
        if i > 1 && i < nx
            % Interior points
            drift_i = (r_loc - omega(i)*lambda(i) - driftSign*0.5*sig2);
            left  = dt_loc * drift_i / (2*dX_loc) ...
                    - 0.5 * sig2 * dt_loc / dX_loc^2;

            right = -dt_loc * drift_i / (2*dX_loc) ...
                    - 0.5 * sig2 * dt_loc / dX_loc^2;

            tridiag_matrix(i, i-1) = left;
            tridiag_matrix(i, i)   = mid_interior;
            tridiag_matrix(i, i+1) = right;
        else
            % Boundary points (Dirichlet-like)
            tridiag_matrix(i, i) = mid_boundary;
        end
    end
end


function res = Omega(kappa)
% Vectorized, array-safe ?(?):
%   ?(?) = ? * (?^2 - ?(?) * (? - r)) / (?^2 + ? * ?(?)^2)
% Works for scalars, vectors, or matrices of kappa.

    global mu r phi sigma

    L   = Lambda(kappa);                             % same size as kappa
    num = phi .* (sigma.^2 - L .* (mu - r));         % .* elementwise
    den = sigma.^2 + phi .* (L.^2);                  % elementwise

    % numerical guard
    den = max(den, 1e-12);
    res = num ./ den;                                % elementwise divide

    % clean up any accidental non-finites (left as in your code)
end


function [GridValueAfter] = ReOrganizeGrid( GridValue, kappa_d , shifted_x_after_execution)
global dk Xgrid dt Xmin Xmax;
nx= size(GridValue,1);
nk = size(GridValue,2);
GridValueAfter =zeros(nx,nk);

% local aliases (micro-optim, no math change)
dk_loc = dk;
dt_loc = dt;

for index_s = 1:nx
    Xlogspot = Xgrid(index_s); %#ok<NASGU> % kept for consistency, even if unused

    % cache current row to reduce repeated indexing
    row_x = GridValue(index_s,:);

    for k = 1:nk
        kappa_idx = kappa_d(index_s,k);
        kappa = kappa_idx * dk_loc;

        % values at current (x,k)
        GridValue_xk = row_x(k);

        shifted_x = shifted_x_after_execution(index_s,k);
        kp = k + kappa_idx;   % target k index

        if shifted_x < 1 || shifted_x > nx || kp < 1 || kp > nk
            test=1;   % MATLAB will stop here and open debug mode
        end
        % values at shifted locations (same as original code)
        GridValue_xpluskplus = GridValue(shifted_x, kp);
        GridValue_xkplus     = row_x(kp);

        omega = Omega(kappa);

        Lchange = omega*dt_loc * (GridValue_xpluskplus - GridValue_xk) ...
                + (1 - omega*dt_loc) * (GridValue_xkplus - GridValue_xk);

        GridValueAfter(index_s,k) = GridValue_xk + Lchange;
    end
end

GridValueAfter = applyBoundary(GridValueAfter);
end



function [kappa_star_number_discretization, UpsilonStar, shifted_x_after_execution] = OptimizeKappa(Upsilon)
    % Uses globals provided by your code
    global dk Xgrid dt dX Xmin
    nx= size(Upsilon,1);
    nk = size(Upsilon,2);

    UpsilonStar = zeros(nx, nk);
    outsideGrid= zeros(nx,nk);
    kappa_star_number_discretization = zeros(nx, nk);
    shifted_x_after_execution = zeros(nx, nk);

    % ---------- speed: precompute all possible kappa-dependent terms ----------
    % index_diff_k ? {-(nk-1), ..., -1, 0, 1, ..., (nk-1)}
    diff_all = int32(-nk+1:nk-1);
    % your kappa definition (unchanged)
    kappa_all = double(diff_all) * (dk / dt);
    LambdaBar_all = real(LambdaBar(kappa_all));   % 1×(2*nk-1)
    Omega_all     = Omega(kappa_all);             % 1×(2*nk-1)
    Lambda_all    = Lambda(kappa_all);            % 1×(2*nk-1)

    inv_dX = 1.0 / dX;  % small micro-optim

    for index_s = 1:nx
        Xlogspot = Xgrid(index_s);

        % ---------- speed: cache the row once ----------
        Upsilon_row = Upsilon(index_s,:);

        % keep your scalar temps (unchanged names)
        omegaH = 0;
        opt_idx_k = 0;
        u_shifted_kandx = 0;
        u_shifted_k = 0;
        tempC = 0;
        shited_index_x = index_s;

        for k = 1:nk
            opt_idx_k = k;
            upsilon_lk = Upsilon_row(k);
            opt_L     = upsilon_lk;

            % --- choose search direction with edge guards (unchanged logic) ---
            if k == 1
                upordown = +1;
            elseif k == nk
                upordown = -1;
            else
                if Upsilon_row(k+1) > Upsilon_row(k-1)
                    upordown = -1;
                else
                    upordown = +1;
                end
            end

            shifted_index_k = k + upordown;
            shited_index_x = index_s;   % default

            found_local_minima = false;
            while shifted_index_k >= 1 && shifted_index_k <= nk && ~found_local_minima
                index_diff_k = shifted_index_k - k;

                % ---------- speed: lookup precomputed ??, ?, ? ----------
                ix = index_diff_k + nk;              % shift to 1..(2*nk-1)
                deltaX = LambdaBar_all(ix);
                omega  = Omega_all(ix);
                lambda = Lambda_all(ix);
                omegaH = omega * (1 + lambda);

                ShiftedLogSpot = Xlogspot + deltaX;
                shited_index_x_bf_test = (ShiftedLogSpot - Xmin)*inv_dX + 1;
                if shited_index_x_bf_test < 1
                    outsideGrid(index_s,k)=1;
                elseif shited_index_x_bf_test >nx
                    outsideGrid(index_s,k)=1;
                end
                
                % snap to a valid integer grid index (unchanged policy)
                if upordown == 1
                    shited_index_x = min(nx, max(1, floor(shited_index_x_bf_test)));
                else
                    shited_index_x = min(nx, max(1, ceil(shited_index_x_bf_test)));
                end

                % ---------- speed: use cached row & direct col pulls ----------
                u_shifted_kandx = Upsilon(shited_index_x, shifted_index_k);
                u_shifted_k     = Upsilon_row(shifted_index_k);

                % (unchanged) your temp cost call/signature
                tempC = temporyCost( kappa_all(ix)*dt, dt );

                % (unchanged) objective update
                newL  = omega*lambda*dt*(u_shifted_kandx) + tempC  + (1- omega*lambda*dt)*(u_shifted_k);

                if newL < opt_L && ShiftedLogSpot
                    opt_L    = newL;
                    opt_idx_k = shifted_index_k;
                    shifted_index_k = shifted_index_k + upordown;   % keep stepping
                else
                    found_local_minima = true;        % stop at local min
                end
            end

            % (unchanged) final write for this node
            UpsilonStar(index_s, k) =  omegaH*dt*(u_shifted_kandx) + (1-omegaH*dt)*(u_shifted_k) + tempC;
            kappa_star_number_discretization(index_s, k) = opt_idx_k - k;
            shifted_x_after_execution(index_s, k) = shited_index_x;
        end
    end

    % (unchanged) post boundary ops
    UpsilonStar = applyBoundary(UpsilonStar);
    %kappa_star_number_discretization = applyBoundary(kappa_star_number_discretization);
    kappa_star_number_discretization = FixShiftedX(kappa_star_number_discretization, outsideGrid);
    %shifted_x_after_execution = FixShiftedX(shifted_x_after_execution, outsideGrid);
end

function shifted_x_after_execution = FixShiftedX(shifted_x_after_execution, outsideGrid)
% shifted_x_after_execution : nx×nk
% outsideGrid               : nx×nk (1 = outside, 0 = inside)
%
% For each column k:
%   - From the top: for all leading i with outsideGrid(i,k)==1, replace
%     shifted_x_after_execution(i,k) with the first shifted_x_after_execution(j,k)
%     below (j>i) such that outsideGrid(j,k)==0. Stop as soon as we hit the
%     first i where outsideGrid(i,k)==0.
%   - From the bottom: same, but for trailing ones, using j<i and stop at
%     the first inside point when scanning upward.

    [nx, nk] = size(shifted_x_after_execution);

    for k = 1:nk

        % ---------- Forward pass: fix leading outside points ----------
        for i = 1:nx
            if outsideGrid(i,k)
                % find first j > i with outsideGrid(j,k)==0
                j = i + 1;
                while j <= nx && outsideGrid(j,k)
                    j = j + 1;
                end
                if j <= nx
                    shifted_x_after_execution(i,k) = shifted_x_after_execution(j,k);
                    outsideGrid(i,k) = 0;
                end
            else
                % we reached the first inside point; stop forward scan
                break;
            end
        end

        % ---------- Backward pass: fix trailing outside points ----------
        for i = nx:-1:1
            if outsideGrid(i,k)
                % find first j < i with outsideGrid(j,k)==0
                j = i - 1;
                while j >= 1 && outsideGrid(j,k)
                    j = j - 1;
                end
                if j >= 1
                    shifted_x_after_execution(i,k) = shifted_x_after_execution(j,k);
                    outsideGrid(i,k) = 0;
                end
            else
                % we reached the first inside point from the bottom; stop
                break;
            end
        end

    end
end


function [resu] = InsiderOptionPayoff()
global CallStrikeInsiderTrader;
global PutStrikeInsiderTrader;
global KOBarrier;
global Xgrid ;
global kgrid;
nx = size(Xgrid,2);
nk= size(kgrid,2);
resu =zeros(nx,nk);
xlogspot = Xgrid;
for i=1:nx
    ss = exp(xlogspot(i));
    for k=1:nk
        resu(i,k)= max(0,PutStrikeInsiderTrader-ss) + max(0,ss-CallStrikeInsiderTrader);
        if ( ss > KOBarrier )
            resu(i,k)= 0.0;
        end
    end
end
end


function [resu] = TerminalShares()
global CallStrikeLargeTrader;
global PutStrikeLargeTrader;
global KOBarrier;
global Xgrid;
nx = size(Xgrid,2);
resu =zeros(nx,1);
xlogspot = Xgrid;
for i=1:nx
    ss = exp(xlogspot(i));
    resu(i)= 0.0;
    if (ss> CallStrikeLargeTrader)
        resu(i)= 1.0;
    end
    if (ss< PutStrikeLargeTrader)
        resu(i)= 0.0;
    end
    if ( ss > KOBarrier )
        resu(i)= 0.0;
    end
end
end


function [Upsilon] =TerminalCost(liquidationTime)
global Xgrid kgrid dk;
global tau ;
global lambda_term ADV tmpCostAlpha;
nx = size(Xgrid,2);
nk = size(kgrid,2);
spot =Xgrid;
res=zeros(nx,nk);
terminalShares = TerminalShares();
for k=1:nk
    for index_s=1:nx
        chi_exec_rate = (terminalShares(index_s)-kgrid(k) )/liquidationTime;
        tmpC= abs(chi_exec_rate)^(tmpCostAlpha+1)/ (abs(chi_exec_rate) + ADV)^tmpCostAlpha;
        %tmpC= exp(spot(index_s))*abs(deltaDiff)^(tmpCostAlpha+1)/ (abs(deltaDiff) + ADV*liquidationTime)^tmpCostAlpha;
        c= tmpC*lambda_term*liquidationTime;  %% be carefull we added 5 in here  abs(deltaDiff)+
        res(index_s,k)= c;
    end
end
Upsilon =res;
end


function [res] = temporyCost(  deltaDiff , liquidationTime)
global nu ADV tmpCostAlpha
chi_exec_rate = deltaDiff/liquidationTime;
res= nu*abs(chi_exec_rate)^(tmpCostAlpha+1)/(abs(chi_exec_rate)+ADV)^tmpCostAlpha*liquidationTime; 
end

function [res] = Lambda( deltaDiff )
global lambda
res= lambda*(deltaDiff);
end

function [res] = LambdaBar( deltaDiff )
res= log( 1+ Lambda(deltaDiff ));
end
