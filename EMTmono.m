function output=EMTmono(ng,nrms,monocond)
% For some reason this doesn't work for vector ng. Have to pass individual
% ng values one at a time.

% tic
% Inputs:
% ng: in 10^10 cm^-2
% nrms: in 10^10cm^-2

% Output: 
% EMT conductivity in units of e^2/hbar.

%%%%%%%%%%%%%%%% 
% Omitted because it makes more sense to pass monocond as an argument
% from sigmadEMT to prevent repeated loads of the monolayersigmas_TXX.mat
% file.

% load('monolayersigmas_T70K.mat','n', 'sigmamono');
% n=n./1e14; %rescale to units of 10^10 cm^-2.
% nplus=n(2:length(n));
% nfull=[-fliplr(nplus) n];
% 
% sigmaplus=sigmamono(2:length(sigmamono));
% sigmafull=[fliplr(sigmaplus) sigmamono];
% 
% monocond=@(x) interp1(nfull,sigmafull,x,'pchip',NaN);
%%%%%%%%%%%%%%%%%%

P=@(n) exp(-((n-ng).^2)./(2.*nrms.^2)) ; 

topintgrnd=@(n,y) P(n) ./ (1+y./monocond(n));
botintgrnd=@(n,y) P(n) ./ (monocond(n)+y);

uppintgrl=@(y) integral(@(n) topintgrnd(n,y), ng-5.*nrms,ng+5.*nrms,'ArrayValued',true,'RelTol',1e-4);
lowintgrl=@(y) integral(@(n) botintgrnd(n,y), ng-5.*nrms,ng+5.*nrms,'ArrayValued',true,'RelTol',1e-4);
fiter=@(y) uppintgrl(y)./lowintgrl(y);

function Ans = solve(ng, delta)    
        %Guess
        y0 = Inf;
        y1 = monocond(ng);

        while abs((y1 - y0) / y1) >= delta
            y0 = y1;
            y1 = fiter(y1);
        end

        Ans = y1; %* echarge^2 /hbar
end

output=solve(ng,1e-6);
% toc
end