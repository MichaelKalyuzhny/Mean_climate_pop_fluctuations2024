%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Population model with fluctuating environment - Gompertz model with power env. response curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Michael Kalyuzhny
%
% Written: 23/06/23
% Modified: 13/07/23
%
%% Generate time series:

% Set parameters:
b = 0.6; %DD
alpha = 1;
betas = [1 0.2 2];
E_means = [0.6 2];
E_half_range = 0.55;
t_tot = 100;

% Preallocate:
sp_num = length(E_means);
betas_num = length(betas);
Ns = nan(t_tot,sp_num,betas_num);

% run dynamics:

for bb = 1:betas_num
    for ss = 1:sp_num %run for every species:

        Ns(1,ss,bb) = (alpha*E_means(ss)^betas(bb))/(1-b); %assign initial abundance

        for tt = 2:t_tot
            E = unifrnd(E_means(ss) - E_half_range, E_means(ss) + E_half_range, 1); %environmental condition this time step
            Ns(tt,ss,bb) = alpha*E^betas(bb) + b*Ns(tt - 1,ss,bb);
        end

    end
end

%% Plot:

figure()

plot_titles = {'Linear', 'Sublinear', 'Superlinear'};
texts = {'(a)', '(b)', '(c)'};
% Plot:
for bb = 1:betas_num 
    s1 = subplot(1,3,bb);
    set(s1,'FontSize',16,'FontWeight','bold','LineWidth',1)
    hold on
    plot(Ns(:,:,bb))
    title(plot_titles{bb})
    
    if bb == 1
        ylabel('Log(abundance)')
    elseif bb == 2
        xlabel('Time')
    end
    
    text(2,2,texts{bb},'FontSize',14,'FontWeight','bold')

end
% subplot(1,3,1)
% 
% plot(Ns')
% subplot(1,2,2)
% plot(exp(Ns'))
% disp(['1-y fluctuations variance:'])
% disp(['species 1: ' num2str(var(Ns(1,2:end) - Ns(1,1:end-1)))])
% disp(['species 2: ' num2str(var(Ns(2,2:end) - Ns(2,1:end-1)))])

%% Examine analytical results:

b = 0.9; %DD
alpha = 1;
beta = 0.7;
E_mean = 0.6;
E_half_range = 0.55;
t_tot = 1000000;
var_e = (1/12)*(2*E_half_range)^2; %variance 


eqilibrium = (alpha*E_mean^beta + 0.5*var_e*alpha*beta*(beta-1)*E_mean^(beta-2))/(1-b);
Ns = nan(1,t_tot);
Ns(1) = eqilibrium; %assign initial abundance

for tt = 2:t_tot
    E = unifrnd(E_mean - E_half_range, E_mean + E_half_range, 1); %environmental condition this time step
    Ns(tt) = alpha*E^beta + b*Ns(tt - 1);
end

figure()
plot(1:t_tot,Ns)

eqilibrium
observed_mean = mean(Ns)

expected_variance = (var_e*(alpha*beta*E_mean^(beta-1))^2)/(1-b^2)
observed_variance = var(Ns)

expected_diff_variance = (var_e*(alpha*beta*E_mean^(beta-1))^2)*(2/(1+b))
observed_diff_variance = var(diff(Ns))