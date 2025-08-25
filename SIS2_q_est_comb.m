% 50 repetitions are combined: 1/22/2015
% SIS2 Quantile estimation
% CMC quantile estimation code is in rechts -- C:\Research\Implementation\CMC_IO\CMC_data\CMC_q_est.m


% Reference codes:
% C:\Users\John_and_Annie\Dropbox\Research\Implementation\Code\ControlledSIS_test_SIS1_POE.m
% C:\Users\John_and_Annie\Dropbox\Research\Implementation\Code\func_quantile_est_IS.m
% C:\Users\John_and_Annie\Dropbox\Research\Implementation\Code\NREL_metamodel_IS_ed5.m
% C:\Users\John_and_Annie\Dropbox\Research\Writing\ExtLoad\Code\SIS1\SIS1_q_est.m

clear
% load('C:\Users\John_and_Annie\Dropbox\Research\Implementation\Code\ControlledSIS_test_SIS2edge9300_0913_numVar.mat')
% %load('C:\Joun\Dropbox\Research\Implementation\Code\ControlledSIS_test_SIS2edge9300_0913_numVar.mat')
% 
% clearvars -except Chosen_Y_edge Chosen_Y_flap Chosen_x LoadType N_long N_reps_of_choice target_quantile Rayleigh_scale_par
load('C:\Research\Implementation\FluxIO2\ChosenData\ChosenData_SIS2edge9300.mat')


%% Obtain individual indicator funtion outputs: I(Y_ij > l )
IO_folder = 'C:\Research\Implementation\CMC_IO\';   

% GEV parameters

%Loading: GEV_loc_par_at_x  GEV_scal_par_at_x   shape_param
if LoadType == 1    %flapwise
        load(horzcat(IO_folder, 'Flap_GAMLSS_param.mat')); 
elseif LoadType ==2 %edgewise     
        load(horzcat(IO_folder, 'Edge_GAMLSS_param.mat'));
end 
        

%% P_meta_given_x
GEV_P_given_x = @(x) 1-gevcdf(target_quantile,shape_param,GEV_scal_par_at_x(x),GEV_loc_par_at_x(x));


%% Calculate the normalizing constant of IS distribution
integrand = (@(x)raylpdf(x,Rayleigh_scale_par).*sqrt(GEV_P_given_x(x))./(raylcdf(25,Rayleigh_scale_par)-raylcdf(3,Rayleigh_scale_par)));
norm_constant_IS_pdf = quadgk(integrand,3,25); %quadgk(integrand,0,Inf);
    %[q,errbnd] =quadgk(integrand,3,25)
    
 
%% Deliver the POE estimate for l_T (=target_quantile)
tic
old_digits = digits;   % For Variable-precision accuracy
digits(50);

N_total = N_long;

%Y ordered
All_Y = zeros([1, N_reps_of_choice*N_total]);
All_X = zeros([1, N_reps_of_choice*N_total]);
for index_repetition = 1:N_reps_of_choice
    if LoadType == 1    %flapwise
        Y_i_j_matrix = Chosen_Y_flap(index_repetition,:);
    elseif LoadType == 2 %edgewise     
        Y_i_j_matrix = Chosen_Y_edge(index_repetition,:);
    end       

    All_Y(((index_repetition-1)*N_total+1):(index_repetition*N_total)) = Y_i_j_matrix;
    All_X(((index_repetition-1)*N_total+1):(index_repetition*N_total)) = Chosen_x(index_repetition,:);    
    
end
%toc
%tic
likelihood_at_x = 1./sqrt(GEV_P_given_x(All_X)); %vpa(1./sqrt(GEV_P_given_x(All_X)));  %without vpa: Tolerable error in total sum: sum(likelihood_at_x) - sum(double(likelihood_at_x))  == 1.3970e-09 
%cumsum_likelihood_at_x = cumsum(likelihood_at_x); %Tolerable error in cumulative sum: cumsum_likelihood_at_x(end) - sum(likelihood_at_x) == -1.8161e-08

%toc %Elapsed time is 12.385740 seconds.
[All_Y_ordered, Ind] = sort(All_Y);
likelihoods_ordered = likelihood_at_x(Ind);
cumsum_likelihoods = fliplr(cumsum(fliplr(likelihoods_ordered)));

All_POE_estimates = [cumsum_likelihoods(2:end) 0].*norm_constant_IS_pdf./(N_reps_of_choice.*N_total);

semilogy(All_Y_ordered, All_POE_estimates)




digits(old_digits); %Recover default accuracy: 32bit
toc  %Elapsed time is 0.159256 seconds.. for  'SIS2_q_est_comb.m' generating 'ExceedancePairs_from_SIS2edge9300_50Combined.mat'.



save('ExceedancePairs_from_SIS2edge9300_50Combined.mat','All_Y_ordered', 'All_POE_estimates')







