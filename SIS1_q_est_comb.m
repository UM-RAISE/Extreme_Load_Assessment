% 50 repetitions are combined: 1/10/2015
% SIS1 Quantile estimation
% CMC quantile estimation code is in rechts -- C:\Research\Implementation\CMC_IO\CMC_data\CMC_q_est.m


% Reference codes:
% C:\Users\John_and_Annie\Dropbox\Research\Implementation\Code\ControlledSIS_test_SIS1_POE.m
% C:\Users\John_and_Annie\Dropbox\Research\Implementation\Code\func_quantile_est_IS.m
% C:\Users\John_and_Annie\Dropbox\Research\Implementation\Code\NREL_metamodel_IS_ed5.m

%load('C:\Users\John_and_Annie\Dropbox\Research\Implementation\Code\ControlledSIS_test_SIS1edge9300.mat')
load('C:\Research\Implementation\FluxIO2\ChosenData\ChosenData_SIS1edge9300.mat')



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
inside_of_sqrt = (@(x)GEV_P_given_x(x)/N_total + (1-1/N_total).*(GEV_P_given_x(x)).^2);
integrand = (@(x)raylpdf(x,Rayleigh_scale_par).*sqrt(inside_of_sqrt(x))./(raylcdf(25,Rayleigh_scale_par)-raylcdf(3,Rayleigh_scale_par)));
norm_constant_IS_pdf = quadgk(integrand,3,25); %quadgk(integrand,0,Inf);
    %[q,errbnd] =quadgk(integrand,3,25)
        
    
 
%% Deliver the POE estimate for l_T (=target_quantile)
tic
old_digits = digits;   % For Variable-precision accuracy
digits(50);


%Y ordered
All_Y_ordered = zeros([1, N_reps_of_choice*N_total]);
for index_repetition = 1:N_reps_of_choice
    if LoadType == 1    %flapwise
        Y_i_j_matrix = squeeze(Chosen_Y_flap(index_repetition,:,:));
    elseif LoadType == 2 %edgewise     
        Y_i_j_matrix = squeeze(Chosen_Y_edge(index_repetition,:,:));
    end       

    Y_nonzero = nonzeros(Y_i_j_matrix);
    All_Y_ordered(((index_repetition-1)*N_total+1):(index_repetition*N_total)) = Y_nonzero';
end
All_Y_ordered = sort(All_Y_ordered);



All_POE_estimates = zeros([1, N_reps_of_choice*N_total]);



for index_repetition = 1:N_reps_of_choice

    
    % Determine Optimal allocation of short-term sample size (Note: rounding error)
    allocation_proportion_for_obs_x = repmat(vpa(0.0),[1,N_long+1]);
    for index_i = 1:N_long  
        P = GEV_P_given_x(Chosen_x(index_repetition,index_i));
        allocation_proportion_for_obs_x(index_i) = vpa(sqrt( (1-P)     / ( 1+(N_total-1)*P )    ) );   %proportional constant for the observed x_i =(f/q)*sqrt(P(Y>y_T | X=x)*(1-P(Y>y_T | X=x))) = sqrt( (1-P)     / ( 1+(N_total-1)*P )    )      : 'Meeting_080212.pptx' p.14.
    end
    alloc_normalizing_constant = vpa(0.0);
    for index_i = 1:N_long  
        alloc_normalizing_constant = vpa(alloc_normalizing_constant + allocation_proportion_for_obs_x(index_i));    
    end


    % Sort the proportions for x_i so that the most important (higest proportion) x_i has the first index.
    [sorted_proportion, original_order]=sort(allocation_proportion_for_obs_x, 'descend');
    temp_N_short_allocation = zeros(1, N_long); 
    adjusted_normalizing_constant = alloc_normalizing_constant;
    for index = N_long:-1:2 %drange(N_long:-1:2)
        remained_sample_size = N_total-sum(temp_N_short_allocation(index+1:N_long));   %N - [  sum_j=(k+1)^(N_long) N_j ]
        adjusted_normalizing_constant =  vpa(adjusted_normalizing_constant - sorted_proportion(index+1));
        unrounded_alloc_for_obs_x = vpa(remained_sample_size * sorted_proportion(index)/adjusted_normalizing_constant);
        if (double(unrounded_alloc_for_obs_x)<1)
            temp_N_short_allocation(index) = 1;
        else    %(double(unrounded_alloc_for_obs_x)>=1)
            temp_N_short_allocation(index) = round(double(unrounded_alloc_for_obs_x));
        end       
    end 
    temp_N_short_allocation(1) = N_total - sum(temp_N_short_allocation(2:N_long));    %N_total - [  sum_j=2^(N_long) N_j ]


    N_short_allocation  = zeros(1, N_long);   % N_i : the allocated sample size for each stratum i
    % Order to recover the original indice
    for index = 1:N_long    %drange(1:N_long)
        N_short_allocation(original_order(index)) = temp_N_short_allocation(index);
    end 

    
    if LoadType == 1    %flapwise
        Y_i_j_matrix = squeeze(Chosen_Y_flap(index_repetition,:,:));
    elseif LoadType ==2 %edgewise     
        Y_i_j_matrix = squeeze(Chosen_Y_edge(index_repetition,:,:));
    end       

    likelihood_at_x = vpa(1./sqrt(inside_of_sqrt(Chosen_x(index_repetition,:))));       

    %----------------------------------------   
    %POE estimates at every observed Y
    Y_ordered = nonzeros(Y_i_j_matrix);
    Y_ordered = sort(Y_ordered);
    
    [~,Y_loc] = ismember(Y_ordered,All_Y_ordered);   
    
    if(All_Y_ordered(1) < 0)
        error('The minimum quantile is less than 0.') %For vectorized computation of Ind_Y_i_j_matrix, I assumed the threshold level >=0, because Y_i_j_matrix coded NA with zeros.
    end    

    if(Y_loc(1) > 1)
        All_POE_estimates(1:(Y_loc(1)-1)) = All_POE_estimates(1:(Y_loc(1)-1)) + sum(likelihood_at_x) / N_long * norm_constant_IS_pdf; %For l less than Y_ordered(1), all(sum_for_short_term==N_short_allocation) ==1, so that (sum_for_long_term_t == sum(likelihood_at_x))
    end
    

    for index = 1:(N_total-1) %Vectorizing this part leads to 'Out of memory'
        Ind_Y_i_j_matrix = (Y_i_j_matrix > All_Y_ordered(Y_loc(index))); %Output Ind_Y_ij = I(Y_ij > y_T) for j-th observation of Y for X = x_i
        sum_for_short_term = sum(Ind_Y_i_j_matrix,2)'; %\sum_(j=1)^(N_i) I(Y_ij > y_T)  
        sum_for_long_term_t = sum(likelihood_at_x./ N_short_allocation.* sum_for_short_term);  %vpa(sum_for_long_term + term_for_short_term(index_i));        
        All_POE_estimates(Y_loc(index):(Y_loc(index+1)-1)) = All_POE_estimates(Y_loc(index):(Y_loc(index+1)-1)) + sum_for_long_term_t / N_long * norm_constant_IS_pdf; %Note: (1/C_q)* N_long = expectation( sum_of_likelihood )  %Differ in order of e-60 from the result using vpa.
    end   

    %All_POE_estimates(Y_loc(N_total):end) = All_POE_estimates(Y_loc(N_total):end) + 0;  %, because sum(sum_for_short_term) = 0
    
%%Legacy code:  Takes too long without using ismember.   1/10/2015
%     for ind_r = 1:N_reps_of_choice
%         start_zero_ind = (ind_r-1)*N_total;
%         for index = 1:N_total %Vectorizing this part leads to 'Out of memory'
%             Ind_Y_i_j_matrix = (Y_i_j_matrix > All_Y_ordered(start_zero_ind+index)); %Output Ind_Y_ij = I(Y_ij > y_T) for j-th observation of Y for X = x_i
%             sum_for_short_term = sum(Ind_Y_i_j_matrix,2)'; %\sum_(j=1)^(N_i) I(Y_ij > y_T)  
%             sum_for_long_term_t = sum(likelihood_at_x./ N_short_allocation.* sum_for_short_term);  %vpa(sum_for_long_term + term_for_short_term(index_i));        
%             All_POE_estimates(start_zero_ind+index) = All_POE_estimates(start_zero_ind+index) + sum_for_long_term_t / N_long * norm_constant_IS_pdf; %Note: (1/C_q)* N_long = expectation( sum_of_likelihood )  %Differ in order of e-60 from the result using vpa.
%         end   
%     end

    %----------------------------------------      

    
    fprintf('Done: %d / %d\n', index_repetition, N_reps_of_choice) 

end

All_POE_estimates = All_POE_estimates ./ N_reps_of_choice;


% 
% fprintf('%g\n', mean(double(CV_POE_estimates)))
% fprintf('%g\n', std(double(CV_POE_estimates)))
% 


digits(old_digits); %Recover default accuracy: 32bit
toc  %Elapsed time is 15836.296501 seconds. for  'SIS1_q_est_comb.m' generating 'ExceedancePairs_from_SIS1edge9300_50Combined.mat'.



save('ExceedancePairs_from_SIS1edge9300_50Combined.mat','All_Y_ordered', 'All_POE_estimates')







%% Plot the pairs (POEs and quantiles)to make the long-term exceedence plots in C:\Users\John_and_Annie\Dropbox\Research\Papers\Wind energy\Database for validation of design load extrapolation techniques.pdf