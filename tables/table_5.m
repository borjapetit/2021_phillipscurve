

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all ; close all ; clc ; tic; format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Working directory

cd ' %% Your Own Directory '

run('matlab/function/parameters.m');

cases = [ 5 2 0 3 4 ] ; 

freq_p     = zeros(size(cases)) ;
aveabsch_p = zeros(size(cases)) ;
skech_p    = zeros(size(cases)) ;
kurch_p    = zeros(size(cases)) ;
std_p      = zeros(size(cases)) ;
ch0_p      = zeros(size(cases)) ;
absch25_p  = zeros(size(cases)) ;
losses_p   = zeros(size(cases)) ;
errors_p   = zeros(size(cases)) ;
setting_p  = zeros(size(cases)) ;
timing_p   = zeros(size(cases)) ;

freq_w     = zeros(size(cases)) ;
aveabsch_w = zeros(size(cases)) ;
skech_w    = zeros(size(cases)) ;
kurch_w    = zeros(size(cases)) ;
std_w      = zeros(size(cases)) ;
ch0_w      = zeros(size(cases)) ;
absch25_w  = zeros(size(cases)) ;
losses_w   = zeros(size(cases)) ;
errors_w   = zeros(size(cases)) ;
setting_w  = zeros(size(cases)) ;
timing_w   = zeros(size(cases)) ;

totalloss  = zeros(size(cases)) ;

for i = 1:size(cases,2)
    
    sol_ss = importdata([pwd '/textfiles/_ss/V1' num2str(cases(i)) '_ss.txt'],' ',0);
    run('matlab/function/extract_ss.m');

    freq_p(i)     = stat.freq_p ;
    aveabsch_p(i) = stat.aveabsch_p ;
    skech_p(i)    = stat.skech_p ;
    kurch_p(i)    = stat.kurch_p ;
    std_p(i)      = stat.std_p ;
    ch0_p(i)      = stat.ch0_p ;
    absch25_p(i)  = stat.absch25_p ;
    losses_p(i)   = stat.losses_p ;
    errors_p(i)   = stat.errors_p ;
    setting_p(i)  = stat.setting_p ;
    timing_p(i)   = stat.timing_p ;

    freq_w(i)     = stat.freq_w ;
    aveabsch_w(i) = stat.aveabsch_w ;
    skech_w(i)    = stat.skech_w ;
    kurch_w(i)    = stat.kurch_w ;
    std_w(i)      = stat.std_w ;
    ch0_w(i)      = stat.ch0_w ;
    absch25_w(i)  = stat.absch25_w ;
    losses_w(i)   = stat.losses_w ;
    errors_w(i)   = stat.errors_w ;
    setting_w(i)  = stat.setting_w ;
    timing_w(i)   = stat.timing_w ;
    
    totalloss(i)   = stat.totalloss ;    
    
    clear sol_ss stat
end

clc ; 

fprintf('\n')
fprintf('\n')
fprintf('    TABLE 5. EVALUATING THE MODEL WITH DIFFERENT TREND INFLATION RATES')
fprintf('\n')
fprintf('\n')
fprintf('    ***************************************************************** \n')
fprintf('    PRICES                        -2%%      0%%      2%%      4%%      8%% \n')
fprintf('    ***************************************************************** \n')
fprintf('    Frequency of changes:      %6.2f  %6.2f  %6.2f  %6.2f  %6.2f      \n', freq_p     )
fprintf('    Average absolute change:   %6.2f  %6.2f  %6.2f  %6.2f  %6.2f      \n', aveabsch_p )
fprintf('    Skewness of change:        %6.2f  %6.2f  %6.2f  %6.2f  %6.2f      \n', skech_p    )
fprintf('    Kurtosis of change:        %6.2f  %6.2f  %6.2f  %6.2f  %6.2f      \n', kurch_p    )
fprintf('    Std.deviation prices:      %6.2f  %6.2f  %6.2f  %6.2f  %6.2f      \n', std_p      )
fprintf('    ----------------------------------------------------------------- \n'             )
fprintf('    %% changes > 0:             %6.2f  %6.2f  %6.2f  %6.2f  %6.2f     \n', ch0_p      )
fprintf('    %% abs. changes < 2.5%%:     %6.2f  %6.2f  %6.2f  %6.2f  %6.2f    \n', absch25_p  )
fprintf('    ----------------------------------------------------------------- \n'             )
fprintf('    Output losses (%%):         %6.2f  %6.2f  %6.2f  %6.2f  %6.2f     \n', losses_p   )
fprintf('       Cost due errors (%%):    %6.2f  %6.2f  %6.2f  %6.2f  %6.2f     \n', errors_p   )
fprintf('       Cost setting    (%%):    %6.2f  %6.2f  %6.2f  %6.2f  %6.2f     \n', setting_p  )
fprintf('       Cost setting    (%%):    %6.2f  %6.2f  %6.2f  %6.2f  %6.2f     \n', timing_p   )
fprintf('    ***************************************************************** \n'             )
fprintf('    WAGES                        Data   Base.      FP      FW    FPFW \n'             )
fprintf('    ***************************************************************** \n'             )
fprintf('    Frequency of changes:      %6.2f  %6.2f  %6.2f  %6.2f  %6.2f      \n', freq_w     )
fprintf('    Average absolute change:   %6.2f  %6.2f  %6.2f  %6.2f  %6.2f      \n', aveabsch_w )
fprintf('    Skewness of change:        %6.2f  %6.2f  %6.2f  %6.2f  %6.2f      \n', skech_w    )
fprintf('    Kurtosis of change:        %6.2f  %6.2f  %6.2f  %6.2f  %6.2f      \n', kurch_w    )
fprintf('    Std.deviation wages:       %6.2f  %6.2f  %6.2f  %6.2f  %6.2f      \n', std_w      )
fprintf('    ----------------------------------------------------------------- \n'             )
fprintf('    %% changes > 0:             %6.2f  %6.2f  %6.2f  %6.2f  %6.2f     \n', ch0_w      )
fprintf('    %% abs. changes < 2.5%%:     %6.2f  %6.2f  %6.2f  %6.2f  %6.2f    \n', absch25_w  )
fprintf('    ----------------------------------------------------------------- \n'             )
fprintf('    Output losses (%%):         %6.2f  %6.2f  %6.2f  %6.2f  %6.2f     \n', losses_w   )
fprintf('       Cost due errors (%%):    %6.2f  %6.2f  %6.2f  %6.2f  %6.2f     \n', errors_w   )
fprintf('       Cost setting    (%%):    %6.2f  %6.2f  %6.2f  %6.2f  %6.2f     \n', setting_w  )
fprintf('       Cost setting    (%%):    %6.2f  %6.2f  %6.2f  %6.2f  %6.2f     \n', timing_w   )
fprintf('    ***************************************************************** \n'             )
fprintf('    Total losses (%%):          %6.2f  %6.2f  %6.2f  %6.2f  %6.2f     \n', totalloss  )
fprintf('    ***************************************************************** \n'             )
fprintf('\n')
fprintf('\n')
fprintf('\n')

