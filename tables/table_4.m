
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all ; close all ; clc ; tic; format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Working directory

cd ' %% Your Own Directory '

run('matlab/function/parameters.m');

for version_kappas=[ 1 3 5 6 ]
    sol_ss = importdata([pwd '/textfiles/_ss/V' num2str(version_kappas) '0_ss.txt'],' ',0);
    run('matlab/function/extract_ss.m');
    eval([ 'statV' num2str(version_kappas) ' = stat ; ' ])
    clear sol_ss stat
end

clc ; 

fprintf('\n')
fprintf('\n')
fprintf('\n')
fprintf('    TABLE 4. EVALUATING THE MODEL WITH DIFFERENT NOISE PARAMETERS')
fprintf('\n')
fprintf('\n')
fprintf('    ***************************************************************** \n')
fprintf('    PRICES                       Data   Base.      FP      FW    FPFW \n')
fprintf('    ***************************************************************** \n')
fprintf('    Frequency of changes:      %6.2f  %6.2f  %6.2f  %6.2f  %6.2f     \n',  10.2 , statV1.freq_p     , statV3.freq_p,statV5.freq_p,statV6.freq_p )
fprintf('    Average absolute change:   %6.2f  %6.2f  %6.2f  %6.2f  %6.2f     \n',  9.90 , statV1.aveabsch_p , statV3.aveabsch_p,statV5.aveabsch_p,statV6.aveabsch_p )
fprintf('    Skewness of change:        %6.2f  %6.2f  %6.2f  %6.2f  %6.2f     \n', -0.42 , statV1.skech_p    , statV3.skech_p,statV5.skech_p,statV6.skech_p )
fprintf('    Kurtosis of change:        %6.2f  %6.2f  %6.2f  %6.2f  %6.2f    \n',  4.81 , statV1.kurch_p    , statV3.kurch_p,statV5.kurch_p,statV6.kurch_p )
fprintf('    Std.deviation prices:          --  %6.2f  %6.2f  %6.2f  %6.2f     \n',  statV1.std_p    , statV3.std_p,statV5.std_p,statV6.std_p )
fprintf('    -----------------------------------------------------------------   \n')
fprintf('    %% changes > 0:             %6.2f  %6.2f  %6.2f  %6.2f  %6.2f    \n',  65.1 , statV1.ch0_p      , statV3.ch0_p,statV5.ch0_p,statV6.ch0_p )
fprintf('    %% abs. changes < 2.5%%:     %6.2f  %6.2f  %6.2f  %6.2f  %6.2f   \n',  12.0 , statV1.absch25_p  , statV3.absch25_p,statV5.absch25_p,statV6.absch25_p)
fprintf('    -----------------------------------------------------------------   \n')
fprintf('    Output losses (%%):             --  %6.2f  %6.2f  %6.2f  %6.2f   \n', statV1.losses_p  , statV3.losses_p  , statV5.losses_p  , statV6.losses_p  )
fprintf('       Cost due errors (%%):        --  %6.2f  %6.2f  %6.2f  %6.2f   \n', statV1.errors_p  , statV3.errors_p  , statV5.errors_p  , statV6.errors_p  )
fprintf('       Cost setting    (%%):        --  %6.2f  %6.2f  %6.2f  %6.2f   \n', statV1.setting_p , statV3.setting_p , statV5.setting_p , statV6.setting_p )
fprintf('       Cost setting    (%%):        --  %6.2f  %6.2f  %6.2f  %6.2f   \n', statV1.timing_p  , statV3.timing_p  , statV5.timing_p  , statV6.timing_p  )
fprintf('    ***************************************************************** \n')
fprintf('    WAGES                        Data   Base.      FP      FW    FPFW \n')
fprintf('    ***************************************************************** \n')
fprintf('    Frequency of changes:      %6.2f  %6.2f  %6.2f  %6.2f  %6.2f     \n',8.30,statV1.freq_w,statV3.freq_w,statV5.freq_w,statV6.freq_w )
fprintf('    Average absolute change:   %6.2f  %6.2f  %6.2f  %6.2f  %6.2f     \n',6.47,statV1.aveabsch_w,statV3.aveabsch_w,statV5.aveabsch_w,statV6.aveabsch_w )
fprintf('    Skewness of change:        %6.2f  %6.2f  %6.2f  %6.2f  %6.2f     \n',0.35,statV1.skech_w,statV3.skech_w,statV5.skech_w,statV6.skech_w )
fprintf('    Kurtosis of change:        %6.2f  %6.2f  %6.2f  %6.2f  %6.2f     \n',4.39,statV1.kurch_w,statV3.kurch_w,statV5.kurch_w,statV6.kurch_w )
fprintf('    Std.deviation wages:           --  %6.2f  %6.2f  %6.2f  %6.2f     \n',statV1.std_w,statV3.std_w,statV5.std_w,statV6.std_w )
fprintf('    -----------------------------------------------------------------   \n')
fprintf('    %% changes > 0:             %6.2f  %6.2f  %6.2f  %6.2f  %6.2f    \n',86.5,statV1.ch0_w,statV3.ch0_w,statV5.ch0_w,statV6.ch0_w )
fprintf('    %% abs. changes < 2.5%%:     %6.2f  %6.2f  %6.2f  %6.2f  %6.2f   \n',11.8,statV1.absch25_w,statV3.absch25_w,statV5.absch25_w,statV6.absch25_w)
fprintf('    -----------------------------------------------------------------   \n')
fprintf('    Output losses (%%):             --  %6.2f  %6.2f  %6.2f  %6.2f   \n', statV1.losses_w,statV3.losses_w,statV5.losses_w,statV6.losses_w)
fprintf('       Cost due errors (%%):        --  %6.2f  %6.2f  %6.2f  %6.2f   \n', statV1.errors_w,statV3.errors_w,statV5.errors_w,statV6.errors_w)
fprintf('       Cost setting    (%%):        --  %6.2f  %6.2f  %6.2f  %6.2f   \n', statV1.setting_w,statV3.setting_w,statV5.setting_w,statV6.setting_w)
fprintf('       Cost setting    (%%):        --  %6.2f  %6.2f  %6.2f  %6.2f   \n', statV1.timing_w,statV3.timing_w,statV5.timing_w,statV6.timing_w)
fprintf('    ***************************************************************** \n')
fprintf('    Total losses (%%):              --  %6.2f  %6.2f  %6.2f  %6.2f   \n', statV1.totalloss,statV3.totalloss,statV5.totalloss,statV6.totalloss)
fprintf('    ***************************************************************** \n')
fprintf('\n')
fprintf('\n')
fprintf('\n')



