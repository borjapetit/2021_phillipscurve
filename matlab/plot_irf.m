
MarkerSize = 5    ;
start      = 2    ;
endp       = 22   ;
time       = [0:1:length(irf.money(start:endp))-1] ;

subplot(2,3,1)
    plot(time,irf.money(start:endp),linescolor,'MarkerSize', MarkerSize)
    ylabel('Money growth')
    xlabel('Months')
    xlim([0 20])
    %ylim([-0.1 1])
    hold on
subplot(2,3,2)
    plot(time,cumsum(irf.pinflation(start:endp)),linescolor,'MarkerSize', MarkerSize)
    ylabel('Price level')
    xlabel('Months')
    %ylim([-0.1 3])
    xlim([0 20])
    hold on
subplot(2,3,3)
    plot(time,irf.consumption(start:endp),linescolor,'MarkerSize', MarkerSize)
    ylabel('Consumption')
    xlabel('Months')
    %ylim([-0.1 3])
    xlim([0 20])
    hold on
subplot(2,3,4)
    plot(time,irf.labor(start:endp),linescolor,'MarkerSize', MarkerSize)
    ylabel('Labor')
    xlabel('Months')
    %ylim([-0.1 3])
    xlim([0 20])
    hold on
subplot(2,3,5)
    plot(time,cumsum(irf.winflation(start:endp)),linescolor,'MarkerSize', MarkerSize)
    ylabel('Wage level')
    xlabel('Months')
    %ylim([-0.5 3.5])
    xlim([0 20])
    hold on
subplot(2,3,6)
    plot(time,irf.wage(start:endp),linescolor,'MarkerSize', MarkerSize)
    ylabel('Real wage')
    xlabel('Months')
    %ylim([-0.1 3])
    xlim([0 20])
    hold on
