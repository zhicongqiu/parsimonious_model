plot(2:500,Log_incomplete(2:end),'linewidth',2 );
hold on;
[temp1 temp2] = min(Log_incomplete(2:end));
plot(temp2+1,temp1,'xr','linewidth',2);
legend('BIC incomplete Log Likelihood','minimum');
xlabel('number of mixture components');
title('BIC cost as a function of number of mixture components,argmin=317');
hold off;
print('BIC_cost.png','-dpng');
close(gcf);

plot(2:500,error(500:-1:2),317,error(500-316),'xr','linewidth',2 );
xlabel('number of mixture components');
title('classification error rate, error(317) = 0.1933');
print('error_rate.png','-dpng');
close(gcf);
