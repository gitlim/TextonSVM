%MATLAB FILE TO TEST THE FAST CLASSIFIER%
%Subhransu Maji (smaji@cs.berkeley.edu)
%Creates a random training and test data set. Expected error is
%close to 50%. Here however we are testing how close to the exact
%predictions do we get using the approximations.

%training data%
numtr = 1000;
numte = 50000;
dim = 100;
trd = rand(numtr,dim);
trl = (rand(numtr,1) > 0.5)+1;
%test data%
ted = rand(numte,dim);
tel = (rand(numte,1) > 0.5)+1;

model = svmtrain(trl,trd,'-s 0 -t 4');
tic;
[p,acc,dec] = svmpredict(tel,ted,model);
svmtime=toc;

nbins = [10 20 30 40 50 80 100 150 200 250 300];
for i = 1:length(nbins)
  teststring = sprintf('-v 0 -n %d',nbins(i));
  [fe,fpwc,fpwl,ft(i,:)] = fastpredict(tel,ted,model,teststring); 
  err_e(i)   = mean(abs(dec-fe));
  err_pwc(i) = mean(abs(dec-fpwc));
  err_pwl(i) = mean(abs(dec-fpwl));
end
figure;
hold on;
plot(nbins,err_e,'--r*');
plot(nbins,err_pwc,'--g*');
plot(nbins,err_pwl,'--b*');
legend('exact (binary search)','approx (piecewise constant)',['approx (piecewise ' ...
                    'linear)']);
title('Error in Approximation of Exact vs. Number of Bins');

xlabel('Number of Bins');
ylabel('Error in Approximation of Exact');
grid on;
set(gca,'XTick',nbins)

figure;
hold on;
plot(nbins,ft(:,1),'--k*');
plot(nbins,ft(:,2),'--r*');
plot(nbins,ft(:,3),'--g*');
plot(nbins,ft(:,4),'--b*');
legend('precomp','exact (binary search)','approx (piecewise constant)',['approx (piecewise ' ...
                    'linear)']);

title(sprintf('Classification Time (libsvm=%.2fs , %d examples, %d sup vec, %d dim)',svmtime,numte,model.totalSV,dim));
xlabel('Number of Bins');
ylabel('Classification Time(s)');
grid on;
set(gca,'XTick',nbins)
