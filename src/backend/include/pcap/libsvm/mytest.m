f = rand(15,10);
l = ones(15,1);
l(6:end) = -1;
model = svmtrain(l,f,'-t 4');
x = precomp_exact_model(model);
