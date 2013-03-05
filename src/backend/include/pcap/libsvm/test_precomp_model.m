%generate random training data and check models
err_cum_alpha=[]; err_cum_alphax=[];err_sort_sv=[];
err_min=[];err_max=[];err_a=[];err_b=[];err_h=[];
err_pwl=[];err_pwc=[]; err_exact=[];
num_iter = 10;
for i = 1:num_iter
  feat_dim = round(max(0.1,rand(1,1))*200);
  num_data = round(max(0.3,rand(1,1))*1000);
  data     = rand(num_data,feat_dim);
  labels   = rand(num_data,1)-0.5;
  labels   = sign(labels);
  model = svmtrain(labels,data,'-t 4');
  exact_model = precomp_model(model,'-m 0');
  %check sorted list
  err_sort_sv(i) = max(max(abs(exact_model.sort_sv-sort(model.SVs)')));
  fprintf('[Iteration #%i] feat_dim:%i , #points:%i, #sup_vec:%i\n',i,feat_dim,num_data,model.totalSV);
  %check cumsum table
      alpha = model.sv_coef;
  clear err_cum_alpha_; clear err_cum_alphax_;
  for j = 1:feat_dim
    [sort_sv,sort_idx] = sort(model.SVs(:,j));
    cum_alpha = cumsum(alpha(sort_idx));
    cum_alphax = cumsum(alpha(sort_idx).* sort_sv);
    err_cum_alpha_(j) = max(abs(cum_alpha'-exact_model.cum_alpha(j,2:end)));
    err_cum_alphax_(j) = max(abs(cum_alphax'-exact_model.cum_alphax(j,2:end)));
  end
  err_cum_alpha(i) = max(err_cum_alpha_);
  err_cum_alphax(i) = max(err_cum_alphax_);
  fprintf('\t-- Exact Model  (sort_sv:%f, cum_alpha:%f, cum_alphax:%f)\n',...
          err_sort_sv(i),err_cum_alpha(i),err_cum_alphax(i));
  
  %check approximate model too
  num_bins = max(1,round(rand(1,1)*100));
  pre_str = sprintf('-m 1 -n %i',num_bins);
  approx_model = precomp_model(model,pre_str);
  
  %check mins and maxs
  
  err_min(i) = max(abs(approx_model.min_sv - min(model.SVs)));
  err_max(i) = max(abs(approx_model.max_sv - max(model.SVs)));
  
  %check the values of h
  clear h;
  for j = 1:feat_dim
    svx=linspace(min(model.SVs(:,j)),max(model.SVs(:,j)),num_bins+1);
    for k=1:length(svx)
      svxk = svx(k)*ones(model.totalSV,1);
      mmx = min(svxk,model.SVs(:,j));
      h(j,k) = sum(alpha'*mmx); 
    end
  end
  % check interpolation coeffs
  a = num_bins./(max(model.SVs) - min(model.SVs));
  b = -min(model.SVs).*a;
  err_a(i) = max(abs(approx_model.a - a));
  err_b(i) = max(abs(approx_model.b - b));
  err_h(i) = max(max(abs(approx_model.h - h)));
  fprintf(['\t-- Approx Model (min_sv:%f, max_sv:%f, err_h:%f, ' ...
           'err_a:%f, err_b:%f)\n'],...
          err_min(i),err_max(i),err_h(i),err_a(i),err_b(i));
  
  %check the prediction values too
  test_data = rand(size(data));
  approx_model = precomp_model(model,'-m 1 -n 100');
  [l,a,libsvm_values] = svmpredict(labels,test_data,model);
  pwc_values = fiksvm_predict(labels,test_data,approx_model,'-e 0 -a 0');
  pwl_values = fiksvm_predict(labels,test_data,approx_model,'-e 0 -a 1');

  err_pwc(i) = median(abs(pwc_values-libsvm_values));
  err_pwl(i) = median(abs(pwl_values-libsvm_values));
  fprintf('\t-- Predictions (pwc_err:%f, pwl_err:%f, exact_err:%f)\n\n',...
           err_pwc(i),err_pwl(i),0);

end

%GLOBAL STATS
fprintf('=============Overall Statistics ===================\n\n');
fprintf('\t-- Exact Model  (sort_sv:%f, cum_alpha:%f, cum_alphax:%f)\n',...
        max(err_sort_sv),max(err_cum_alpha),max(err_cum_alphax));
fprintf(['\t-- Approx Model (min_sv:%f, max_sv:%f, err_h:%f, ' ...
         'err_a:%f, err_b:%f)\n'],...
        max(err_min),max(err_max),max(err_h),max(err_a),max(err_b));
  fprintf('\t-- Predictions (pwc_err:%f, pwl_err:%f, exact_err:%f)\n\n',...
           median(err_pwc),median(err_pwl),0);
