function [K1,K2,K3]=find_multipli(NN,K_MIN,K_MAX);
K1=0;K2=0;K3=0;
for N=NN+20:-1:NN
K = 1:ceil(sqrt(N));D = K(rem(N,K)==0);D = [D sort(N./D)];
ii=find(D>K_MIN&D<K_MAX);
  if(numel(ii)>0)
  fprintf('N:%5.0f ',N);
  for jj=ii, 
       fprintf('%6.0f(%6.0f) ',[N./D(jj) ; D(jj)]);
       K1=D(jj);
       K2=N./D(jj);
       K3=N;
  end 
  fprintf('\n');
  end
end


end
