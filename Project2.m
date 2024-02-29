%%
clear all
load('iddata-06.mat');
u_id=id.InputData;
y_id=id.OutputData;
t=0:2:3998;
% plot(t,u_id);figure;plot(t,y_id)
% figure
% plot(t,val.y);figure;plot(t,val.u );

%%
  for na = 1:1
      nb = na;
   for m = 1:2


l = 1;
c = 1;
clear combinations 
combinations(l,:) = zeros(1,na+nb);
v = zeros(1,na+nb);
len = height(combinations);
while c<=(na+nb)
    for nr = 1:m
        v = zeros(1,na+nb);
        v(c) = nr;
        for k = 1:len
            Vprot = combinations(k,:)+v;
            sum = 0;
            for number = 1:(na+nb)
                sum = sum+Vprot(number);
            end
            flag = 0;
            for contr = 1:height(combinations)
                if combinations(contr,:) == Vprot
                    flag = 1;
                end
            end
            if (sum <= m) && (flag == 0)
                l = l+1;
                combinations(l,:) = Vprot;
            end
            Vprot = 0;
        end
        len = height(combinations);
    end
  c = c+1;
end

data_matrix_id = zeros(length(id.u),na+nb);

% data matrix on id
for k = 1:length(id.u)
    for j = 1:na
            if k-j<=0
                data_matrix_id(k,j) = 0;
            else
                data_matrix_id(k,j) = id.y(k-j);
            end
    end

    for j = 1:nb
            if k-j <=0
                 data_matrix_id(k,na+j) = 0;
            else
                 data_matrix_id(k,na+j) = id.u(k-j);
            end
    end
end

% data matrix pe val
data_matrix_val = zeros(length(val.u), na+nb);
for k = 1:length(val.u)
    for j = 1:na
            if k-j<=0
                data_matrix_val(k,j) = 0;
            else
                data_matrix_val(k,j) = val.y(k-j);
            end
    end

    for j = 1:nb
            if k-j <=0
                 data_matrix_val(k,na+j) = 0;
            else
                 data_matrix_val(k,na+j)= val.u(k-j);
            end
    end
end

lines_c = length(combinations(:,end));
% PHI, Y_hat and theta on id
PHI_id = ones(length(id.y),lines_c);

for i = 1:length(id.y) % nr of lines - one line corresponds to one polynomial equation
   for j = 1:lines_c   % for the number of columns of PHI (which is = with the number of lines of combinations) - one column corresponds to one product of the polinomial equation
      for k = 1:na+nb  % for taking the corresponding data from the data_matrix and the corresponding powers from the combinations matrix -each k corresponds to one term of the j'th product 
            PHI_id(i,j) = PHI_id(i,j)*data_matrix_id(i,k)^combinations(j,k); 
      end
   end
end

THETA  = PHI_id\id.y;

Y_hat_id = PHI_id*THETA; % y prediction for identification set

% PHI, Y_HAT ON VAL
PHI_val = ones(length(val.y),lines_c);

for i = 1:length(val.y) % nr of lines
   for j = 1:lines_c   % for the number of columns of PHI which is the number of lines of combinations
      for k = 1:na+nb % for taking the corresponding data from the data_matrix and the corresponding powers from the combinations matrix
            PHI_val(i,j) = PHI_val(i,j)*data_matrix_val(i,k)^combinations(j,k);
      end
   end
end


Y_hat_val = PHI_val*THETA; % y prediction for validation set


% Ysim id
Y_sim_id = zeros(length(id.y),1); % y simulation for identification

for k = 1:length(id.y)
for p = 1:lines_c
    prod = 1;
      for i = 1:na
           if(k-i<=0)
               prod = 0;
           else
               prod = prod * (Y_sim_id(k-i)^combinations(p,i));
           end
      end

      for i = 1:nb
            if(k-i<=0)
               prod = 0;
            else
              prod = prod * (id.u(k-i)^combinations(p,na+i));
            end
      end
  Y_sim_id(k) = Y_sim_id(k) + prod * THETA(p);
end

end

% Ysim val

Y_sim_val = zeros(length(val.y),1); % y simulation for identification

for k = 1:length(val.y)
for p = 1:lines_c
    prod = 1;
      for i = 1:na
           if(k-i<=0)
               prod = 0;
           else
               prod = prod * (Y_sim_val(k-i)^combinations(p,i));
           end
      end

      for i = 1:nb
            if(k-i<=0)
               prod = 0;
            else
              prod = prod * (val.u(k-i)^combinations(p,na+i));
            end
      end
  Y_sim_val(k) = Y_sim_val(k)+ prod * THETA(p);
end

end


 % errors
 % For identification
 s1 = 0;
 for i = 1:length(id.y)
     s1 = s1+(id.y(i)-Y_hat_id(i)).^2;
 end
 s1 = 1/length(id.y)*s1;
 V_MSE_approx_id(m,na)=s1;


 s2 = 0;
 for i = 1:length(id.y)
     s2 = s2+(id.y(i)-Y_sim_id(i)).^2;
 end
 s2 = 1/length(id.y)*s2;
 V_MSE_sim_id(m,na)=s2;


% For validation

 s3 = 0;
 for i = 1:length(val.y)
     s3 = s3+(val.y(i)-Y_hat_val(i)).^2;
 end
 s3 = 1/length(val.y)*s3;
 V_MSE_approx_val(m,na)=s3;


 s4 = 0;
 for i = 1:length(val.y)
     s4 = s4+(val.y(i)-Y_sim_val(i)).^2;
 end
 s4 = 1/length(val.y)*s4;
 V_MSE_sim_val(m,na)=s4;


    end
  end

 



%% plot the prediction and real data for identification
 figure
 plot(Y_hat_id,'r');
 hold on
 plot(id.y,'b');
 legend('approximation-Y.hat.id','system-id.y')
 title ('Approximation for identification set ')

 %% plot the prediction and real data for validation
 figure
 plot(Y_hat_val,'r');
 hold on
 plot(val.y,'b');
 legend('approximation-Y.hat.val','system-val.y');
 title ('Approximation for validation set ');

%% plot the simulation and real data for identification
 figure
 plot(Y_sim_id,'r');hold;
 plot(id.y,'b'); 
 legend('simulation-Y.sim.id','system-id.y');
 title ('Simulation for identification set');
 %% plot the simulation and real data for validation
 figure
 plot(Y_sim_val,'r');hold;
 plot(val.y,'b');
 legend('simulation-Y.sim.val','system-val.y')
 title ('Simulation for validation set ')