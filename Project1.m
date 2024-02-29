%%
%%we looked at the data that we where given
clear all
load('proj_fit_15.mat')
mesh(id.X{1},id.X{2},id.Y)

%%
m=5;

for e=1:30
%we calculate an array that contains the powers of each term
c=1;
for i=1:e+1
    for j=1:e+1
        if (i+j)<=(e+2)
p(c)=i-1;
p(c+1)=j-1;
c=c+2;
        end
    end
end
%we compute the matrix phi with the identification data
phi=ones(length(id.X{1})*length(id.X{2}),length(p)/2);
nr=1;
for i=1:length(id.X{1})
    for j=1:length(id.X{2})
        for k=2:2:length(p)
            phi(nr,k/2)=(id.X{1}(i).^p(k-1))*(id.X{2}(j).^p(k));
        end  
        nr=nr+1;
    end
end

%we transform the identification matrix into an array
nr2=1;
v_yid=ones(length(id.Y)*length(id.Y),1);
for i=1:length(id.Y)
    for j=1:length(id.Y)
   v_yid(nr2)=id.Y(i,j);
   nr2=nr2+1;
    end 
end


theta=phi\v_yid;

%we calculate the phi matrix with the validation data
phi1=ones(length(val.X{1})*length(val.X{2}),length(p)/2);
nr3=1;
for i=1:length(val.X{1})
    for j=1:length(val.X{2})
        for k=2:2:length(p)
            phi1(nr3,k/2)=(val.X{1}(i).^p(k-1))*(val.X{2}(j).^p(k));
        end  
        nr3=nr3+1;
    end
end

%we calculate the approximation for both phi and phi1
y_val=phi1*theta;
y_id=phi*theta;
%the approximation in matrix form
nr4=1;
y_hat=ones(length(val.Y));
for i=1:71
    for j=1:71
y_hat(i,j)=y_val(nr4);
nr4=nr4+1;
    end
end

nr5=1;
v_yval=ones(length(val.Y)*length(val.Y),1);
for i=1:length(val.Y)
    for j=1:length(val.Y)
   v_yval(nr5)=val.Y(i,j);
   nr5=nr5+1;
    end 
end
%we calculate the error arrays
E1(e)=1/(length(val.X{1})*length(val.X{2}))*sum((v_yval-y_val).^2);
E2(e)=1/(length(id.X{1})*length(id.X{2}))*sum((v_yid-y_id).^2);

%when we reach the optimal degree, we're plotting the approximation and the
%validation set
if e==m
subplot(121)
mesh(val.X{1},val.X{2},val.Y);title('System model');
subplot(122)
mesh(val.X{1},val.X{2},y_hat);
end

end

%calculate the MSE and m_optimal is the optimal degree
[MSE,m_optimal]=min(E1);title('Approximated model.MSE=',MSE)
%we plot the error arrays for later analysis
figure;
plot(E1);hold;plot(E2);