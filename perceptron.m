%input vectors
n=100;
x_mu=[0,0;0,1;1,0;1,1];
y_t=[0,1,1,0];
w= zeros(1,n);

p=50;
alpha=1;
b=1;
T = 50;
E_list=zeros(1,10);
for i = 1:T
    E=0;
    for mu = 1:p
        chi=w*x_mu(mu,:)'-b;
        if chi<0
            y_mu=0;
        else 
            y_mu=1;
        end
        y_d=y_t(mu)-y_mu;
        w=w+alpha*x_mu(mu,:)*(y_d);
        E=E+sum(y_d^2);  
    end
    E_list(i)=E;
end

figure;
plot(E_list)
xlabel('iteration number','fontsize',20)
ylabel('error','fontsize',20)
