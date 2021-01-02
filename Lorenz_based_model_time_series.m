

tau=10;
emb_dim=5;

tspan=[0 50000];
y0=[0; 1; 1; 0; 1.1; 1.1];

[t,y]=ode45(@f,tspan,y0);

[r_y r_t]=resample(y,t,15.,'spline');

x1=r_y(:,1);
x2=r_y(:,4);


org_x1=x1;
org_x2=x2;

R=corrcoef(x1,x2);

%con_later_w=(con_later_w_i-1)*0.01;
con_later_w=0.15;
%con_later_w=(con_later_w_i-1)*0.1;
w11=0.3;
w22=0.3;

% w12=0.1;
% w21=0.1;

w12=con_later_w;
w21=con_later_w;

b1=1.8+3.0;
b2=1.8+3.0;

%b1=0+(b_i-1)*0.2;
%b2=0+(b_i-1)*0.2;

%para_b(b_i)=b1;

s1=0.3;
s2=0.3;

beta1=2.0;
beta2=2.0;

LC_ac=1.5;

theta1=0.;
theta2=0.;

base_diami=3;


x1=LC_ac*zscore(x1);
x2=LC_ac*zscore(x2);



x1=x1+b1;
x2=x2+b2;

lc_x1=x1;
lc_x2=x2;


sp_m1=activate_func_ew(-w11*x1-w21*x2+beta1-theta1);
sp_m2=activate_func_ew(-w22*x2-w12*x1+beta2-theta2);

di_m1=s1*x1;
di_m2=s2*x2;

pup1=di_m1-sp_m1+base_diami;

pup2=di_m2-sp_m2+base_diami;


set(0,'defaultAxesFontSize',20);
set(0,'defaultAxesFontName','Arial');
set(0,'defaultTextFontSize',20);
set(0,'defaultTextFontName','Arial');
set(gca,'color',[1 1 1])

figure;

plot(r_t,org_x1,'-b','Linewidth',2);
hold on;

grid on;

plot(r_t,org_x2,'-r','Linewidth',2);

xlim([200 300]);
xlabel('time: t');
ylabel('Variable X_i of Lorenz model');

legend('i=1 (left)','i=2 (right)');

figure;

plot(r_t,lc_x1,'-b','Linewidth',2);
hold on;

grid on;

plot(r_t,lc_x2,'-r','Linewidth',2);

xlim([200 300]);
xlabel('time: t');
ylabel('LC activity: x_i');

legend('i=1 (left)','i=2 (right)');






figure;

plot(r_t,sp_m1,'-b','Linewidth',2);
hold on;

grid on;

plot(r_t,sp_m2,'-r','Linewidth',2);


xlim([200 300]);
xlabel('time: t');
ylabel('Output EWN to ciliary ganglion: S_i');

legend('i=1 (left)','i=2 (right)');



figure;

plot(r_t,di_m1,'-b','Linewidth',2);
hold on;

grid on;

plot(r_t,di_m2,'-r','Linewidth',2);

xlim([200 300]);

xlabel('time: t');
ylabel('Inputs to the superior cervical ganglion: D_i');

legend('i=1 (left)','i=2 (right)');



figure;

plot(r_t,pup1,'-b','Linewidth',2);
hold on;

grid on;

plot(r_t,pup2,'-r','Linewidth',2);

xlim([200 300]);

xlabel('time: t');
ylabel('Pupil diameter: P_i');

legend('i=1 (left)','i=2 (right)');



function y=activate_func_ew(x)
   
    y=tanh(x)+1;
    
end

function dydt=f(t,y)
a1=10.;
a2=10.;
b1=8./3.;
b2=8./3.;
c1=28.;
c2=28.;
J=0.7;
%J=0.1;

dydt =[
    a1*(y(2)-y(1))+J*(y(4)-y(1));
    c1*y(1)-y(1)*y(3)-y(2);
    y(1)*y(2)-b1*y(3);
    a2*(y(5)-y(4))+J*(y(1)-y(4));
    c2*y(4)-y(4)*y(6)-y(5);
    y(4)*y(5)-b2*y(6);    
    ];

end

