X=sum(A,2);

%thresh = 10;
%phi = length(X(X>thresh))/length(X);
%X=X(X<=thresh);

CV = std(X)/mean(X);
R0 = Pe/Pr*mean(X);
F = (1+CV^2);
adjR0 = R0*F;

% mean degree, R0, mean saturation, mean herd immunity
[mean(X)  adjR0 mean(NUM_INF(end,:)) mean(HIT)]


zz=[0:0.02:0.1];
figure; hold on;
scatter(zz,mean(AT,2),'k')
scatter(zz,mean(BT,2),'b')
CI=[];
CIB = [];
for i = 1:length(zz)
    ci = bootci(10000,@mean,AT(i,:));
    CI = [CI ci];
    ci = bootci(10000,@mean,BT(i,:));
    CIB = [CIB ci];
end
errorbar(zz,mean(AT,2),mean(AT,2)-CI(1,:)',CI(2,:)'-mean(AT,2),'k', 'LineStyle', 'none')
errorbar(zz,mean(BT,2),mean(BT,2)-CIB(1,:)',CIB(2,:)'-mean(BT,2),'b', 'LineStyle', 'none')
%xlim([0 1])
%ylim([0 1])
f=fit(zz',mean(BT,2),'poly1');
fitlm(mean(BT,2),mean(AT,2))
%plot(zz,f(zz),'b')
%plot(zz,1.5*f(zz),'k')
box on
set(gcf, 'Position',  [  393   375   802/2   496/2])
mean(AT,2)./mean(BT,2)

set(gcf, 'Position',  [  493   450   164   124])

[0 0;0.0085  0.0068; 0.0329 0.0305; 0.0602 0.0602]

figure; hold on;
xxx = [mean(A1) mean(A1p) mean(A2) mean(A2p) mean(A3)];
bar([1 2 3 4 5],xxx,0.5,'g')
c1 = bootci(10000,@mean,A1);
c2 = bootci(10000,@mean,A1p);
c3 = bootci(10000,@mean,A2);
c4 = bootci(10000,@mean,A2p);
c5 = bootci(10000,@mean,A3);
errorbar([1 2 3 4 5],xxx,xxx-[c1(1) c2(1) c3(1) c4(1) c5(1)],[c1(2) c2(2) c3(2) c4(2) c5(2)]-xxx,'k', 'LineStyle', 'none')
ylim([0 1])



HT = vertcat(zeros(1,2),HT);

zz=[1.0:0.5:3.0];
figure; hold on;
scatter(zz,mean(HT,2),'k')
CI=[];
for i = 1:length(zz)
    ci = bootci(10000,@mean,HT(i,:));
    CI = [CI ci];
end
errorbar(zz,mean(HT,2),mean(HT,2)-CI(1,:)',CI(2,:)'-mean(HT,2),'k', 'LineStyle', 'none')
zz=[1.0:0.01:3.0];
plot(zz,1-1./zz)



figure; hold on; 
plot(diff(NUM_INF,1))
plot(diff(mean(NUM_INF,2),1),'g','LineWidth',5)
ylim([0 0.02])
set(gcf, 'Position',  [  19   762   961   229])

min(find(mean(NUM_INF,2)>mean(HIT)))

clear all
n = 100;
A(1:n,1:n) = 1;
for i = 1:n
    A(i,i) = 0;
end

figure;plot(mean(NUM_INF,2))

%if t == 250 || t == 400 || t == 550 || t == 700 || t == 850 || t == 1000
[mean((NUM_INF(399,:)-NUM_INF(249,:))/(NUM_INF(250,:)-NUM_INF(249,:))) mean((NUM_INF(549,:)-NUM_INF(399,:))/(NUM_INF(400,:)-NUM_INF(399,:))) mean((NUM_INF(699,:)-NUM_INF(549,:))/(NUM_INF(550,:)-NUM_INF(549,:))) mean((NUM_INF(849,:)-NUM_INF(699,:))/(NUM_INF(700,:)-NUM_INF(699,:))) mean((NUM_INF(999,:)-NUM_INF(849,:))/(NUM_INF(850,:)-NUM_INF(849,:))) mean((NUM_INF(1150,:)-NUM_INF(999,:))/(NUM_INF(1000,:)-NUM_INF(999,:)))]

y = totalI-totalIbase;
z = totalI-totalIbase;
[(y(219)-y(69))/(z(70)-z(69)) (y(369)-y(219))/(z(220)-z(219)) (y(519)-y(369))/(z(370)-z(369)) (y(669)-y(519))/(z(520)-z(519)) (y(819)-y(669))/(z(670)-z(669)) (y(1150)-y(819))/(z(1150)-z(819))]


% t == 55 || t == 205 || t == 355 || t == 505 || t == 655 || t == 805
y = mean(NN,2);
z = mean(NUM_INF,2);
id = [55,205,355,505,655,805,1150];
for i = 1:length(id)-1
    ((y(id(i+1)-1)-y(id(i)-1))-(z(id(i+1)-1)-z(id(i)-1)))/((y(id(i))-y(id(i)-1))-(z(id(i))-z(id(i)-1)))
end

% t == 85 || t == 235 || t == 385 || t == 535 || t == 685 || t == 835
y = mean(NN,2);
z = mean(NUM_INF,2);
id = [85,235,385,535,685,835,1150];
for i = 1:length(id)-1
    ((y(id(i+1)-1)-y(id(i)-1))-(z(id(i+1)-1)-z(id(i)-1)))/((y(id(i))-y(id(i)-1))-(z(id(i))-z(id(i)-1)))
end

% t == 70 || t == 220 || t == 370 || t == 520 || t == 670 || t == 820
y = totalI;
z = reftotal;
id = [70,220,370,520,670,820,1150];
for i = 1:length(id)-1
    ((y(id(i+1)-1)-y(id(i)-1))-(z(id(i+1)-1)-z(id(i)-1)))/((y(id(i))-y(id(i)-1))-(z(id(i))-z(id(i)-1)))
end



% WS final 85
AT=[ 0.0241    0.0570 -0.0013    0.0819 0.0168    0.0722  0.0114    0.0035    0.0034    0.0699    0.0353 0.0817    0.0475    0.0074    0.0316    0.0126    0.0270    0.0078    0.0741    0.0441];
BT=[ 0.0427    0.0432  0.0375    0.0468 0.0367    0.0437 0.0370    0.0391    0.0407    0.0443    0.0374 0.0462    0.0390    0.0349    0.0402    0.0436    0.0417    0.0373    0.0474    0.0430];

AT = [
       -0.0158   -0.0081          0
   -0.0151   -0.0071    -0.0242
    0.0080    0.0109     0.0218
    0.0089   -0.0236    -0.0010
    0.0406    0.0019    -0.0196
   -0.0013    0.0819 0.0241 ];
BT = [    0.0021   -0.0002          0
    0.0072    0.0086     0.0076
    0.0176    0.0164     0.0201
    0.0247    0.0214     0.0223
    0.0316    0.0286     0.0289
    0.0375    0.0468  0.0427];

