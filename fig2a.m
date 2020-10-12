tic

AT = [];
BT = [];

for FF=1:6
f = (FF-1)*0.02;
NUM_INF = [];
tend = 700;
sum_S = zeros(tend,n);

d = mean(sum(A,2));
X=sum(A,2);
CV = std(X)/mean(X);

ADDTRANS = [];
BOUNCE = [];

%f = 0.3;
HERD = [];
HIT = [];
for J = 1:3
num_infected = [];
Pi = 0.2; % delta (1/5 days)
Pr = 0.0667; % gamma
Pe = 3*Pr/d/(1+CV^2); % beta

S(1,:) = zeros(1,n);
H = zeros(1,tend);
init = 0;
while init < n/100
    ix = rand;
    if S(1,floor(ix*n)+1) ~=2
        S(1,floor(ix*n)+1) = 2; % Seed the initial infected
        init = init + 1;
    end
end
num_infected(1) = n/100;
allA{1} = A;
% States: 0, susceptible, 1, exposed, 2, infected, 3 recovered

this_A = A;
deg_dist = sum(A,2);
%sspreaders = find(deg_dist>10);

for t = 2:tend
    %S(t-1,sspreaders) = 3; % Superspreaders are removed
    exposed = [];
    infected = [];
    recovered = [];
    herd = 0;
    expect = [];
    for i = 1:n
        if S(t-1,i) == 0 % Susceptible
             if max(S(t-1,this_A(i,:)==1)==2)
                 herd = herd + 1;
             end
        end
            
        if S(t-1,i) == 1 % Exposed
            if rand < Pi
                infected = [infected i];
            end
        end

        if S(t-1,i) == 2 % Infected
            R0S = length(find(S(t-1,this_A(i,:)==1)==0))/sum(this_A(i,:));
            if isnan(R0S)
                R0S = 0;
            end
            expect = [expect R0S];
            if rand < Pr
                recovered = [recovered i];

            end
            for j = 1:n
                if this_A(i,j) && (S(t-1,j) == 0) % Connected and susceptible
                    if rand < Pe
                        exposed = [exposed j];
                    end
                end
            end
        end
        
        if t == 85
            if (S(t-1,i) == 0) && (rand<f)
                infected = [infected i];
            end
        end
    end
    
    S(t,:) = S(t-1,:);
    S(t,exposed) = 1;
    S(t,infected) = 2;
    S(t,recovered) = 3;
    H(t) = mean(expect); % fraction of at-risk individuals
    if isnan(H(t))
        H(t) = 0;
    end
    num_infected(t) = length(infected);
    allA{t} = this_A;
end

HERD(:,J) = H;
xx=find(H*3<1);
yy=cumsum(num_infected)/n;
HIT(J) = yy(xx(2));
%plot(cumsum(num_infected)/n)
NUM_INF(:,J) = cumsum(num_infected)/n;
sum_S = sum_S + S;
BOUNCE(J) = yy(85)-refI(85)-(yy(84)-refI(84));
ADDTRANS(J) = yy(700)-refI(700)-(yy(84)-refI(84));
end

AT(FF,:)=ADDTRANS;
BT(FF,:)=BOUNCE;
end

toc