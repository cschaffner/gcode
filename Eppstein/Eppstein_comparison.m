%% performs Eppstein "inverse optimization" for 20 amino properties

% created: June 28, 2011
% by Christian Schaffner, c.schaffner@uva.nl


%% clear workspace
clear;

%% reading in things
geneticcode;
% 
% polar=[ 50 49 49 53 56 75 66 66 70 54 84 86 100 101 130 125 48 52 91 79]/10;
% polarvals=polar;
% 
% % theoretical PR (polar requirement)
% theoreticalPR = [ 4.5 4.4 5.0 5.0 6.2 7.5 6.1 6.2 6.5 7.7 7.9 8.9 9.6 10.2 12.2 13.6 4.3 4.9 8.6 9.0];

% sets of 20 numerical values for amino-acid properties which make the SGC most error robust

% using the all-1 weights
w1=[
0.571	-1	-1.11	-0.826
0.473	-0.55	-0.58	-0.583
0.156	-0.73	-0.97	-0.597
0	-1.62	-1.43	-0.742
0.193	-0.07	-0.17	-0.593
0.57	-0.49	-0.43	-0.636
0.489	-0.19	-0.01	-0.453
0.212	-0.31	-0.23	-0.504
0.233	0.28	0.4	-0.561
1.131	0.57	0.46	-0.477
0.498	1.08	0.92	1.014
0.491	1.29	1.21	2.058
0.178	0.72	0.74	0.904
0.078	0.86	0.83	1.99
0.188	1.61	1.64	0.897
0.124	2.04	1.76	1.629
1.501	-1.36	-1.29	-0.848
1.834	-1.83	-1.62	-0.855
0.568	-0.35	-0.28	-0.292
0.511	0.04	0.15	-0.525 ];

% using FH weights
FHw=[-1.46	-1.33	-1.26	0.24	0.25
-1.1	-1.25	-1.07	0.18	0.19
-1.34	-1.28	-1.37	0.01	0.01
-1.49	-1.36	-1.4	0	0
-0.68	-0.63	-1.12	0.06	0.08
-0.3	-0.33	-0.14	0.37	0.39
-0.47	-0.49	-0.17	0.34	0.33
-0.65	-0.57	-0.62	0.26	0.26
-0.39	-0.35	-0.23	0.31	0.31
0.3	1.11	0.89	0.81	0.81
0.33	1.65	2.04	0.79	0.79
0.57	1.98	2.1	0.72	0.72
1.34	0.7	0.39	0.56	0.57
1.67	0.81	0.31	0.46	0.44
1.43	0.89	0.87	0.7	0.7
1.76	1.03	0.45	0.59	0.59
-0.03	-0.16	0	0.83	0.82
0.18	-0.22	-0.09	1.53	1.5
0.09	-0.08	0.28	0.64	0.64
0.24	-0.14	0.17	0.6	0.6];

%% comparing things

count=size(aaindex1,1);
A=vertcat(aaindex1,w1',FHw');

% compute correlations
C=corrcoef(A');

%interesting part, without self-correlations
CC=C(count+1:count+9,1:count);
for i=1:9
    [rows,cols]=find(abs(C(count+i,1:count)>0.8));
    if ~isempty(cols)
        for j=1:size(cols,2)
 %    fprintf('%i correlates well with %i (coeff: %f)\n',i,char(name(cols(j))),C(count+1,cols(j)))
            fprintf('%i correlates well with %i (coeff: %f)\n',i,cols(j),C(count+i,cols(j)))
        end
        fprintf('\n');
    end
end


break;


%% plotting things
M=vertcat(aa_polar,aa_theoreticalPR,w1',FHw');
M=normalize(M);

nmin=1;
nmax=10;
step=1;
figure;
hold on;
ycount=0;
for i=nmin:step:nmax
    ycount=ycount+1;
    plot(M(i,:),ones(1,20)*ycount,'x')
    % we want the aminonames flipping up and down
    [~,index]=sort(M(i,:));
    iindex=[index ; (-1) .^[1:20] ];
    iiindex=sortrows(iindex')';
    text(M(i,:)-0.01,ones(1,20)*ycount  - iiindex(2,:)*0.2,aminoshort);
%    text(M(i,:)-0.01,ones(1,20)*ycount - ones(1,20)*0.05 - iiindex(2,:)*0.1,aminoshort);
end
axis([-2 3 0.5 10.5]);

hold off;



