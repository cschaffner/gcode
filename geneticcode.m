%% Definitions of genetic code matrices and amino acid values
%
% this MATLAB module does the following:
% * clear the whole workspace
% * defines short names of amino acids
% * reads in various values of amino acid properties
% * computes corresponding A matrices
% * defines various A matrices for pairwise differences between amino acids
% * defines permutations for some optimal codes
% * defines the block structure of the standard genetic code
% * creates various B matrices (containing the number of neighbors reachable
% by 1 mutation)

% created: Feb 3, 2011
% cleanup: Oct 21, 2011
% by Christian Schaffner, c.schaffner@uva.nl

%%
% clear workspace
clear all;

% declare some variable global
global Code ASerge BlockSize probSerge aminos aa_theoreticalPR;

% suppression switch
suppression = 0;

%% names of amino acids
aminos = {'Phe','Leu','Ile','Met','Val','Ser','Pro','Thr','Ala','Tyr', 'His', 'Gln', 'Asn', 'Lys', 'Asp', 'Glu', 'Cys', 'Trp', 'Arg', 'Gly', 'STOP'};

%% various amino acid values (variable names starting with aa_)
% from the original Haig and Hurst
aa_polar=[ 50 49 49 53 56 75 66 66 70 54 84 86 100 101 130 125 48 52 91 79]/10;

%polarvals=polar;

% theoretical polar requirement
aa_theoreticalPR = [ 4.5 4.4 5.0 5.0 6.2 7.5 6.1 6.2 6.5 7.7 7.9 8.9 9.6 10.2 12.2 13.6 4.3 4.9 8.6 9.0];

% first column:  Grantham's volume                         
% second column: Hydropathy
% third column:  Isoelectric point
aa_otherValues=[ 132                                               2.8                                          5.48
111                                              3.8                                          5.98
111                                              4.5                                          6.02
105                                              1.9                                          5.74
84                                               4.2                                          5.96
32                                             -0.8                                          5.68
32.5                                          -1.6                                          6.30
61                                            -0.7                                           6.16
31                                               1.8                                          6.00
136                                              -1.3                                         5.66
96                                              -3.2                                         7.59
85                                              -3.5                                        5.65
56                                              -3.5                                        5.41
119                                              -3.9                                        9.74
54                                              -3.5                                        2.77
83                                               -3.5                                        3.22
55                                                2.5                                        5.07
170                                              -0.9                                        5.89
124                                              -4.5                                       10.76
3                                               -0.4                                         5.97 ];

aa_Grantham = aa_otherValues(:,1)';
aa_hydropathy = aa_otherValues(:,2)';
aa_isoelectric = aa_otherValues(:,3)';

% Apolar=zeros(20,20);
% for i=1:20
%     for j=1:20
%         Apolar(i,j)=abs(polar(i)-polar(j))^2;
%     end
% end
% short form:

% define A matrices based on aa_ values
p=2;
Apolar = mypdist(aa_polar') .^ p;
Atheoreticpolar = mypdist(aa_theoreticalPR') .^ p;
AGrantham = mypdist(aa_Grantham') .^ p;
Ahydropathy = mypdist(aa_hydropathy') .^ p;
Aisoelectric = mypdist(aa_isoelectric') .^ p;

% certified numbers:
% all-1 weights
aa_w1=[
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

% FH weights
aa_FH=[-1.46	-1.33	-1.26	0.24	0.25
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


% Wimley & White: nature structural biology, 1996
aa_WWValues= [538 481 452 448 418 412 380 411 408 519 408 367 383 377 302 223 449 610 391 424];

AWW=mypdist(aa_WWValues') .^ 2;

% reading in aaindex1
count=0;
fid = fopen('aaindex1');
% initialize array
aaindex1=single(zeros(544,20));

while ~feof(fid)
    C = textscan(fid, 'H %s\n D %[^\n]');
    if isempty(C{1}) 
           break ;
    end
    D = textscan(fid,'%c %*[^\n]',1);
    while ~(D{1} == 'I' && strcmp(D{2},'A/L'))
        D = textscan(fid, '%c %s %*[^\n]',1);
    end
    F= textscan(fid, '%f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 \n',2,'TreatAsEmpty','NA');

    count=count+1;

    % assign values
    aaind_name(count)=C{1};
    aaind_desc(count)=C{2};
    aaindex1(count,1)=F{4}(2);
    aaindex1(count,2)=F{1}(2);
    aaindex1(count,3)=F{10}(1);
    aaindex1(count,4)=F{3}(2);
    aaindex1(count,5)=F{10}(2);
    aaindex1(count,6)=F{6}(2);
    aaindex1(count,7)=F{5}(2);
    aaindex1(count,8)=F{7}(2);
    aaindex1(count,9)=F{1}(1);
    aaindex1(count,10)=F{9}(2);
    aaindex1(count,11)=F{9}(1);
    aaindex1(count,12)=F{6}(1);
    aaindex1(count,13)=F{3}(1);
    aaindex1(count,14)=F{2}(2);
    aaindex1(count,15)=F{4}(1);
    aaindex1(count,16)=F{7}(1);
    aaindex1(count,17)=F{5}(1);
    aaindex1(count,18)=F{8}(2);
    aaindex1(count,19)=F{2}(1);
    aaindex1(count,20)=F{8}(1);

    textscan(fid,'////');
end
fclose(fid);
clear C D F count;


%% define other A matrices:
% from Benner, PAM 74-100
Benner=[     2.4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    -0.8     4.8  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    -0.2     0.3     3.6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    -0.3    -0.5     2.2     4.8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
     0.3    -2.2    -1.8    -3.2    11.8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    -0.3     1.6     0.7     0.8    -2.6     3.0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    -0.1     0.3     1.0     2.9    -3.2     1.7     3.7 0 0 0 0 0 0 0 0 0 0 0 0 0 
     0.6    -1.0     0.4     0.2    -2.0    -1.1    -0.5     6.6 0 0 0 0 0 0 0 0 0 0 0 0 
    -1.0     1.0     1.2     0.4    -1.3     1.4     0.2    -1.6     6.1 0 0 0 0 0 0 0 0 0 0 0 
    -0.8    -2.6    -2.8    -3.9    -1.2    -2.0    -2.9    -4.3    -2.3     4.0 0 0 0 0 0 0 0 0 0 0
    -1.4    -2.4    -3.1    -4.2    -1.6    -1.7    -3.1    -4.6    -1.9     2.8     4.2 0 0 0 0 0 0 0 0 0
    -0.4     2.9     0.9     0.4    -2.9     1.7     1.2    -1.1     0.6    -2.3    -2.4     3.4 0 0 0 0 0 0 0 0
    -0.8    -1.8    -2.2    -3.2    -1.2    -1.0    -2.2    -3.5    -1.5     2.6     2.9    -1.5     4.5 0 0 0 0 0 0 0 
    -2.6    -3.5    -3.2    -4.7    -0.7    -2.8    -4.3    -5.4     0.0     0.9     2.1    -3.6     1.3     7.2 0 0 0 0 0 0
     0.4    -1.0    -1.0    -1.0    -3.1    -0.2    -0.7    -1.7    -1.0    -2.6    -2.2    -0.8    -2.4    -3.8     7.5 0 0 0 0 0
     1.1    -0.2     0.9     0.4     0.1     0.1     0.1     0.4    -0.3    -1.8    -2.2     0.0    -1.4    -2.6     0.5     2.1 0 0 0 0
     0.7    -0.3     0.4    -0.2    -0.6    -0.1    -0.2    -1.0    -0.5    -0.3    -1.1     0.1    -0.4    -2.2     0.1     1.4     2.5 0 0 0
    -4.1    -1.6    -4.0    -5.5    -0.9    -2.8    -4.7    -4.1    -1.0    -2.3    -0.9    -3.6    -1.3     3.0    -5.2    -3.4    -3.7    14.7 0 0
    -2.6    -2.0    -1.4    -2.8    -0.4    -1.8    -3.0    -4.3     2.5    -1.0    -0.1    -2.4    -0.5     5.3    -3.4    -1.9    -2.1     3.6     8.1 0
     0.1    -2.2    -2.2    -2.9    -0.2    -1.7    -2.1    -3.1    -2.1     3.2     1.9    -1.9     1.8     0.1    -1.9    -1.0     0.2    -2.9    -1.4     3.4];
 
 Benner=Benner+Benner'-diag(diag(Benner));
 
% change the ordering to how we need it
 p_Benner=[14 11 10 13 20 16 15 17 1 19 9 6 3 12 4 7 5 18 2 8];
 Benner=Benner(p_Benner,p_Benner);

 ABenner= (10 .^ (-Benner/10) )^1;
 ABenner=ABenner-diag(diag(ABenner));

 clear p_Benner Benner;
 

% Serge Massar's matrix:
Serge=[ +7 0 0 0 0 0 0 0 0 0             0 0 0 0 0 0 0 0 0 0 
 -3 +7 0 0 0 0 0 0 0 0             0 0 0 0 0 0 0 0 0 0
 0 -4 +7 0 0 0 0 0 0 0            0 0 0 0 0 0 0 0 0 0
 0 -5 +2 +7 0 0 0 0 0 0           0 0 0 0 0 0 0 0 0 0
 0 -2 -2 -3 +7 0 0 0 0 0          0 0 0 0 0 0 0 0 0 0
 -2 -4 -1 -2 -3 +7 0 0 0 0        0 0 0 0 0 0 0 0 0 0
 +1 -2 +1 +1 -1 0 +7  0 0 0       0 0 0 0 0 0 0 0 0 0
 -2 -3 -4 -4 0 -4 -2 +7 0 0       0 0 0 0 0 0 0 0 0 0
 +1 -5 +2 +3 -2 -1 +1 -3 +7 0     0 0 0 0 0 0 0 0 0 0
 -1 -3 -3 -3 0 -4 -2 0 -2 +7      0 0 0 0 0 0 0 0 0 0
 +1 -3 -1 -1 0 -2 0 0 -1 +1 +7      0 0 0 0 0 0 0 0 0 
  0 -3 +2 +2 -2 0 +2 -3 +2 -2 -1 +7   0 0 0 0 0 0 0 0 
 -3 -6 -1 -2 -4 -4 -2 -5 -1 -5 -4 -1 +7 0 0 0 0 0 0 0
 +1 -4 +2 +3 -1 -1 +2 -2 +3 -1 0 +2 -2 +7 0  0 0 0 0 0 
 +1 -4 +2 +2 -2 -1 +1 -3 +3 -2 0 +2 -2 +2 +7 0 0 0 0 0 
 +1 -2 +2 +1 -1 0 +2 -2 +2 -1 0 +2 -1 +2 +2 +7 0 0 0 0 
 0 -2 +1 0 0 -1 +2 -1 +1 -1 0 +1 -1 +1 +1 +2 +7  0 0 0 
 -1 -3 -3 -4 0 -4 -2 +1 -3 0 0 -3 -4 -2 -2 -1 0 +7 0 0
 +1 -2 0 -1   0 -2 +1 -1 0 -1 +1 0 -3 0 0 +1 0 -1 +7 0
 0  -1 -1 -2 +1 -1 +1 -1 -1 -1 0 0 -3 0 -1 +1 +1 0 +1 +7 ];

Serge=Serge+Serge'-diag(diag(Serge));
p_Serge=[5 10 8 11 18 16 13 17 1 20 7 14 12 9 3 4 2 19 15 6];
Serge=Serge(p_Serge,p_Serge);

ASerge=-Serge;
%ASerge= (10 .^ (-Serge/10) )^1;
%ASerge=ASerge-diag(diag(ASerge));
clear p_Serge Serge;

% probabilities also defined in Serge's paper
probSerge=[ 7.80
 5.23
 5.19
 4.37
 1.10 
 6.72
 3.45 
 6.77
 2.03
 6.95
 10.15
 6.32
 2.28
 4.39
 4.26
 6.46
 5.12 
 1.09
 3.30 
 7.01 ];

pord2=[14 11 10 13 20 16 15 17 1 19 9 7 4 12 3 6 5 18 2 8];
probSerge=probSerge(pord2);
clear pord2;

% Peter's substitution matrix
Peter = [
0 , 0 , 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0
15, 0 , 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0 
21, 10, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0
21, 14, 14, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0
22, 15, 5 , 11, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0
17, 12, 14, 10, 11, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0
17, 8 , 8 , 10, 11, 10, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0
20, 13, 9 , 9,  6,  5,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0
16, 11, 13, 9,  10, 3,  9,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0
3 , 16, 22, 22, 23, 18, 18, 21, 17, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0
18, 15, 17, 17, 18, 13, 13, 16, 12, 19, 0,  0,  0,  0,  0,  0,  0,  0,  0, 0
20, 13, 13, 11, 12, 11, 9,  10, 10, 21, 12, 0,  0,  0,  0,  0,  0,  0,  0, 0
19, 14, 16, 12, 13, 8,  12, 11, 7,  20, 13, 13, 0,  0,  0,  0,  0,  0,  0, 0
17, 12, 12, 14, 15, 14, 8,  13, 13, 19, 17, 13, 16, 0,  0,  0,  0,  0,  0, 0
18, 13, 15, 11, 12, 7,  11, 10, 6,  19, 14, 12, 5,  15, 0,  0,  0,  0,  0, 0
19, 12, 12, 10, 11, 10, 8,  9,  9,  20, 15, 5,  12, 12, 11, 0  ,0,  0,  0, 0
17, 12, 14, 10, 11, 4,  10, 9,  3,  18, 13, 11, 8,  14, 7,  10, 0,  0,  0, 0
10, 23, 27, 27, 28, 23, 23, 26, 22, 13, 18, 22, 25, 23, 24, 25, 23, 0 , 0, 0
24, 15, 15, 17, 18, 17, 11, 16, 16, 25, 10, 12, 19, 15, 18, 15, 17, 24, 0, 0
19, 14, 14, 12, 11, 6,  12, 9,  5,  20, 15, 13, 10, 16, 9,  12, 6,  25,19, 0
];

Peter=Peter+Peter'-diag(diag(Peter)); % diagonal is already 0

APeter= (10 .^ (Peter/100) )^1;
APeter=APeter-diag(diag(APeter));
MatrixPeter=Peter;

%APeter = Peter;
clear Peter;

% Higgs matrix:
Higgs=[
 0    38    39    30    56   134   111    97   117    61    98   120   137   146   160   144    81    59   131   172
38     0    12    46    20   123    98    85    98    82   107   124   135   152   158   148    66    90   147   159
39    12     0    50    22   131   109    95   104    90   115   134   145   161   167   157    68    94   155   165
30    46    50     0    59   108    91    75    94    52    76    99   114   130   139   126    63    69   116   145
56    20    22    59     0   119    98    84    90    96   113   129   136   158   156   149    61   107   156   152
134   123   131   108   119     0    57    48    41   116    84    77    57   112    80    91    91   157   126    54
111    98   109    91    98    57     0    23    68    85    71    56    58    89    89    86    96   123   105   108
97    85    95    75    84    48    23     0    51    79    62    60    60    98    90    88    75   116   107    97
117    98   104    94    90    41    68    51     0   119    97   102    88   136   107   116    61   153   146    64
61    82    90    52    96   116    85    79   119     0    60    73   100    98   129   107   105    44    80   162
98   107   115    76   113    84    71    62    97    60     0    47    59    79    87    72   101    93    64   124
120   124   134    99   129    77    56    60   102    73    47     0    35    49    66    48   124   112    61   126 
137   135   145   114   136    57    58    60    88   100    59    35     0    68    40    40   122   139    86   100
146   152   161   130   158   112    89    98   136    98    79    49    68     0    89    67   163   130    46   158
160   158   167   139   156    80    89    90   107   129    87    66    40    89     0    32   144   163   110   111
144   148   157   126   149    91    86    88   116   107    72    48    40    67    32     0   144   140    85   131
  81    66    68    63    61    91    96    75    61   105   101   124   122   163   144   144     0   126   157   109
   59    90    94    69   107   157   123   116   153    44    93   112   139   130   163   140   126     0   104   200
  131   147   155   116   156   126   105   107   146    80    64    61    86    46   110    85   157   104     0   169
  172   159   165   145   152    54   108    97    64   162   124   126   100   158   111   131   109   200   169     0 ];

AHiggs=Higgs;
clear Higgs;

combHiggs= [ 
135 19.80 0.35 5.48 2.8 3.7 218 0.88 5.0 4.39
124 21.40 0.13 5.98 3.8 2.8 180 0.85 4.9 10.15
124 21.40 0.13 6.02 4.5 3.1 182 0.88 4.9 6.95
124 16.25 1.43 5.74 1.9 3.4 204 0.85 5.3 2.28
105 21.57 0.13 5.96 4.2 2.6 160 0.86 5.6 7.01
73 9.47 1.67 5.68 -0.8 0.6 122 0.66 7.5 6.46
90 17.43 1.58 6.30 -1.6 -0.2 143 0.64 6.6 4.26
93 15.77 1.66 6.16 -0.7 1.2 146 0.70 6.6 5.12
67 11.50 0.00 6.00 1.8 1.6 113 0.74 7.0 7.80
141 18.03 1.61 5.66 -1.3 -0.7 229 0.76 5.4 3.30
118 13.69 51.60 7.59 -3.2 -3.0 194 0.78 8.4 2.03
114 14.45 3.53 5.65 -3.5 -4.1 189 0.62 8.6 3.45
96 12.28 3.38 5.41 -3.5 -4.8 158 0.63 10.0 4.37
135 15.71 49.50 9.74 -3.9 -8.8 211 0.52 10.1 6.32
91 11.68 49.70 2.77 -3.5 -9.2 151 0.62 13.0 5.19
109 13.57 49.90 3.22 -3.5 -8.2 183 0.62 12.5 6.72
86 13.46 1.48 5.07 2.5 2.0 140 0.91 4.8 1.10
163 21.67 2.10 5.89 -0.9 1.9 259 0.85 5.2 1.09
148 14.28 52.00 10.76 -4.5 -12.3 241 0.64 9.1 5.23
48 3.40 0.00 5.97 -0.4 1.0 85 0.72 7.9 6.77 ];

% normalize columns
for i=1:9
    combHiggs(:,i) = (combHiggs(:,i) - mean(combHiggs(:,i)) ) / std(combHiggs(:,i));
end

combweightsHiggs = [0.000 0.155 0.000 0.028 0.218 0.000 0.277 0.179 0.142 ];

AHiggs2=mypdistweights(combHiggs(:,1:9),combweightsHiggs);
clear combHiggs combweightsHiggs;

% Higgs = round(Higgs2 / (sum(sum(Higgs2))/sum(sum(Higgs)) ))


A= Apolar;


%% optimal codes 

% Goldman's minimum code (GMC)
% optimal for Apolar, weights [1 1 1 1 1 1]
%pGoldman = [13 12 19 16 11 6 7 8 5 10 1 3 18 4 2 17 14 15 20 9 ];
pGoldman=[11 15 12 14 9 6 7 8 20 10 5 2 1 17 18 4 16 13 3 19];

% optimal for Apolar, FH weights [1 0.5 0.5 0.1 1 1]
pFHpolar=[15 11 12 13 18 6 17 20 19 14 9 8 2 5 4 3 10 16 1 7];

%% Block structure of genetic code

% U=0, C=1, A=2, G=3

Block{1}= [0 0 0 ; 0 0 1];
Block{2}= [0 0 2 ; 0 0 3 ; 1 0 0 ; 1 0 1 ; 1 0 2 ; 1 0 3];
Block{3}= [2 0 0 ; 2 0 1 ; 2 0 2 ];
Block{4}= [2 0 3 ];
Block{5}= [3 0 0 ; 3 0 1 ; 3 0 2 ; 3 0 3];
Block{6}= [0 1 0 ; 0 1 1 ; 0 1 2 ; 0 1 3 ; 2 3 0 ; 2 3 1];
Block{7}= [1 1 0 ; 1 1 1 ; 1 1 2 ; 1 1 3];
Block{8}= [2 1 0 ; 2 1 1 ; 2 1 2 ; 2 1 3];
Block{9}= [3 1 0 ; 3 1 1 ; 3 1 2 ; 3 1 3];
Block{10}=[0 2 0 ; 0 2 1];
Block{11}=[1 2 0 ; 1 2 1];
Block{12}=[1 2 2 ; 1 2 3];
Block{13}=[2 2 0 ; 2 2 1];
Block{14}=[2 2 2 ; 2 2 3];
Block{15}=[3 2 0 ; 3 2 1];
Block{16}=[3 2 2 ; 3 2 3];
Block{17}=[0 3 0 ; 0 3 1];
Block{18}=[0 3 3 ];
Block{19}=[1 3 0 ; 1 3 1 ; 1 3 2 ; 1 3 3 ; 2 3 2 ; 2 3 3];
Block{20}=[3 3 0 ; 3 3 1 ; 3 3 2 ; 3 3 3];
Block{21}=[0 2 2 ; 0 2 3 ; 0 3 2 ];  % STOP codon

% fill in code and compute sizes of blocks
Code=zeros(4,4,4);
BlockSize=zeros(21,1);
for r=1:21
    BlockSize(r)=size(Block{r},1);
    for i=1:BlockSize(r)
        Code(Block{r}(i,1)+1,Block{r}(i,2)+1,Block{r}(i,3)+1)=r;
    end
end

if (suppression==1)
 % iterate over blocks
    for t1=1:4
        for t2=1:4
            if (Code(t1,t2,3)==21 && Code(t1,t2,4)==21)
                Code(t1,t2,3)=Code(t1,t2,1);
                Code(t1,t2,4)=Code(t1,t2,1);
            elseif (Code(t1,t2,3)==21)
                Code(t1,t2,3)=Code(t1,t2,4);
            elseif (Code(t1,t2,4)==21)
                error('this should not happen');
            end
        end
    end
end


%% compute B-matrices

% U=0, C=1, A=2, G=3

% transitions (edges between U-C and A-G)
Btransit1=zeros(21,21);
Btransit2=zeros(21,21);
Btransit3=zeros(21,21);
% transversions (edges between U-{AG} and A-{UC})
Btransver1=zeros(21,21);
Btransver2=zeros(21,21);
Btransver3=zeros(21,21);

for i=1:4
    for j=1:4
        for k=1:4            
            r=Code(i,j,k);
            for m=1:3  % three possible mutations
                inew=mod(i+m-1,4)+1;
                jnew=mod(j+m-1,4)+1;               
                knew=mod(k+m-1,4)+1;
                if istransit(i,inew)
                    Btransit1(r,Code(inew,j,k))=Btransit1(r,Code(inew,j,k))+1;
                else
                    Btransver1(r,Code(inew,j,k))=Btransver1(r,Code(inew,j,k))+1;
                end
                
                if istransit(j,jnew)
                    Btransit2(r,Code(i,jnew,k))=Btransit2(r,Code(i,jnew,k))+1;
                else
                    Btransver2(r,Code(i,jnew,k))=Btransver2(r,Code(i,jnew,k))+1;
                end
                
                if istransit(k,knew)
                    Btransit3(r,Code(i,j,knew))=Btransit3(r,Code(i,j,knew))+1;
                else
                    Btransver3(r,Code(i,j,knew))=Btransver3(r,Code(i,j,knew))+1;
                end
            end
        end
    end
end

%% number of neighbors between blocks
B1 = Btransit1 + Btransver1;
B2 = Btransit2 + Btransver2;
B3 = Btransit3 + Btransver3;
B= B1 + B2 + B3;

%% define B variants
% Bst: B with standard all-1 weights
Bst1=Btransit1 + Btransver1;
Bst2=Btransit2 + Btransver2;
Bst3=Btransit3 + Btransver3;
Bst=Bst1+Bst2+Bst3;

% BFH: B with Freeland & Hurst weights
BFH1=    Btransit1 + 0.5*Btransver1;
BFH2=0.5*Btransit2 + 0.1*Btransver2;
BFH3=    Btransit3 +     Btransver3;
BFH=BFH1+BFH2+BFH3;

% trim those matrices to 20 x 20 (get rid of the STOP codon row / column)
Bst = Bst(1:20,1:20);
Bst1 = Bst1(1:20,1:20);
Bst2 = Bst2(1:20,1:20);
Bst3 = Bst3(1:20,1:20);

BFH = BFH(1:20,1:20);
BFH1 = BFH1(1:20,1:20);
BFH2 = BFH2(1:20,1:20);
BFH3 = BFH3(1:20,1:20);


%% cleanup
clear i j k m r inew jnew knew 