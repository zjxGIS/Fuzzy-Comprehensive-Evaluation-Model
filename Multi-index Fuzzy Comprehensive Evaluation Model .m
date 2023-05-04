clc;
clear;
%%Multi-index Fuzzy Comprehensive Evaluation Model with In-formation Entropy of Alfalfa Salt Tolerance Based on LiDAR Data and Hyperspectral Image Data
%the evaluation index system, X=[[X11,X12,...,X1m],...,[Xn1,Xn2,...,Xnm]]
%the evaluation set, U=[Very tolerant,Intermidiate,Susceptible,Very susceptible]
%Decide whether the evaluation index is positive or negative。
%Weight-the entropy weight method
%The membership function(trapezoidal membership function) is selected according to the larger, smaller or intermediate type of the factor, 
%The fuzzy comprehensive evaluation matrix R of indicators of each sample was calculated.
%For each sample, the fuzzy operator of multiplication and addition was used for fuzzy transformation, 
%and the weight and fuzzy comprehensive evaluation matrix were synthesized into a fuzzy vector b.
%The maximum membership principle was used to judge the salt tolerance of the samples.

%% Input parameters
%'row' is the sample, 'column' is the factor set
F=xlsread('.\table.xlsx','Sheet_name');
X=F(1:24,8:26);

%% Decide whether the evaluation index is positive(1) or negative(-1)。
global A 
A=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
%the evaluation set, U=[Very tolerant(100),Intermidiate(75),Susceptible(50),Very susceptible(25)]
global U 
U=[100 75 50 25];
yinsu=size(A,2);
pingjia=size(U,2);

%% sortdata
global sortdata 
sortdata=zeros(size(X,1),size(X,2));
for a = 1:size(X,2)
    if A(a)==1
        raw=sort(X,'ascend');
        sortdata(:,a)=raw(:,a);
    elseif A(a)==-1
        raw=sort(X,'descend');
        sortdata(:,a)=raw(:,a);
    end
end
sortdata;


%% weight
w=quanzhong(X);

%% fuzzy vector
normalized_data=normalized(X);
B=zeros(size(X,1),size(U,2));
for n=1:size(X,1)
    rnm=zeros(size(X,2),size(U,2));
    for m=1:size(X,2)
        rnm(m,:)=lishudu(normalized_data(n,m),m);
    end
    B(n,:)=fuzzymm(w,rnm);
end
%%judge the salt tolerance of the samples
fce=max_memb(B);



%% function
%% the entropy weight method
function W = quanzhong(X)
[n,m]=size(X); 
[X,ps]=mapminmax(X',0,1);
ps.min=min(min(X)); 
ps.max=max(max(X)); 
ps.range=ps.max-ps.min; 
X=mapminmax(X,ps);
X=mapminmax('reverse',X,ps); 
X=X';  
p=zeros(n,m);
for i=1:n
    for j=1:m
        p(i,j)=X(i,j)./sum(X(:,j));
    end
end
e=zeros(1,m);
k=1/log(n);
for j=1:m
    e(j)=-k*sum(p(:,j).*log(p(:,j)));
end
d=ones(1,m)-e;  
W=d./sum(d);   
end

%% normalized_data
function [normalized_data] = normalized(source_data)
[source_data,ps]=mapminmax(source_data',0,1);
normalized_data=source_data';
end

%% the membership degree
function R = lishudu(X,a)
global A
global U
R=zeros(1,size(U,2));
if A(a)==1 
    R(1)=trapmf(X,[0.7,0.8,1,1.2]);
    R(2)=trapmf(X,[0.45,0.55,0.7,0.8]);
    R(3)=trapmf(X,[0.2,0.3,0.45,0.55]);
    R(4)=trapmf(X,[0,0,0.2,0.3]);
end
if A(a)==-1 
    R(1)=trapmf(X,[0,0,0.2,0.3]);
    R(2)=trapmf(X,[0.2,0.3,0.45,0.55]);
    R(3)=trapmf(X,[0.45,0.55,0.7,0.8]);
    R(4)=trapmf(X,[0.7,0.8,1,1.2]);
end
end

%% the fuzzy operator of multiplication and addition
function B=fuzzymm(w,r)
[m,n]=size(r);
B=zeros(1,n);
for j=1:n
    tmp=0;
    for i=1:m
        tmp=tmp+(w(i)*r(i,j));
    end
    B(j)=tmp;
end
sum_b=sum(B);
B=B./sum_b;
end

%% The maximum membership principle
function fce=max_memb(b)
global sortdata
global U
fce=zeros(size(sortdata,1),1);
for k=1:size(b,1)
    [~,p]=max(b(k,:));
    fce(k,1)=U(p);
end
end
    
