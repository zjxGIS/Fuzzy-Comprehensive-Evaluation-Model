clc;
%%#模糊综合评价
%需要设定正向指标、负向指标和中间指标。
%#因素集A=[X1,X2,...,Xm],A=[3个结构参数，20个植被指数]
%#评价集U=[V1,V2,...,Vn]，U=[优质，良好，中等，较差]
%#权重-熵值法w
%#根据因素的偏大型、偏小型或偏中间型选择隶属函数（梯形隶属函数），计算模糊评价矩阵R
%#利用模糊合成算子进行模糊运算并归一化B=W*R，*是合成算子,每个样本都生成一个模糊向量，所有都在[0,1]区间。
%#系统总得分

%% 输出参数，行为样本数，列为因素集
F=xlsread('.\xiu2\1_3.xlsx','19因素');
%X=F(1:12,5:11);
%X=F(1:24,4:7);
X=F(1:24,8:26);

%% 指明因素是正向指标(1)还是负向指标(-1)
global A %全局变量A-因素集的正负向矩阵
A=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
%A=[1 1 1 1];
%评价集-4个等级
global U %全局变量U-评价集-（1优质，2良好，3中等，4较差）
U=[100 75 50 25];
yinsu=size(A,2);
pingjia=size(U,2);

%% 梯形隶属函数上标下标的确定
%对矩阵按列升序排序（负向指标-1）sort(X,1,'descend')
%对矩阵按列降序排序（正向指标1）sort(X,1,'ascend')
%size(X,1)矩阵的行数,样本数；
%size(X,2)矩阵的列数,因素数；
global sortdata %对正向指标和负向指标进行排序
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


%% 权重-熵值法
w=quanzhong(X);

%% 计算单因素的模糊合成值B
%n为样本
%m为因素
%normalized_data=zeros(size(X,1),size(X,2));
normalized_data=normalized(X);
B=zeros(size(X,1),size(U,2));
for n=1:size(X,1)
    %rnm是单一样本的隶属度矩阵，是一个因素集个数×评价集个数大小的矩阵
    rnm=zeros(size(X,2),size(U,2));
    %一个样本一个样本的计算
    for m=1:size(X,2)
        %X(n,m):一个样本的某一因素
        %lishudu(X(n,m),m):单因素的隶属矩阵（1×评价集个数）
        rnm(m,:)=lishudu(normalized_data(n,m),m);
    end
    %模糊算子合成
    B(n,:)=fuzzymm(w,rnm);
end
fce=max_memb(B);



%自定义函数
%% 熵值法权重
function W = quanzhong(X)
[n,m]=size(X); % X中有n个样本, m个指标

%% 数据的归一化处理,异质指标同质化
%由于各项指标的计量单位并不统一，因此在用它们计算综合指标前，先对它们进行标准化处理，即
%把指标的绝对值转化为相对值，并取同向，从而解决各项不同质指标值的同质化问题，而且，由于
%正向指标和负向指标数值代表的含义不同（正向指标数值越高越好，负向指标数值越低越好），
%因此，对于高低指标我们用不同的算法进行数据标准化处理。
%注：
%Min-Max标准化，离差标准化
%mapminmax(X,0,1)每一行做一个归一化处理
%对单因素做归一化，需对原矩阵做转置
[X,ps]=mapminmax(X',0,1);
ps.min=min(min(X)); % 归一化后的最小值
ps.max=max(max(X)); % 归一化后的最大值
ps.range=ps.max-ps.min; % 归一化后的极差,若不调整该值, 则逆运算会出错
X=mapminmax(X,ps);
X=mapminmax('reverse',X,ps); % 反归一化, 回到原数据
X=X';  % X为归一化后的数据
%% 求各评价对象在各指标下的比值
%%计算第j个指标下，第i个记录占该指标的比重p(i,j)
p=zeros(n,m);
for i=1:n
    for j=1:m
        %i是样本-行，j是因素-列
        %sum(X(:,j)代表一列的总和
        p(i,j)=X(i,j)./sum(X(:,j));
    end
end
%% 计算第j个指标的熵值e(j)
e=zeros(1,m);
k=1/log(n);
for j=1:m
    e(j)=-k*sum(p(:,j).*log(p(:,j)));
end

%% 通过熵值计算各指标的权重
d=ones(1,m)-e;  % 计算信息熵冗余度
W=d./sum(d);    % 求权值w
end

%% 数据标准化
function [normalized_data] = normalized(source_data)
[source_data,ps]=mapminmax(source_data',0,1);
%ps.min=min(min(source_data)); % 归一化后的最小值
%ps.max=max(max(source_data)); % 归一化后的最大值
%ps.range=ps.max-ps.min; % 归一化后的极差,若不调整该值, 则逆运算会出错
%source_data=mapminmax(source_data,ps);
%source_data=mapminmax('reverse',source_data,ps); % 反归一化, 回到原数据
normalized_data=source_data';  % X为归一化后的数据
end

%% 单因素对评价集的隶属度
function R = lishudu(X,a)
%X为一个样本中某一因素
%a为第几个因素
global A
global U
%R矩阵，行数-因素集个数，列数-评价集个数
R=zeros(1,size(U,2));
if A(a)==1 %正向指标
    %trapmf(X,[a,b,c,d]),其中，a，d为梯形下标，b，c为梯形上标，X为待计算数值或矩阵
    %四舍五入取整函数round(1.7)=2
    %计算梯形函数的上标和下标
    R(1)=trapmf(X,[0.7,0.8,1,1.2]);
    R(2)=trapmf(X,[0.45,0.55,0.7,0.8]);
    R(3)=trapmf(X,[0.2,0.3,0.45,0.55]);
    R(4)=trapmf(X,[0,0,0.2,0.3]);
end
if A(a)==-1 %负向指标
    R(1)=trapmf(X,[0,0,0.2,0.3]);
    R(2)=trapmf(X,[0.2,0.3,0.45,0.55]);
    R(3)=trapmf(X,[0.45,0.55,0.7,0.8]);
    R(4)=trapmf(X,[0.7,0.8,1,1.2]);
end
end

%% 模糊矩阵合成运算
function B=fuzzymm(w,r)
%运算规则，加权平均型，
%输入，w为权重矩阵（1×因素个数），r单因素隶属度矩阵（因素个数m×评价集个数n）
%输出，B为单因素模糊综合评价向量（1×评价集个数）
[m,n]=size(r);
B=zeros(1,n);
for j=1:n
    tmp=0;
    for i=1:m
        tmp=tmp+(w(i)*r(i,j));
    end
    B(j)=tmp;
end
%归一化处理
sum_b=sum(B);
B=B./sum_b;
end

%% 评价规则-最大隶属度
function fce=max_memb(b)
global sortdata
global U
fce=zeros(size(sortdata,1),1);
for k=1:size(b,1)
    [~,p]=max(b(k,:));
    fce(k,1)=U(p);
end
end
    