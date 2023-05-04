clc;
%%#ģ���ۺ�����
%��Ҫ�趨����ָ�ꡢ����ָ����м�ָ�ꡣ
%#���ؼ�A=[X1,X2,...,Xm],A=[3���ṹ������20��ֲ��ָ��]
%#���ۼ�U=[V1,V2,...,Vn]��U=[���ʣ����ã��еȣ��ϲ�]
%#Ȩ��-��ֵ��w
%#�������ص�ƫ���͡�ƫС�ͻ�ƫ�м���ѡ��������������������������������ģ�����۾���R
%#����ģ���ϳ����ӽ���ģ�����㲢��һ��B=W*R��*�Ǻϳ�����,ÿ������������һ��ģ�����������ж���[0,1]���䡣
%#ϵͳ�ܵ÷�

%% �����������Ϊ����������Ϊ���ؼ�
F=xlsread('.\xiu2\1_3.xlsx','19����');
%X=F(1:12,5:11);
%X=F(1:24,4:7);
X=F(1:24,8:26);

%% ָ������������ָ��(1)���Ǹ���ָ��(-1)
global A %ȫ�ֱ���A-���ؼ������������
A=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
%A=[1 1 1 1];
%���ۼ�-4���ȼ�
global U %ȫ�ֱ���U-���ۼ�-��1���ʣ�2���ã�3�еȣ�4�ϲ
U=[100 75 50 25];
yinsu=size(A,2);
pingjia=size(U,2);

%% �������������ϱ��±��ȷ��
%�Ծ������������򣨸���ָ��-1��sort(X,1,'descend')
%�Ծ����н�����������ָ��1��sort(X,1,'ascend')
%size(X,1)���������,��������
%size(X,2)���������,��������
global sortdata %������ָ��͸���ָ���������
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


%% Ȩ��-��ֵ��
w=quanzhong(X);

%% ���㵥���ص�ģ���ϳ�ֵB
%nΪ����
%mΪ����
%normalized_data=zeros(size(X,1),size(X,2));
normalized_data=normalized(X);
B=zeros(size(X,1),size(U,2));
for n=1:size(X,1)
    %rnm�ǵ�һ�����������Ⱦ�����һ�����ؼ����������ۼ�������С�ľ���
    rnm=zeros(size(X,2),size(U,2));
    %һ������һ�������ļ���
    for m=1:size(X,2)
        %X(n,m):һ��������ĳһ����
        %lishudu(X(n,m),m):�����ص���������1�����ۼ�������
        rnm(m,:)=lishudu(normalized_data(n,m),m);
    end
    %ģ�����Ӻϳ�
    B(n,:)=fuzzymm(w,rnm);
end
fce=max_memb(B);



%�Զ��庯��
%% ��ֵ��Ȩ��
function W = quanzhong(X)
[n,m]=size(X); % X����n������, m��ָ��

%% ���ݵĹ�һ������,����ָ��ͬ�ʻ�
%���ڸ���ָ��ļ�����λ����ͳһ������������Ǽ����ۺ�ָ��ǰ���ȶ����ǽ��б�׼��������
%��ָ��ľ���ֵת��Ϊ���ֵ����ȡͬ�򣬴Ӷ�������ͬ��ָ��ֵ��ͬ�ʻ����⣬���ң�����
%����ָ��͸���ָ����ֵ����ĺ��岻ͬ������ָ����ֵԽ��Խ�ã�����ָ����ֵԽ��Խ�ã���
%��ˣ����ڸߵ�ָ�������ò�ͬ���㷨�������ݱ�׼������
%ע��
%Min-Max��׼��������׼��
%mapminmax(X,0,1)ÿһ����һ����һ������
%�Ե���������һ�������ԭ������ת��
[X,ps]=mapminmax(X',0,1);
ps.min=min(min(X)); % ��һ�������Сֵ
ps.max=max(max(X)); % ��һ��������ֵ
ps.range=ps.max-ps.min; % ��һ����ļ���,����������ֵ, ������������
X=mapminmax(X,ps);
X=mapminmax('reverse',X,ps); % ����һ��, �ص�ԭ����
X=X';  % XΪ��һ���������
%% ������۶����ڸ�ָ���µı�ֵ
%%�����j��ָ���£���i����¼ռ��ָ��ı���p(i,j)
p=zeros(n,m);
for i=1:n
    for j=1:m
        %i������-�У�j������-��
        %sum(X(:,j)����һ�е��ܺ�
        p(i,j)=X(i,j)./sum(X(:,j));
    end
end
%% �����j��ָ�����ֵe(j)
e=zeros(1,m);
k=1/log(n);
for j=1:m
    e(j)=-k*sum(p(:,j).*log(p(:,j)));
end

%% ͨ����ֵ�����ָ���Ȩ��
d=ones(1,m)-e;  % ������Ϣ�������
W=d./sum(d);    % ��Ȩֵw
end

%% ���ݱ�׼��
function [normalized_data] = normalized(source_data)
[source_data,ps]=mapminmax(source_data',0,1);
%ps.min=min(min(source_data)); % ��һ�������Сֵ
%ps.max=max(max(source_data)); % ��һ��������ֵ
%ps.range=ps.max-ps.min; % ��һ����ļ���,����������ֵ, ������������
%source_data=mapminmax(source_data,ps);
%source_data=mapminmax('reverse',source_data,ps); % ����һ��, �ص�ԭ����
normalized_data=source_data';  % XΪ��һ���������
end

%% �����ض����ۼ���������
function R = lishudu(X,a)
%XΪһ��������ĳһ����
%aΪ�ڼ�������
global A
global U
%R��������-���ؼ�����������-���ۼ�����
R=zeros(1,size(U,2));
if A(a)==1 %����ָ��
    %trapmf(X,[a,b,c,d]),���У�a��dΪ�����±꣬b��cΪ�����ϱ꣬XΪ��������ֵ�����
    %��������ȡ������round(1.7)=2
    %�������κ������ϱ���±�
    R(1)=trapmf(X,[0.7,0.8,1,1.2]);
    R(2)=trapmf(X,[0.45,0.55,0.7,0.8]);
    R(3)=trapmf(X,[0.2,0.3,0.45,0.55]);
    R(4)=trapmf(X,[0,0,0.2,0.3]);
end
if A(a)==-1 %����ָ��
    R(1)=trapmf(X,[0,0,0.2,0.3]);
    R(2)=trapmf(X,[0.2,0.3,0.45,0.55]);
    R(3)=trapmf(X,[0.45,0.55,0.7,0.8]);
    R(4)=trapmf(X,[0.7,0.8,1,1.2]);
end
end

%% ģ������ϳ�����
function B=fuzzymm(w,r)
%������򣬼�Ȩƽ���ͣ�
%���룬wΪȨ�ؾ���1�����ظ�������r�����������Ⱦ������ظ���m�����ۼ�����n��
%�����BΪ������ģ���ۺ�����������1�����ۼ�������
[m,n]=size(r);
B=zeros(1,n);
for j=1:n
    tmp=0;
    for i=1:m
        tmp=tmp+(w(i)*r(i,j));
    end
    B(j)=tmp;
end
%��һ������
sum_b=sum(B);
B=B./sum_b;
end

%% ���۹���-���������
function fce=max_memb(b)
global sortdata
global U
fce=zeros(size(sortdata,1),1);
for k=1:size(b,1)
    [~,p]=max(b(k,:));
    fce(k,1)=U(p);
end
end
    