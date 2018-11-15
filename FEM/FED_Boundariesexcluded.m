clear all
format long

tab=zeros(7,8);

call=zeros(7,1);



%%Defining Variables%%
for bm=1:4
    
clearvars -except bm tab call
OC = load(strcat("ele",num2str(bm),".txt"));
OA = load(strcat("node",num2str(bm),".txt"));
totcount=0;
numtri=OC(1,1);
numnodes=OA(1,1);
trimat=zeros(3,3);
finA=zeros(3,3,numtri);
intA=zeros(3,3,numtri);
area=zeros(1,numtri);

triunalt=zeros(3,3,numtri);
dummy=0;
dummy1=0;

%The following for loop will count the 
%number of unique interior nodes
for i=2:numnodes+1
    if OA(i,4)==0
        dummy=dummy+1;
        storedintnodes(dummy,:)=OA(i,2:3);
        numintnodes=dummy;
    end
end
%Exterior nodes
for i=2:numnodes+1
    if OA(i,4)==1
        dummy1=dummy1+1;
        storedexnodes(dummy1,:)=OA(i,2:3);
        
    end
end
N=numintnodes;

%%Triangle Calculation%%
for i=1:numtri
    
    for p=1:3
        trimat(p,:)=OA(OC(i+1,p+1)+2,2:end);
        
    end
        a=sqrt((trimat(1,1)-trimat(2,1))^2+(trimat(1,2)-trimat(2,2))^2);
        b=sqrt((trimat(2,1)-trimat(3,1))^2+(trimat(2,2)-trimat(3,2))^2);
        c=sqrt((trimat(3,1)-trimat(1,1))^2+(trimat(3,2)-trimat(1,2))^2);
        s=(a+b+c)/2;
        area(i)=sqrt(s*(s-a)*(s-b)*(s-c));
        
    triunalt(:,:,i)=trimat;
    triunaltex(:,:,i)=trimat;
    for o=1:3
        
        if trimat(o,3)==1
            triunalt(o,3,i)=-5;
        
            for k=1:(numnodes-numintnodes)
                if triunalt(o,1:2,i)==storedexnodes(k,1:2)
                    triunaltex(o,3,i)=k;
                end
            end
        end
        
        if trimat(o,3)~=1
            for k=1:numintnodes
                if triunalt(o,1:2,i)==storedintnodes(k,1:2)
                    triunalt(o,3,i)=k;
                end
            end
        end
     end
    


    for k=1:3
        trimat(:,3)=1;
        B=transpose(zeros(1,3));
        B(k,1)=1;
        finA(k,:,i)=transpose(trimat\B);
    end
end


%%Phi Matrix Calculation%%
for i=1:numtri
    for p=1:3
        for k=1:3
            a1=finA(p,1,i);
            b1=finA(p,2,i);
            a2=finA(k,1,i);
            b2=finA(k,2,i);
            
            intA(p,k,i)=(a1*a2+b1*b2)*area(i);
        end
        
    end
end

t=0;
%Matching Interior Nodes to IntA 
%and Sparse Matrix Calc.
    for i=1:numtri
        la=0;
        cas=0;
        for j=1:3
            if triunalt(j,3,i)~=-5
                cas=cas+1;
                la(1,cas)=j;  
            end
        end
        
        if cas~=0
            for j=1:length(la)
                for k=1:length(la)
                    t=t+1;
                    globali(1,t)=triunalt(la(1,j),3,i);
                    globalj(1,t)=triunalt(la(1,k),3,i);
                    globals(1,t)=intA(la(1,j),la(1,k),i);
                end
            end
        
        end
    end

    

%Creating The A Sparse Matrix
Globala=sparse(globali,globalj,globals,N,N);



%Creating The b Sparse Vector
t=0;
preu_h=zeros(numintnodes,3);
for i=1:numtri
    la=0;
    cas=0;
        for j=1:3
            if triunalt(j,3,i)~=-5
                cas=cas+1;
                la(1,cas)=j;
            end
        end
        if cas~=0
            for j=1:length(la)
                
                    t=t+1;
                    bi(1,t)=triunalt(la(1,j),3,i);
                    bj(1,t)=1;
                    A=[triunalt(1,1,i),triunalt(1,2,i)];
                    B=[triunalt(2,1,i),triunalt(2,2,i)];
                    C=[triunalt(3,1,i),triunalt(3,2,i)];
                    h=10^-14;
                    bo=la(1,j);
                    MAT=finA(:,:,i);
                    
%________________________________TRIANGLE FUNCTION_____________________________________________________                    
n=7;
W=[27/60,3/60,3/60,3/60,8/60,8/60,8/60];
sum=0;
Xt=[0,h,-h/2,-h/2,-h/2,h/4,h/4];
Yt=[0,0,(h/2)*sqrt(3),(-h/2)*(sqrt(3)),0,(h/4)*sqrt(3),(-h/4)*sqrt(3)];
X=[0,0,0,0];
Y=[0,0,0,0];
fnc=[0,0,0,0];

b=(A(1)-B(1))/(sqrt(3)*h);
d=(A(2)-B(2))/(sqrt(3)*h);

a=(A(1)-C(1)-sqrt(3)/2*h*b)/((-3/2)*h);
c=(A(2)-C(2)-sqrt(3)/2*h*d)/((-3/2)*h);

f=C(2)-c*h;
e=C(1)-a*h;

ja=abs(a*d-b*c);
sum1=0;
for it=1:n
    X(it)=a*Xt(it)+b*Yt(it)+e;
    Y(it)=c*Xt(it)+d*Yt(it)+f;
    
    fnc(it)=func(X(it),Y(it))*(((MAT(bo,1)*X(it))+(MAT(bo,2)*Y(it))+MAT(bo,3)))*ja*W(it);
    
    sum=sum+(fnc(it));
    
end

bs(1,t)=sum*(3/4*h^2*sqrt(3));
%______________________________END OF TRIANGLE FUNCTION________________________


            end
        end
        
 
end

Globalb=sparse(bi,bj,bs,N,1);
Globalx=Globala\((Globalb));

%______________________________ Boundary Calculations________________


for i=1:numnodes
    
        if OA(i+1,4)==1
            C_g(i,1)=gunc(OA(i+1,2),OA(i+1,3));
        end
        
            
end


%Graphical Representation
for i=1:(numtri)
    
graphx=triunalt(:,1,i);
graphy=triunalt(:,2,i);
graphz=zeros(1,3);
for k=1:3
    if triunalt(k,3,i)==-5
        graphz(1,k)=0;
    else
        graphz(1,k)=Globalx(triunalt(k,3,i),1);
    end
end

patch(graphx,graphy,graphz,[i*(1/numtri),0.75,0]);

%hold on
end
view(3)
%______________________________ Calculating u-u_h Error __________________________
u_h=0;
sum2=0;
finalerrorsum=0;
aftergq=0;
ca=0;
for i=1:numtri

            A= triunalt(1,1:2,i);
            B= triunalt(2,1:2,i);
            C= triunalt(3,1:2,i);

            W=[27/60,3/60,3/60,3/60,8/60,8/60,8/60];
            sum=0;
            Xt=[0,h,-h/2,-h/2,-h/2,h/4,h/4];
            Yt=[0,0,(h/2)*sqrt(3),(-h/2)*(sqrt(3)),0,(h/4)*sqrt(3),(-h/4)*sqrt(3)];
            X=[0,0,0,0];
            Y=[0,0,0,0];
            fnc=[0,0,0,0];

            b=(A(1)-B(1))/(sqrt(3)*h);
            d=(A(2)-B(2))/(sqrt(3)*h);

            a=(A(1)-C(1)-sqrt(3)/2*h*b)/((-3/2)*h);
            c=(A(2)-C(2)-sqrt(3)/2*h*d)/((-3/2)*h);

            f=C(2)-c*h;
            e=C(1)-a*h;

            ja=abs(a*d-b*c);

            for j=1:7
                u_h=0;
                X(j)=a*Xt(j)+b*Yt(j)+e;
                Y(j)=c*Xt(j)+d*Yt(j)+f;
    
  
                for p=1:3
                    if triunalt(p,3,i)~=-5
                        
                        u_h=u_h+(Globalx(triunalt(p,3,i))*(finA(p,1,i)*X(j)+finA(p,2,i)*Y(j)+finA(p,3,i)));
                        
                    end
                end

                
                sum=sum+((u(X(j),Y(j))-u_h)^2*ja*W(j));
            end
            if u_h==0
                sum=0;
            end
            aftergq=aftergq+(sum*(3/4*h^2*sqrt(3)));    
        
    
 end

%______________________________ end of u-u_h error calculation __________________________


%Tab Creation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tab(bm,1)=(sqrt(aftergq));




if bm>1
    tab(bm,2)=tab(bm,1)/tab(bm-1,1);
 

end
call(bm,1)=bm;
end

