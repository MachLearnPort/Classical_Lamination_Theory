% Classical Lamination Theory for determining the Stiffness of the
% composite laminates
% Sept 2, 2016
% Jeff Schwartzentruber

%***** ONLY VALID FOR PLY LENGTHS UNDER 1mm, other wise, you need to change
%the constant SConst in order to proper magnitude, i.e. for mm:1000, for cm:
%100 for m:11

clc;
clear;

SConst=1000;
%% Prompts
prompt={'1. How many plies in the composite being analyzed:','2. Is the laminate symmetric? [1= Yes, 2= No]: ','3. Are all plies the same thickness? [1= Yes, 2= No]:','4. Are all plies of the laminate made of the same material? [1= Yes, 2= No]:','5. Would you like to calculate the macomechanical properties based on direct values (1) or fiber/matrix constituants (2)?'};
def={'6','1','1','1','1'};
TITLE='Preliminary Questions';
line=1;
ANSWER=inputdlg(prompt,TITLE,line,def);
convertc=char(ANSWER);
prec=str2num(convertc);
n=prec(1);
Quest2=prec(2);
Quest3=prec(3);
Quest4=prec(4);
Quest5=prec(5);


%% Thickness Dialog(s)
if Quest3==2
    for i=1:1:n;
        sQ4=['Height of lamina ',num2str(i),' in <m>:'];
        prompt={sQ4};
        TITLE='Lamina Heights';
        line=1;
        ANSWER=inputdlg(prompt,TITLE,line);
        convertc=char(ANSWER);
        prec=str2num(convertc);
        h(i)=prec;
    end
else
    sQ4=['Height of lamina used in <m>:'];
    prompt={sQ4};
    def={'0.0002'};
    TITLE='Lamina Height';
    line=1;
    ANSWER=inputdlg(prompt,TITLE,line,def);
    convertc=char(ANSWER);
    prec=str2num(convertc);
    h=prec;
    h(1:n)=h;
end

%% Mat props Dialog(s) based on consituants
if Quest4==2 && Quest5==1
    for i=1:1:n;
        %Long Modulus
        sQ4=['Longitudinal Modulus of lamina ',num2str(i),' in <GPa>:'];
        prompt={sQ4};
        TITLE='Lamina Material Properties';
        line=1;
        ANSWER=inputdlg(prompt,TITLE,line);
        convertc=char(ANSWER);
        prec=str2num(convertc);
        EL(i)=prec*10^9; %in Pa
        %Transvers Modulus
        sQ4=['Transvers Modulus of lamina ',num2str(i),' in <GPa>:'];
        prompt={sQ4};
        TITLE='Lamina Material Properties';
        line=1;
        ANSWER=inputdlg(prompt,TITLE,line);
        convertc=char(ANSWER);
        prec=str2num(convertc);
        ET(i)=prec*10^9; % in Pa
        %Long Possions Ratio
        sQ4=['Logitudinal Possions Ratio ',num2str(i),' :'];
        prompt={sQ4};
        TITLE='Lamina Material Properties';
        line=1;
        ANSWER=inputdlg(prompt,TITLE,line);
        convertc=char(ANSWER);
        prec=str2num(convertc);
        vL(i)=prec;
        %Transever Possions Ratio
        sQ4=['Transverse Possions Ratio ',num2str(i),' :'];
        prompt={sQ4};
        TITLE='Lamina Material Properties';
        line=1;
        ANSWER=inputdlg(prompt,TITLE,line);
        convertc=char(ANSWER);
        prec=str2num(convertc);
        vT(i)=prec;
        %Shear Modulus
        sQ4=['Shear Modulus ',num2str(i),' in <GPa>:'];
        prompt={sQ4};
        TITLE='Lamina Material Properties';
        line=1;
        ANSWER=inputdlg(prompt,TITLE,line);
        convertc=char(ANSWER);
        prec=str2num(convertc);
        GLT(i)=prec*10^9;
        %Orientation
        sQ4=['Lamina orientation ',num2str(i),' in <deg>:'];
        prompt={sQ4};
        TITLE='Lamina Material Properties';
        line=1;
        ANSWER=inputdlg(prompt,TITLE,line);
        convertc=char(ANSWER);
        prec=str2num(convertc);
        theta(i)=prec;
    end
end

if Quest4==1 && Quest5==1
    prompt={'Longitudinal Modulus in <GPa>','Transverse Modulus in <GPa>','Longitudinal Possions ratio','Shear Modulus in <GPa>'};
    %def={'122','8','0.31','3.74'};
    def={'122','8','0.31','3.74'};
    TITLE='Preliminary Questions';
    line=1;
    ANSWER=inputdlg(prompt,TITLE,line,def);
    convertc=char(ANSWER);
    prec=str2num(convertc);
    EL=prec(1)*10^9;
    ET=prec(2)*10^9;
    vL=prec(3);
    GLT=prec(4)*10^9;
    for i=1:1:n
        %Orientation
        sQ4=['Lamina ',num2str(i),' orientation  in <deg>:'];
        prompt={sQ4};
        TITLE='Lamina Material Properties';
        line=1;
        ANSWER=inputdlg(prompt,TITLE,line);
        convertc=char(ANSWER);
        prec=str2num(convertc);
        theta(i)=prec(1);
    end
    %     theta=[0 0];
    EL(1:n)=EL;
    ET(1:n)=ET;
    vL(1:n)=vL;
    vT(1:n)=(ET(1:n).*vL(1:n))./EL(1:n);
    vT(1:n)=vT; %Possions ratio in the Traverse direction (90 direction) - equ 5.79 Agarwal
    GLT(1:n)=GLT;
end

if Quest4==2 && Quest5==2
    for i=1:1:n
        prompt={'Ef (Modulus of Fibre)(GPa)','Em Modulus of Matrix (GPa)','Gf Shear Modulus of Fibre (GPa)','Gm Shear Modulus of Matrix (GPa)','vf Poisson of fibre','vm Poisson of Matrix','Vf Volume fraction','Df Density of fibre (kg/m^3)','Dm Density of Matrix (kg/m^3)'};
        TITLE='Consitant Properties of lamina' ,num2str(i),'in the laminate stack';
        line=1;
        ANSWER=inputdlg(prompt,TITLE,line);
        convertc=char(ANSWER);
        prec=str2num(convertc);
        double(prec);
        Ef(i)=prec(1)*10^9; %Elastic Modulus of Fiber
        Em(i)=prec(2)*10^9; %Elastic Modulus of Matrix
        Gf(i)=prec(3)*10^9; %Shear Modulus of Fiber
        Gm(i)=prec(4)*10^9; %Elastic Modulus of Matrix
        vf(i)=prec(5); %Poison's Ratio of Fiber
        vm(i)=prec(6); %Poison's Ratio of Matrix
        V(i)=prec(7);  %Volume Fiber Fraction
        Df(i)=prec(8); %Density of Fiber
        Dm(i)=prec(9);  %Density Ratio of Matrix
        
        % Calc material proeprties
        EL(i)=Ef*V+Em*(1-V);
        ET(i)=Em*(Ef+Em+(Ef-Em)*V)/(Ef+Em-(Ef-Em)*V);
        vL(i)=vf*V+vm*(1-V);
        vT(i)=vf*V+vm*(1-V)*(1+vm-vL*Em/EL)/(1-vm^2+vm*vL*Em/EL);
        GLT(i)=Gm*(Gf+Gm+(Gf-Gm)*V)/(Gf+Gm-(Gf-Gm)*V);
        GTL(i)=ET/(2*(1+vT));
        Den(i)=Df*V+Dm*(1-V);
    end
end

if Quest4==1 && Quest5==2
    prompt={'Ef (Modulus of Fibre)(GPa)','Em Modulus of Matrix (GPa)','Gf Shear Modulus of Fibre (GPa)','Gm Shear Modulus of Matrix (GPa)','vf Poisson of fibre','vm Poisson of Matrix','Vf Volume fraction (NOT PERCENTAGE)','Df Density of fibre (kg/m^3)','Dm Density of Matrix (kg/m^3)'};
    TITLE='Consitant Properties of lamina' ;num2str(i),'in the laminate stack';
    line=1;
    ANSWER=inputdlg(prompt,TITLE,line);
    convertc=char(ANSWER);
    prec=str2num(convertc);
    double(prec);
    Ef=prec(1)*10^9; %Elastic Modulus of Fiber
    Em=prec(2)*10^9; %Elastic Modulus of Matrix
    Gf=prec(3)*10^9; %Shear Modulus of Fiber
    Gm=prec(4)*10^9; %Elastic Modulus of Matrix
    vf=prec(5); %Poison's Ratio of Fiber
    vm=prec(6); %Poison's Ratio of Matrix
    V=prec(7);  %Volume Fiber Fraction
    Df=prec(8); %Density of Fiber
    Dm=prec(9);  %Density Ratio of Matrix
    
    % Calc material proeprties
    EL=Ef*V+Em*(1-V);
    ET=Em*(Ef+Em+(Ef-Em)*V)/(Ef+Em-(Ef-Em)*V);
    vL=vf*V+vm*(1-V);
    vT=vf*V+vm*(1-V)*(1+vm-vL*Em/EL)/(1-vm^2+vm*vL*Em/EL);
    GLT=Gm*(Gf+Gm+(Gf-Gm)*V)/(Gf+Gm-(Gf-Gm)*V);
    GTL=ET/(2*(1+vT));
    Den=Df*V+Dm*(1-V);
    for i=1:1:n
        %Orientation
        sQ4=['Lamina ',num2str(i),' orientation  in <deg>:'];
        prompt={sQ4};
        TITLE='Lamina Material Properties';
        line=1;
        ANSWER=inputdlg(prompt,TITLE,line);
        convertc=char(ANSWER);
        prec=str2num(convertc);
        theta(i)=prec(1);
    end
    EL(1:n)=EL;
    ET(1:n)=ET;
    vL(1:n)=vL;
    vT(1:n)=vT;
    GLT(1:n)=GLT;
    GTL(1:n)=GTL;
    Den(1:n)=Den;
end

%% CLT

%Stiffness Matrix Constants for ply - equation 5.78 Argwal
for i=1:1:n
    Q11(i)=EL(i)/(1-vL(i)*vT(i));
    Q22(i)=ET(i)/(1-vL(i)*vT(i));
    Q12(i)=((vL(i)*ET(i))/(1-vL(i)*vT(i))); %Or Q12=((vT*EL)/(1-vL*vT)) - see which one gives a better answer, the should be the same
    Q66(i)=GLT(i);
    
    qbar(1,1,i)=Q11(i);
    qbar(1,2,i)=Q12(i);
    qbar(1,3,i)=0;
    qbar(2,1,i)=Q12(i);
    qbar(2,2,i)=Q22(i);
    qbar(2,3,i)=0;
    qbar(3,1,i)=0;
    qbar(3,2,i)=0;
    qbar(3,3,i)=Q66(i);
end

%Stiffness Matrix for Arbitrary Axis - Equation 5.95 Agarwal and Generation
%of Q matrix for laminate stack
for i=1:1:n
    Q11b(i)=(Q11(i)*(cosd(theta(i)))^4)+(Q22(i)*(sind(theta(i)))^4)+(2*(Q12(i)+2*Q66(i)))*((sind(theta(i)))^2)*((cosd(theta(i)))^2);
    Q22b(i)=(Q11(i)*(sind(theta(i)))^4)+(Q22(i)*(cosd(theta(i)))^4)+(2*(Q12(i)+2*Q66(i)))*((sind(theta(i)))^2)*((cosd(theta(i)))^2);
    Q12b(i)=(Q11(i)+Q22(i)-4*Q66(i))*((sind(theta(i)))^2)*((cosd(theta(i)))^2)+Q12(i)*(((cosd(theta(i)))^4)+((sind(theta(i)))^4));
    Q66b(i)=(Q11(i)+Q22(i)-2*Q12(i)-2*Q66(i))*((sind(theta(i)))^2)*((cosd(theta(i)))^2)+Q66(i)*(((sind(theta(i)))^4)+((cosd(theta(i)))^4));
    Q16b(i)=(Q11(i)-Q12(i)-2*Q66(i))*(((cosd(theta(i)))^3))*(sind(theta(i)))-(Q22(i)-Q12(i)-2*Q66(i))*(cosd(theta(i)))*((sind(theta(i)))^3);
    Q26b(i)=(Q11(i)-Q12(i)-2*Q66(i))*(cosd(theta(i)))*((sind(theta(i)))^3)-(Q22(i)-Q12(i)-2*Q66(i))*((cosd(theta(i)))^3)*(sind(theta(i)));
    
    %Construct Q bar matrix
    Qbar(1,1,i)=Q11b(i);
    Qbar(1,2,i)=Q12b(i);
    Qbar(1,3,i)=Q16b(i);
    Qbar(2,1,i)=Q12b(i);
    Qbar(2,2,i)=Q22b(i);
    Qbar(2,3,i)=Q26b(i);
    Qbar(3,1,i)=Q16b(i);
    Qbar(3,2,i)=Q26b(i);
    Qbar(3,3,i)=Q66b(i);
end

%Determination of heights for A B D matrix generation
MP=sum(h)/2; %midplan location
%depth of plys
H(1)=h(1); %Inatlize H and give its first height the height on the top ply (1st ply(
for i=2:1:n
    H(i)=H(i-1)+h(i);
end
Hdiff(1:n)=abs(H(1:n)-MP);
[Min,Id] = min(Hdiff);

% % number of segments for heights (hk's)
% if Min==0
%     Seg=n;
% else
%     Seg=n+1;
% end

% %upward h determination hk(1) is hk(0) when comparing Agarwal
% if Id<2
%     hku(1)=-MP;
% end
% for i=2:Id
%     hku(i)=-(MP-H(i-1));
%     hku(1)=-MP;
%
% end
% %downward
% for i=Id:n
%     hkd(i)=H(i)-MP;
% end
% hkd(hkd==0)=[];
% hk=[hku,hkd];

% hk(1)=-MP;
% for i=1:n
%     hk(i+1)=H(i)-MP;
% end
t=h(1);
hk(n+1)=((n*t)/2)*SConst;
for i=1:n
    hk(i)=(i-n/2-1)*t*SConst; %*10^3 becuase this is suppose to be more of a ratio, so that we can get proper values at the end
end

A=zeros(3);
B=zeros(3);
D=zeros(3);

for k =1:n
    for i=1:3
        for j=1:3
            A(i,j) = round(Qbar(i,j,k) * (hk(k+1) - hk(k)) + A(i,j));
            Arec(i,j,k)=A(i,j);
            B(i,j) = round((Qbar(i,j,k) * (hk(k+1)^2 - hk(k)^2)) + B(i,j));
            Brec(i,j,k)=B(i,j);
            D(i,j) = round((Qbar(i,j,k) * (hk(k+1)^3 - hk(k)^3)) + D(i,j));
            Drec(i,j,k)=D(i,j);
        end
    end
end
%The fractions outfront of the B and D martix on equation 6.20 in Agarwal
B=B/2;
D=D/3;

disp(A);
E11=A(1,1);
E22=A(2,2);

%% Calculation of Elastic Modulus, possions ratio and shear modulus for Blister code

