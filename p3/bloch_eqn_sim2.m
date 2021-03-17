% Bloch Equation Simulation, Excercise B-1b
% -----------------------------------------
% 
clc
df = 0;		% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.
TE = 1;		% ms. Echo time
TR = 500;	% ms. excitation
flip = pi/3;	% radians. 60 degrees


% ----------- B1-a -------------

M = [0;0;1];
Rflip = yrot(flip);
[Ate,Bte] = freeprecess(TE,T1,T2,df);
M = Rflip*M;	% Magnetization after tip.
M = Ate*M+Bte	% Magnetization at TR


% ----------- B1-b -------------

M = [0;0;1];
[Atr,Btr] = freeprecess(TR,T1,T2,df);
M = Rflip*M;	% Magnetization after tip.
M = Atr*M+Btr;	% Magnetization at TR

%%
% Bloch Equation Simulation, Excercise B-1c
% -----------------------------------------
% 

df = 0;		% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.
TE = 1;		% ms.
dT = 1;
TR = 500;	% ms.
flip = pi/3;	% radians.
Ntr = round(TR/dT);
Nex = 10;	% 20 excitations.

M = [0;0;1];
Rflip = yrot(flip);
[A1,B1] = freeprecess(dT,T1,T2,df);


M(1,Nex*Ntr)=0;	%	Allocate to record all M's.
		% 	Not necessary, but makes program faster.

Mcount=1;
for n=1:Nex
	M(:,Mcount) = Rflip*M(:,Mcount);	

	for k=1:Ntr
		Mcount=Mcount+1;
		M(:,Mcount)=A1*M(:,Mcount-1)+B1;
	end;
end;

time = [0:Mcount-1]*dT;
plot(time,M(1,:),'b-',time,M(2,:),'r--',time,M(3,:),'g-.');
legend('M_x','M_y','M_z');
xlabel('Time (ms)');
ylabel('Magnetization');
axis([min(time) max(time) -1 1]);
grid on;


%%
% Bloch Equation Simulation, Excercise B-1d
% -----------------------------------------
% 

df = 0;		% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.
TE = 1;		% ms.
TR = 500;	% ms.
flip = pi/3;	% radians.


M = [0;0;1];
Rflip = yrot(flip);
[Atr,Btr] = freeprecess(TR,T1,T2,df);

% M1 = Atr * Rflip * M + Btr.
%
% But M1=M in steady state, so
%
% 	M = Atr*Rflip * M + Btr.
%	(I-Atr*Rflip)*M = Btr.

Mss = inv(eye(3)-Atr*Rflip)*Btr



%%
% Bloch Equation Simulation, Excercise B-1e
% -----------------------------------------
% 

df = 0;		% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.
TE = 1;		% ms.
TR = 500;	% ms.
flip = pi/3;	% radians.


M = [0;0;1];
Rflip = yrot(flip);
[Atr,Btr] = freeprecess(TR,T1,T2,df);
[Ate,Bte] = freeprecess(TE,T1,T2,df);
[Atetr,Btetr] = freeprecess(TR-TE,T1,T2,df);

%	Calculation using B-1d first:

Mss = inv(eye(3)-Atr*Rflip)*Btr;
Mte1 = Ate*Rflip*Mss+Bte %ss magnetization just before excitation


% 	Direct calculation at TE

% 	Starting at TE, M=M1
%	At TR, M=M2, and M2=Atetr*M1+Btetr.
%	At TE, M=M3, and M3=Ate*Rflip*M2+Bte.
%			M3=Ate*Rflip*(Atetr*M1+Btetr)+Bte.
%
%	But M3=M1=Mte2 in steady state:

Mte2 = inv(eye(3)-Ate*Rflip*Atetr)* (Ate*Rflip*Btetr+Bte) %ss magnetization at echotime TE





