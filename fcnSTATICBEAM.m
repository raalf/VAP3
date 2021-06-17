function [SURF, OUTP, COND, INPU] = fcnSTATICBEAM(SURF, OUTP, COND, INPU, VEHI, FLAG, valTIMESTEP)

valDY = 0.5*INPU.vecSPAN/(INPU.valNSELE-1);

temp_y = (0:valDY:0.5*INPU.vecSPAN)';

[SURF, OUTP, INPU] = fcnBEAMFORCE(SURF, OUTP, COND, INPU, FLAG, valTIMESTEP);

% n = 15;
% L = 16;
% valDY = L/(n-1);
OUTP.vecBEAMFORCE = 5*ones(INPU.valNSELE,1);
OUTP.vecBEAMMOM = zeros(INPU.valNSELE,1);

% INPU.matEIx = [60000*ones(n,1),zeros(n,2)];
% INPU.matGJt = [30000*ones(n,1),zeros(n,1)];

%% Deflection
def_mat = zeros(size(INPU.matEIx,1)-1);

def_mat(1,1) = -2*INPU.matEIx(2,3)*valDY*valDY + 6*INPU.matEIx(2,2)*valDY + 7*INPU.matEIx(2,1);
def_mat(1,2) = INPU.matEIx(2,3)*valDY*valDY - 6*INPU.matEIx(2,2)*valDY - 4*INPU.matEIx(2,1);
def_mat(1,3) = 2*INPU.matEIx(2,2)*valDY + INPU.matEIx(2,1);
def_mat(2,1) = INPU.matEIx(3,3)*valDY*valDY - 2*INPU.matEIx(3,2)*valDY - 4*INPU.matEIx(3,1);
def_mat(2,2) = -2*INPU.matEIx(3,3)*valDY*valDY + 6*INPU.matEIx(3,2)*valDY + 6*INPU.matEIx(3,1);
def_mat(2,3) = INPU.matEIx(3,3)*valDY*valDY - 6*INPU.matEIx(3,2)*valDY - 4*INPU.matEIx(3,1);
def_mat(2,4) = 2*INPU.matEIx(3,2)*valDY + INPU.matEIx(3,1);

for i = 3:size(INPU.matEIx,1)-3
     
    def_mat(i,i-2) = INPU.matEIx(i+1,1);  
    def_mat(i,i-1) = INPU.matEIx(i+1,3)*valDY*valDY - 2*INPU.matEIx(i+1,2)*valDY - 4*INPU.matEIx(i+1,1);
    def_mat(i,i) = -2*INPU.matEIx(i+1,3)*valDY*valDY + 6*INPU.matEIx(i+1,2)*valDY + 6*INPU.matEIx(i+1,1);
    def_mat(i,i+1) = INPU.matEIx(i+1,3)*valDY*valDY - 6*INPU.matEIx(i+1,2)*valDY - 4*INPU.matEIx(i+1,1);
    def_mat(i,i+2) = 2*INPU.matEIx(i+1,2)*valDY + INPU.matEIx(i+1,1);
    
end

def_mat(end-1,end-3) = INPU.matEIx(end-1,1);
def_mat(end-1,end-2) = INPU.matEIx(end-1,3)*valDY*valDY - 2*INPU.matEIx(end-1,2)*valDY - 4*INPU.matEIx(end-1,1);
def_mat(end-1,end-1) = -2*INPU.matEIx(end-1,3)*valDY*valDY + 4*INPU.matEIx(end-1,2)*valDY + 5*INPU.matEIx(end-1,1);
def_mat(end-1,end) = INPU.matEIx(end-1,3)*valDY*valDY - 2*INPU.matEIx(end-1,2)*valDY - 2*INPU.matEIx(end-1,1);

def_mat(end,end-2) = INPU.matEIx(end,1);
def_mat(end,end-1) = -2*INPU.matEIx(end,1);
def_mat(end,end) = INPU.matEIx(end,1);
def_mat = (def_mat)./(valDY^4);

% force_vec = 10*ones(n,1);

OUTP.matDEFGLOB(valTIMESTEP,:) = [0; def_mat\OUTP.vecBEAMFORCE(2:end)]';

%% Torsion
tor_mat = zeros(size(INPU.matGJt,1)-1);

tor_mat(1,1) = -2*INPU.matGJt(2,1)/valDY;
tor_mat(1,2) = INPU.matGJt(2,1) - INPU.matGJt(2,2)*valDY/2;

for i = 2:size(INPU.matGJt,1)-2
    
    tor_mat(i,i-1) = INPU.matGJt(i+1,1) + INPU.matGJt(i+1,2)*valDY/2;
    tor_mat(i,i) = -2*INPU.matGJt(i+1,1);
    tor_mat(i,i+1) = INPU.matGJt(i+1,1) - INPU.matGJt(i+1,2)*valDY/2;
    
end

tor_mat(end,end-1) = 2*INPU.matGJt(end,1);
tor_mat(end,end) = -2*INPU.matGJt(end,1);

% mom_vec = 5*ones(n,1);
tor_mat = (tor_mat)./(valDY*valDY);

OUTP.matTWISTGLOB(valTIMESTEP,:) = [0; -tor_mat\OUTP.vecBEAMMOM(2:end)]';

OUTP.matDEF = [];
OUTP.matTWIST = [];

OUTP.matDEF(1:2,3:size(OUTP.matDEFGLOB,2)+2) = repmat(OUTP.matDEFGLOB(valTIMESTEP,:),2,1);
OUTP.matTWIST(1:2,3:size(OUTP.matTWISTGLOB,2)+2) = repmat(OUTP.matTWISTGLOB(valTIMESTEP,:),2,1);

OUTP.matDEF(1:2,end+1) = repmat(2*OUTP.matDEFGLOB(valTIMESTEP,end) - OUTP.matDEFGLOB(valTIMESTEP,end-1),2,1);
OUTP.matDEF(1:2,end+1) = repmat(3*OUTP.matDEFGLOB(valTIMESTEP,end) - 2*OUTP.matDEFGLOB(valTIMESTEP,end-1),2,1);
OUTP.matDEF(1:2,2) = repmat(OUTP.matDEFGLOB(valTIMESTEP,2),2,1);

OUTP.matTWIST(1:2,end+1) = repmat(OUTP.matTWISTGLOB(valTIMESTEP,end-1),2,1);

% y = [0; cumsum(valDY*ones(size(u,1)-1,1))];
% 
% figure(10)
% plot(y,u,'-k','LineWidth',1.5)
% grid on
% grid minor
% 
% figure(11)
% plot(y,rad2deg(theta),'-k','LineWidth',1.5)
% grid on
% grid minor