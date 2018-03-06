function [WAKE, OUTP, INPU] = fcnINITIALIZE(COND, INPU)

% Preallocating for a turbo-boost in performance
OUTP.vecCL = nan(COND.valMAXTIME,INPU.valVEHICLES);
OUTP.vecCLF = nan(COND.valMAXTIME,INPU.valVEHICLES);
OUTP.vecCLI = nan(COND.valMAXTIME,INPU.valVEHICLES);
OUTP.vecCDI = nan(COND.valMAXTIME,INPU.valVEHICLES);
OUTP.vecE = nan(COND.valMAXTIME,INPU.valVEHICLES);
OUTP.vecCT = nan(COND.valMAXTIME,max(INPU.vecPANELROTOR));
OUTP.vecCP = nan(COND.valMAXTIME,max(INPU.vecPANELROTOR));
OUTP.vecCTCONV = nan(COND.valMAXTIME, max(INPU.vecPANELROTOR));

INPU.vecPANELROTOR = uint16(INPU.vecPANELROTOR);
INPU.vecPANELWING = uint16(INPU.vecPANELWING);
INPU.vecN = uint8(INPU.vecN);
INPU.vecM = uint8(INPU.vecM);


% Initializing wake parameters
WAKE.matWAKEGEOM = [];
WAKE.matNPWAKEGEOM = [];
WAKE.vecWDVEHVSPN = [];
WAKE.vecWDVEHVCRD = [];
WAKE.vecWDVEROLL = [];
WAKE.vecWDVEPITCH = [];
WAKE.vecWDVEYAW = [];
WAKE.vecWDVELESWP = [];
WAKE.vecWDVEMCSWP = [];
WAKE.vecWDVETESWP = [];
WAKE.vecWDVEAREA = [];
WAKE.matWDVENORM = [];
WAKE.matWVLST = [];
WAKE.matWDVE = uint32([]);
WAKE.valWNELE = 0;
WAKE.matWCENTER = [];
WAKE.matWCOEFF = [];
WAKE.vecWK = [];
WAKE.matWADJE = uint32([]);
WAKE.vecWDVEPANEL = uint16([]);
WAKE.valLENWADJE = 0;
WAKE.vecWKGAM = [];
WAKE.vecWDVESYM = uint8([]);
WAKE.vecWDVETIP = uint8([]);
WAKE.vecWDVESURFACE = uint8([]);
WAKE.vecWDVETRI = [];
WAKE.vecWPLOTSURF = uint8([]);

end
