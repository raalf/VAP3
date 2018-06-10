function [ ] = fcnFIG2LATEX( varargin )
%SAVEFIG2LATEX Save Matlab figure to EPS or PDF for latex document
%   saveFig2Latex(FILENAME) saves current active figure as 6"x4" pdf to
%   default folder location 'Proposal/figures/' with FILENAME.
%
%   saveFig2Latex(FILENAME, [WIDTH, HEIGHT]) saves current active figure
%   with specific WIDTH and HEIGHT
%
%   saveFig2Latex(HFIG, FILENAME) saves HFIG
%
%   saveFig2Latex(HFIG, FILENAME, [WIDTH, HEIGHT]) saves HFIG with specific
%   WIDTH and HEIGHT
%
%   saveFig2Latex(FILENAME, PATH) saves current active figure as 6"x4" pdf to
%   default folder location 'Proposal/figures/' with FILENAME.
%
%   saveFig2Latex(FILENAME, PATH, [WIDTH, HEIGHT]) saves current active figure
%   with specific WIDTH and HEIGHT
%
%   saveFig2Latex(HFIG, FILENAME, PATH) saves HFIG
%
%   saveFig2Latex(HFIG, FILENAME, PATH, [WIDTH, HEIGHT]) saves HFIG with specific
%   WIDTH and HEIGHT


%
PATH = ''; % default PATH
WH = [6,4]; % default WIDTH and HEIGHT

% Determine if HFIG exist
if isa(varargin{1},'matlab.ui.Figure') % saveFig2Latex(HFIG, ...)
    hFig = varargin{1};
    varargin = varargin(2:end); % take hFig out of varargin
else
    hFig = gcf;% figure is not parsed, run gcf
end
%
% now varargin{1} must hold the FILENAME
if isa(varargin{1},'char')
    FILENAME = varargin{1};
else
    warning('FILENAME not found.');
end
%
% check if PATH, [WIDTH,HEIGHT] exist in varargin
% if both of them exist, varargin should have length of 3
if length(varargin) == 3
    PATH = varargin{2}; %replace default PATH
    WH = varargin{3}; %replace default WIDTH HEIGHT
elseif length(varargin) == 2
    WH = varargin{2};
end
%


%% Actual function starts here:


% Verify Inputs

% check if WIDTH and HEIGHT are stored in 1x2 array
if min(size(WH) == [1,2]) ~= 1
    WH = [6,4];% if the array is not 1x2
    warning('SIZE input error. Using default size.');
end
% check if PATH has '/' at the end
if ~isempty(PATH) && strcmp(PATH(end),'/') ~= 1 % if not
%     disp('no')
    PATH = strcat(PATH,'/'); %add '/' at the end of PATH
end

% assemble full file path
fullfilepath = strcat(PATH,FILENAME);




set(hFig,'Renderer','painters'); % ensure vector output
set(hFig,'Units','Inches');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',WH)

% '-fillpage' will ignore the aspect ratio of the axis
print(hFig,fullfilepath,'-dpdf','-r0','-fillpage');




end

