function OptimumTE=simContrastvTE(T2pair,varargin)
%% usage function OptimumTE=simContrastvTE(T2pair,varargin)
% T2pair should be a vector with 2 T2 values
%
% this function plots the signal as a function of the echo time for two
% samples with different T2s 
% additionally it calculates the contrast as a function of TE and finds the
% optimum echo time


%% check if the correct variables exist
if ~exist('T2pair','var') 
    disp('ERROR: T2pair must be specified')
    return
end

if length(T2pair)~=2
    disp('ERROR: T2pair must contain 2 values')
    return
end

if ~exist('tColours','var')
    tColours = {[66 122 223]/255,[108 158 80]/255,[223 76 76]/255};
end
if nargin>1
    ShowFigure=varargin{1};
else
    ShowFigure=1;
end;

%% compute signal evolution as function of echo time

TE_vals = linspace(0,200,1000);
sig = zeros(length(TE_vals),2);
for iT2 = 1:2
    sig(:,iT2) = (exp(-TE_vals/T2pair(iT2)));
end
contrast=abs(sig(:,1)-sig(:,2));
[maximum position]=max(contrast);
OptimumTE=TE_vals(position);




%% plot results
if ShowFigure
figure
% set(gcf,'Position',[    50   679   740   409])
set(gcf,'Position',[ 75   306   400   400])
set(gcf,'Color',[1 1 1]);

plot(TE_vals,sig(:,1),'linewidth',2,'color',tColours{1}); hold all
plot(TE_vals,sig(:,2),'linewidth',2,'color',tColours{2}); 
h1 = plot(TE_vals,contrast);
set(h1,'linewidth',2,'color','k')
xlabel('TE (ms)')
ylabel('Relative signal')
grid on
legend(['T2 = ' num2str(T2pair(1)) 'ms'], ['T2 = ' num2str(T2pair(2)) 'ms'], 'Contrast');
text(OptimumTE ,0.5,['Optimum TE=',num2str(round(OptimumTE)),'ms'])

fontScale(1.2)
end;