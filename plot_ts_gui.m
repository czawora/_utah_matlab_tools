function [hf, ha, hl] = plot_ts_gui(X, varargin)
%PLOT_TS_GUI
%Written by: Julio I. Chapeton
%
% Plots multiple time series with the option to scroll or click to navigate through time.
% By default scrolling is used only if there are less than 1.5 x 10^7 total samples (i.e. numel(X) < 1.5e7). This corresponds to 2.5 minutes for 100 channels...
% with a 1000Hz sampling rate. This can be overridden using the 'scroll' option, however, there will be a lot of lag if the number of samples/channels is too large.
%
%Usage:
%PLOT_TS_GUI(X,varargin)
%PLOT_TS_GUI(...,'PropertyName',PropertyValue,...)
%
%Inputs:
%    Required:
%        X          = matrix  (Input array, this function operates columnwise, i.e. each column is a time series)
%
%    Optional:
%        (name-value pairs, not case sensitve)
%        scroll     = logical (Use this option to override the default scroll setting. There will be a lot of lag if there is too much data)
%        Fs         = scalar  (Sampling rate in Hz. If used the time axis will be in seconds)
%        y_labels   = cell    (Cell array of labels for the y axis, e.g. channel names)
%        npoints    = scalar  (If scroll mode is disabled, this option will define, in number of samples, how much data will be displayed at a time for each time series...
%                             The default is the minimum between the length of each time seires and 150,000 samples. i.e. 'min(size(X,1),1.5e5)' )
%        title      = string  (Add a plot title)
%        disp       = logical (provide a value of 1 to display the figure, or 0 to hide them...
%                              This is useful in the case that the figure is being saved, but the user does not want it displayed.)
%Outputs:
%    Required:
%        hf         = handle  (figure handle)
%        ha         = handle  (axis handle)
%        hl         = handle  (line handles)

%% parse inputs
inp_pars = inputParser;

[nsamples,nchanns]=size(X);
N=nsamples*nchanns;
if N <= 1.5e7
    def_scroll = 1;
else
    def_scroll = 0;
end
def_Fs = 1;
def_ylabels=num2cell([1:size(X,2)]);
def_npoints=min(nsamples-1,1.5e5);
def_title='';
def_disp = 1;

%def_saving = 0;
addParameter(inp_pars,'scroll',def_scroll,@(x) (x==1 || x==0));
addParameter(inp_pars,'Fs',def_Fs,@(x) (isscalar(x)));
addParameter(inp_pars,'y_labels',def_ylabels,@(x) assert(length(x)==numel(x) && length(x)==nchanns,'There are more labels than columns'));
addParameter(inp_pars,'npoints',def_npoints,@(x) (isscalar(x) && (x <= nsamples)) );
addParameter(inp_pars,'title',def_title,@(x) ischar(x) );
addParameter(inp_pars,'disp',def_disp,@(x) (x==1 || x==0));

parse(inp_pars,varargin{:})

scroll=inp_pars.Results.scroll;
Fs=inp_pars.Results.Fs;
y_labels=inp_pars.Results.y_labels;
npoints=inp_pars.Results.npoints;
plot_title=inp_pars.Results.title;
disp_ind=inp_pars.Results.disp;

if ~any(strcmp('npoints',inp_pars.UsingDefaults)) && ~any(strcmp('scroll',inp_pars.UsingDefaults)) && scroll
    warning('Both the ''scroll'' and ''npoints'' options have been turned on. The ''npoints'' option is only allowed if scrolling is off and will be ignored')
end
%%

if disp_ind
    hf=figure;
else
    hf=figure('visible','off');
end
%         set(h,'PaperOrientation','landscape'); % these commands are usefull if pdfs are going to be saved or if plots are going to be printed
%         set(h,'PaperUnits','normalized');
%         set(h,'PaperPosition', [0 0 1 1]);
ha=axes('Parent',hf);
set(ha,'LooseInset',get(ha,'TightInset'));

if any(strcmp(inp_pars.UsingDefaults,'Fs'))
    xlabel('sample number')%,ylabel('electrode');
    xscale = 1000;
else
    xlabel('time (s)')%,ylabel('electrode');
    xscale = 1;
end
set(ha,'YGrid','off')
set(ha,'XGrid','off')
title(plot_title)
drawnow;
switch scroll
    case 0
        %         hl=plot([0,1],ones(size(X,2),2));
        set(hf,'units','normalized','outerposition',[1 1 1 1]);
        set(ha,'LooseInset',get(ha,'TightInset'));
        btn1 = uicontrol('Style', 'pushbutton', 'String', 'next','Callback', @next_callback,'units','normalized','position',[0.95 0.95 0.025 0.025],'Fontsize',16,'Interruptible','on','UserData',0);
        btn2 = uicontrol('Style', 'pushbutton', 'String', 'previous','Callback',@previous_callback ,'units','normalized','position',[0.925 0.95 0.025 0.025],'Fontsize',16,'Interruptible','on','UserData',0);
        btn3 = uicontrol('Style', 'pushbutton', 'String', 'exit','Callback',@close_exit ,'units','normalized','position',[0.9 0.95 0.025 0.025],'Fontsize',16,'Interruptible','on','UserData',0);
        drawnow;
        
        %         ts=nan(npoints+1,size(X,2)); % should check for segments that are all nan or zeros
        %ts=nan(npoints,size(X,2)); % should check for segments that are all nan or zeros
        
        ind=1;
        while ind < size(X,1)
            ts=nan(npoints,size(X,2)); % should check for segments that are all nan or zeros
            if (ind+npoints)>size(X,1)
                ts(1:size(X,1)-ind+1,:)=X(ind:end,:);
                ts(size(X,1)-ind+1:end,:)=NaN;
            else
                ts=X(ind:ind+npoints,:);
            end
            %             rm_chans=all(isnan(ts)) | all(ts==0); % maybe check for rows that are constant?
            rm_chans=all(isnan(ts) | ts==0);
            ts(:,rm_chans)=[];
            
            if isempty(ts)
                warning('This segment consists of all NaNs or all 0s')
                commandwindow;
                inp=input('input 1 to go forward, -1 to go back, and 0 to exit');
                if inp==1
                    if ind + npoints <= size(X,1)
                        ind=ind+npoints;
                    else
                        beep
                        warning('You have reached the end of the clip')
                    end
                elseif inp==-1
                    if ind - npoints >= 1
                        ind=ind-npoints;
                    else
                        beep
                        warning('You have reached the beginning of the clip')
                    end
                elseif ~inp
                    return
                end
                continue
            end
            
            if disp_ind
                figure(hf)
            end
            hl=plot(ha,[0,1],ones(size(ts,2),2));
            y_labels_temp=y_labels(~rm_chans);
            
            T_secs = [ind-1:(ind-1)+(size(ts,1)-1)]./Fs;
            ub = max(ts); lb = min(ts);
            Sig=nanstd(ts);
            spacing=cumsum([0, abs(ub(1:end-1))+abs(lb(2:end))]);%+[0;Sig(1:end-1)]);
            spacing = repmat(spacing,numel(T_secs),1);
            ts_spaced=ts+spacing;
            ts_spaced=ts_spaced-min(ts_spaced(:));
            
%             set(hl,'xdata',T_secs,{'ydata'},num2cell(ts_spaced,1)','color','k','linewidth',1)
            set(hl,'xdata',T_secs,{'ydata'},num2cell(ts_spaced,1)','linewidth',1)
            %             temp=nanmean(ts_spaced)-min(ts_spaced(:));
            temp=nanmean(ts_spaced);
            
            % set ticks/labels
            %set(ha,'yaxislocation','right','ytick',temp,'yticklabel',y_labels,'ylim',[min(ts_spaced(:))-Sig(1) max(ts_spaced(:))+Sig(end)])
            % set ticks/labels
            if all(diff(temp))
                set(ha,'yaxislocation','right','ytick',temp(~isnan(temp)),'yticklabel',y_labels_temp(~isnan(temp)),'ylim',[min(ts_spaced(:))-Sig(1) max(ts_spaced(:))+Sig(end)])
            else
                % in case there are columns with all zeros add a small shift
                run_starts = [0; find(diff(temp)~=0)'] + 1;
                run_lengths = [diff(run_starts); numel(temp) - run_starts(end) + 1];
                
                rem_starts=run_starts(run_lengths>1);
                rem_ends=rem_starts+run_lengths(run_lengths>1);
                rem_lengths=run_lengths(run_lengths>1) - 1;
                for k=1:length(rem_starts)
                    temp(rem_starts(k)+1:rem_ends(k)-1)=temp(rem_starts(k)+1:rem_ends(k)-1)+cumsum(ones(1,rem_lengths(k))*mean(diff(temp))/length(temp));
                end
                set(ha,'yaxislocation','right','ytick',temp(~isnan(temp)),'yticklabel',y_labels_temp(~isnan(temp)),'ylim',[min(ts_spaced(:))-Sig(1) max(ts_spaced(:))+Sig(end)])
            end
            
            %             set(ha,'xtick',[0:15*xscale:max(T_secs)],'xticklabel',[0:15*xscale:max(T_secs)],'Ticklength', [0 0],'xlim',[min(T_secs) max(T_secs)],'xgrid','on','xminorgrid','on')
            set(ha,'xtick',[0:15*xscale:max(T_secs)],'xticklabel',[0:15*xscale:max(T_secs)],'Ticklength', [0 0],'xlim',[min(T_secs) max(T_secs)])
            
            drawnow;
            if npoints==size(X,1)
                delete(btn1)
                delete(btn2)
                delete(btn3)
                return
            end
            
            % enable moving through plots using buttons **
            uiwait
            try
                b1=get(btn1,'UserData');
                b2=get(btn2,'UserData');
            catch
                %warning('Figure was closed without the exit button')
                return
            end
            bb=b1+b2; % one should always be zero
            
            set(btn1,'UserData',0); % reset buttons
            set(btn2,'UserData',0);
            
            if get(btn3,'UserData')
                %close all
                delete(btn1)
                delete(btn2)
                delete(btn3)
                return %*** There is a known major bug using 'return' with 'parfor' and/or anonymous functions R2015b and R2016a
                %break
            elseif (ind+npoints*bb) >= 1 && (ind+npoints*bb)<size(X,1)
                ind=ind+npoints*bb;
            elseif (ind+npoints*bb) < 1
                %                 ind=1;
                beep
                warning('You are at the beginning of the clip')
            elseif (ind+npoints*bb)>=size(X,1)
                %                 ind=size(X,1)-npoints*bb-1;
                beep
                warning('You have reached the end of the clip')
            end % **
        end %while
        
    case 1
        
        ts=X;
        rm_chans=all(isnan(ts)) | all(ts==0);% maybe check for rows that are constant?
        ts(:,rm_chans)=[];
        if isempty(ts)
            warning('Data consists of all NaNs or all 0s')
            return
        end
        T_secs = [0:size(ts,1)-1]./Fs;
        ub = max(ts); lb = min(ts);
        Sig=nanstd(ts);
        spacing=cumsum([0, abs(ub(1:end-1))+abs(lb(2:end))]);%+[0;Sig(1:end-1)]);
        spacing = repmat(spacing,numel(T_secs),1);
        ts_spaced=ts+spacing;
        ts_spaced=ts_spaced-min(ts_spaced(:));
        
        hl=plot([0,1],ones(size(ts,2),2));
        set(hl,'xdata',T_secs,{'ydata'},num2cell(ts_spaced,1)')
        temp=nanmean(ts_spaced);
        ylim([min(ts_spaced(:),[],1,'omitnan')-Sig(1) max(ts_spaced(:),[],1,'omitnan')+Sig(end)])
        xlim([min(T_secs) max(T_secs)])
        
        % make figure scrollable
        hold on
        warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
        warning('off','MATLAB:modes:mode:InvalidPropertySet');
        yscroll=min(size(ts,2),10);
        scrollplot(ha,'MinY',min(ts_spaced(:))-Sig(1),'MaxY',max(ts_spaced(:,yscroll))+Sig(yscroll),'MaxX',min(30*xscale,T_secs(end)),'Axis','xy');
        set(ha,'LooseInset',get(ha,'TightInset'))
        
        % setting outerposition messes up scrollplot, this acomplishes the same thing but may be deprecated in future matlab releases
        jFrame = get(handle(hf), 'JavaFrame');
        drawnow;
        jFrame.setMaximized(1);
        
        % set ticks/labels
        if all(diff(temp))
            set(ha,'yaxislocation','right','ytick',temp(~isnan(temp)),'yticklabel',y_labels(~isnan(temp)))
        else
            % in case there are columns with all zeros add a small shift
            run_starts = [0; find(diff(temp)~=0)'] + 1;
            run_lengths = [diff(run_starts); numel(temp) - run_starts(end) + 1];
            
            rem_starts=run_starts(run_lengths>1);
            rem_ends=rem_starts+run_lengths(run_lengths>1);
            rem_lengths=run_lengths(run_lengths>1) - 1;
            for k=1:length(rem_starts)
                temp(rem_starts(k)+1:rem_ends(k)-1)=temp(rem_starts(k)+1:rem_ends(k)-1)+cumsum(ones(1,rem_lengths(k))*mean(diff(temp))/length(temp));
            end
            set(ha,'yaxislocation','right','ytick',temp(~isnan(temp)),'yticklabel',y_labels(~isnan(temp)))
        end
        set(ha,'xtick',[0:15*xscale:max(T_secs)],'xticklabel',[0:15*xscale:max(T_secs)],'Ticklength', [0 0]) % in minutes
        set(ha,'YGrid','on')
        set(ha,'XGrid','on')
        %grid on; box on
        %title(str)
        %xlabel('time in seconds')%,ylabel('electrode');
end
end
%% callback functions for ui objects
function close_exit(src,event)
drawnow
uiresume
set(src,'UserData',1)
end

function next_callback(src,event)
drawnow
uiresume
set(src,'UserData',1)
end

function previous_callback(src,event)
drawnow
uiresume
set(src,'UserData',-1)
end


%% OLD
% scroll=1;
% Fs=1000;
% str='title';
% npoints=min(size(X,1)-1,2.5*60*Fs); % If there is too much data then 2.5 minute chunks will be plotted and you can click the 'previous' or 'next' buttons to navigate
%
% figure
% hf=gcf;
% %         set(h,'PaperOrientation','landscape'); % these commands are usefull if pdfs are going to be saved or if plots are going to be printed
% %         set(h,'PaperUnits','normalized');
% %         set(h,'PaperPosition', [0 0 1 1]);
% ha=gca;
% set(ha,'LooseInset',get(gca,'TightInset'))
% hl=plot([0,1],ones(size(X,2),2));
% xlabel('time in seconds')%,ylabel('electrode');
% set(gca,'YGrid','on')
% set(gca,'XGrid','on')
% drawnow;
% switch scroll
%     case 0
%
%         set(hf,'units','normalized','outerposition',[1 1 1 1]);
%         set(gca,'LooseInset',get(gca,'TightInset'));
%         btn1 = uicontrol('Style', 'pushbutton', 'String', 'next','Callback', @next_callback,'units','normalized','position',[0.95 0.95 0.025 0.025],'Fontsize',16,'Interruptible','on','UserData',0);
%         btn2 = uicontrol('Style', 'pushbutton', 'String', 'previous','Callback',@previous_callback ,'units','normalized','position',[0.925 0.95 0.025 0.025],'Fontsize',16,'Interruptible','on','UserData',0);
%         btn3 = uicontrol('Style', 'pushbutton', 'String', 'exit','Callback',@close_exit ,'units','normalized','position',[0.9 0.95 0.025 0.025],'Fontsize',16,'Interruptible','on','UserData',0);
%         drawnow;
%
%         ind=1;
%         while ind < size(X,1)
%             if (ind+npoints)>=size(X,1)
%                 npoints=size(X,1)-ind;
%             end
%
%             ts=X(ind:ind+npoints,:)';
%             T_secs = [ind-1:(ind-1)+(size(ts,2)-1)]./Fs;
%             ub = max(ts,[],2); lb = min(ts,[],2);
%             Sig=nanstd(ts,0,2);
%             spacing=cumsum([0; abs(ub(1:end-1))+abs(lb(2:end))]);%+[0;Sig(1:end-1)]);
%             spacing = repmat(spacing,1,numel(T_secs));
%             ts_spaced=ts+spacing;
%
%             set(hl,'xdata',T_secs,{'ydata'},num2cell(ts_spaced,2))
%             temp=nanmean(ts_spaced,2);
%
%             % set ticks/labels
%             %             set(gca,'yaxislocation','right','ytick',temp,'yticklabel',tal_list(chan_ind),'ylim',[min(ts_spaced(:))-Sig(1) max(ts_spaced(:))+Sig(end)])
%             set(gca,'yaxislocation','right','ytick',temp,'ylim',[min(ts_spaced(:))-Sig(1) max(ts_spaced(:))+Sig(end)])
%             set(gca,'xtick',[0:5:max(T_secs)],'xticklabel',[0:5:max(T_secs)],'Ticklength', [0 0],'xlim',[min(T_secs) max(T_secs)]) % in minutes
%             %grid on; box on
%             title(str)
%             drawnow;
%
%             % enable moving through plots using buttons **
%             uiwait
%             b1=get(btn1,'UserData');
%             b2=get(btn2,'UserData');
%             bb=b1+b2; % one should always be zero
%
%             set(btn1,'UserData',0); % reset buttons
%             set(btn2,'UserData',0);
%
%             if get(btn3,'UserData')
%                 close all
%                 return %*** There is a known major bug using 'return' with 'parfor' and/or anonymous functions R2015b and R2016a
%                 %break
%             elseif (ind+npoints*bb) >= 1 && (ind+npoints*bb)<size(X,1)
%                 ind=ind+npoints*bb;
%             elseif (ind+npoints*bb) < 1
%                 %ind=1;
%                 beep
%                 warning('You are at the beginning of the clip')
%             elseif (ind+npoints*bb)>=size(X,1)
%                 %ind=ind-npoints*bb;
%                 beep
%                 warning('You have reached the end of the clip')
%             end % **
%         end %while
%
%     case 1
%
%         ts=X';
%         T_secs = [0:size(ts,2)-1]./Fs;
%         % T_secs = [ind-1:(ind-1)+(size(ts,2)-1)]./Fs;
%         ub = max(ts,[],2); lb = min(ts,[],2);
%         Sig=nanstd(ts,0,2);
%         spacing=cumsum([0; abs(ub(1:end-1))+abs(lb(2:end))]);%+[0;Sig(1:end-1)]);
%         spacing = repmat(spacing,1,numel(T_secs));
%         ts_spaced=ts+spacing;
%
%         set(hl,'xdata',T_secs,{'ydata'},num2cell(ts_spaced,2))
%         temp=nanmean(ts_spaced,2);
%         ylim([min(ts_spaced(:))-Sig(1) max(ts_spaced(:))+Sig(end)])
%         xlim([min(T_secs) max(T_secs)])
%
%         % make figure scrollable
%         hold on
%         warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
%         warning('off','MATLAB:modes:mode:InvalidPropertySet');
%         yscroll=min(size(X,2),10);
%         %scrollplot(h,'MinY',min(eegspaced(:))-Sig(1),'MaxY',max(eegspaced(yscroll,:))+Sig(yscroll),'MinX',0,'MaxX',min(300,T_secs(end)),'Axis','xy');
%         scrollplot(gca,'MinY',min(ts_spaced(:))-Sig(1),'MaxY',max(ts_spaced(yscroll,:))+Sig(yscroll),'MaxX',min(300,T_secs(end)),'Axis','xy');
%         set(gca,'LooseInset',get(gca,'TightInset'))
%
%         % setting outerposition messes up scrollplot, this acomplishes the same thing but may be deprecated in future matlab releases
%         jFrame = get(handle(hf), 'JavaFrame');
%         drawnow;
%         jFrame.setMaximized(1);
%
%         % set ticks/labels
%         %         set(gca,'yaxislocation','right','ytick',temp,'yticklabel',tal_list(chan_ind))
%         set(gca,'yaxislocation','right','ytick',temp)
%         set(gca,'xtick',[0:5:max(T_secs)],'xticklabel',[0:5:max(T_secs)],'Ticklength', [0 0]) % in minutes
%         set(gca,'YGrid','on')
%         set(gca,'XGrid','on')
%         %grid on; box on
%         title(str)
%         xlabel('time in seconds')%,ylabel('electrode');
% end

%                 % To enable moving through plots using arrows, replace the block between ** with the code below. Not as smooth.
%                 [~,~,usr_inp] = ginput(1);
%                 while isempty(ismember(usr_inp,[27,28,29]))
%                     [~,~,usr_inp] = ginput(1);
%                 end
%                 % right arrow
%                 % Use this for left click
%                 % if b == 1
%                 if usr_inp == 29
%                     if (ind+npoints) <= size(X,1)
%                         ind=ind+npoints; %// Move to the right
%                     elseif (ind+npoints)>size(X,1)
%                         %close(setdiff(get(0,'children'),[50,150]))
%                         break
%                     end
%                     % Left arrow
%                     % Use this for right click
%                     % elseif b == 3
%                 elseif usr_inp == 28
%                     if (ind-npoints) >= 1 %// Again check for out of bounds
%                         ind=ind-npoints; %// Move to the left
%                     end
%                     % Check for escape
%                 elseif usr_inp == 27
%                     %close(setdiff(get(0,'children'),[])) % if you want to keep a set of previous figures open then fill in [] with an array of figure handels
%                     break;
%                 end
