function DG_SetTraceColor(Color,Trace,FigNum,SubFigNum)
%//
%// Set the trace color for a plot
%//
%// Parameters 
%//		- Color				: New color of the trace 'black','blue','red','green','magenta','yellow'
%//										if Trace = '*', a cell array of colors can be given, one color per trace
%//		- Trace				: index of the trace to be changed ('*' = all traces)
%//		- FigNum			: Pointer to the figure as returned by gcf
%//										if empty -> use the current figure
%//		- SubFigNum		: Index of the subplot as is used by subplot
%//										if empty  -> apply to current plot on figure
%//										if '*'    -> apply to all plots on figure

switch nargin
	case 1,
		error('  >> DG_SetTraceColor: At least a symbol and a trace are needed')

	case 2,
		FigNum			= gcf;
		SubFigNum		= gca;
		DoTraceColor(Color,Trace,FigNum,SubFigNum)

	case 3,
		SubFigNum		= gca;
		DoTraceColor(Color,Trace,FigNum,SubFigNum)

	case 4,
		ch = get(FigNum,'children');
		ch = flipud(ch(:));
		figure(FigNum)
		if isstr(SubFigNum)
			for ind_ch = 1:length(ch)
				DoTraceColor(Color,Trace,FigNum,ch(ind_ch))
			end
		else
			if (SubFigNum <= length(ch))
				DoTraceColor(Color,Trace,FigNum,ch(SubFigNum))
			else
				error(' >> DG_SetTraceColor : illegal number of subplots');
			end
		end
end

%/////////////////////////////////////////////////////////////////////////////////////
function DoTraceColor(Color,Trace,FigNum,SubFigNum)

	figure(FigNum);
	AllTraces = get(SubFigNum,'Children');
	AllTraces = flipud(AllTraces(:));
  if isstr(Trace)
		if (~iscell(Color))
			for ind_tr = 1:length(AllTraces)
				set(AllTraces(ind_tr),'Color',GetColorVal(Color));
				set(AllTraces(ind_tr),'UserData',[]);
			end;
		elseif (length(Color)==length(AllTraces))
			for ind_tr = 1:length(AllTraces)
				set(AllTraces(ind_tr),'Color',GetColorVal(Color{ind_tr}));
				set(AllTraces(ind_tr),'UserData',[]);
			end;
		else
			error(' >> DG_SetTraceColor : illegal trace specified');
		end
	elseif (Trace<=length(AllTraces))
		set(AllTraces(Trace),'Color',GetColorVal(Color));
		set(AllTraces(Trace),'UserData',[]);
	else
		error(' >> DG_SetTraceColor : illegal trace specified');
	end

%/////////////////////////////////////////////////////////////////////////////////////
%// Color definitions
function ColorVal = GetColorVal(Color)
	switch upper(Color(1:3))
		case 'BLU'
			ColorVal = [0 0 1];
		case 'RED'
			ColorVal = [1 0 0];
		case 'GRE'
			ColorVal = [0 1 0];
		case 'CYA'
			ColorVal = [0 1 1];
		case 'MAG'
			ColorVal = [0 1 1];
		case 'YEL'
			ColorVal = [1 1 0];
		case 'BLA'
			ColorVal = [0 0 0];
		otherwise
			ColorVal = [0 0 0];
	end
