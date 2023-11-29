function DG_SetTraceWidth(LineWidth,Trace,FigNum,SubFigNum)
%//
%// Set the trace width for a plot
%//
%// Parameters 
%//		- LineWidth		: New width of the trace
%//		- Trace				: index of the trace to be changed ('*' = all traces)
%//		- FigNum			: Pointer to the figure as returned by gcf
%//										if empty -> use the current figure
%//		- SubFigNum		: Index of the subplot as is used by subplot
%//										if empty  -> apply to current plot on figure
%//										if '*'    -> apply to all plots on figure

switch nargin
	case 1,
		error('  >> DG_SetTraceWidth: At least a symbol and a trace are needed')

	case 2,
		FigNum			= gcf;
		SubFigNum		= gca;
		DoTraceWidth(LineWidth,Trace,FigNum,SubFigNum)

	case 3,
		SubFigNum		= gca;
		DoTraceWidth(LineWidth,Trace,FigNum,SubFigNum)

	case 4,
		ch = get(FigNum,'children');
		ch = flipud(ch(:));
		figure(FigNum)
		if isstr(SubFigNum)
			for ind_ch = 1:length(ch)
				DoTraceWidth(LineWidth,Trace,FigNum,ch(ind_ch))
			end
		else
			if (SubFigNum <= length(ch))
				DoTraceWidth(LineWidth,Trace,FigNum,ch(SubFigNum))
			else
				error('  >> DG_SetTraceWidth: illegal number of subplots')
			end
		end
end

%///////////////////////////////////////////////////////////////////////
function DoTraceWidth(LineWidth,Trace,FigNum,SubFigNum)

	figure(FigNum);
	AllTraces = get(SubFigNum,'Children');
	AllTraces = flipud(AllTraces(:));
  if isstr(Trace)
		if (length(LineWidth)==1)
			for ind_tr = 1:length(AllTraces)
				set(AllTraces(ind_tr),'LineWidth',LineWidth);
			end;
		elseif (length(LineWidth)==length(AllTraces))
			for ind_tr = 1:length(AllTraces)
				set(AllTraces(ind_tr),'LineWidth',LineWidth(ind_tr));
			end;
		else
			error(' >> DG_SetSymbolSize : illegal trace specified');
		end
	elseif (Trace<=length(AllTraces))
		set(AllTraces(Trace),'LineWidth',LineWidth);
	else
		error(' >> DG_SetSymbolSize : illegal trace specified');
	end
 