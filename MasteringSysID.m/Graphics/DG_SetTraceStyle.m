function DG_SetTraceStyle(TraceStyle,Trace,FigNum,SubFigNum)
%//
%// Set the trace style for a plot
%//
%// Parameters 
%//		- TraceStyle	: [ {-} | -- | : | -. | none ]
%//		- Trace				: index of the trace to be changed ('*' = all traces)
%//		- FigNum			: Pointer to the figure as returned by gcf
%//										if empty -> use the current figure
%//		- SubFigNum		: Index of the subplot as is used by subplot
%//										if empty  -> apply to current plot on figure
%//										if '*'    -> apply to all plots on figure

switch nargin
	case 1,
		error('  >> DG_SetTraceStyle: At least a symbol and a trace are needed')

	case 2,
		FigNum			= gcf;
		SubFigNum		= gca;
		DoTraceStyle(TraceStyle,Trace,FigNum,SubFigNum)

	case 3,
		SubFigNum		= gca;
		DoTraceStyle(TraceStyle,Trace,FigNum,SubFigNum)

	case 4,
		ch = get(FigNum,'children');
		ch = flipud(ch(:));
		figure(FigNum)
		if isstr(SubFigNum)
			for ind_ch = 1:length(ch)
				DoTraceStyle(TraceStyle,Trace,FigNum,ch(ind_ch))
			end
		else
			if (SubFigNum <= length(ch))
				DoTraceStyle(TraceStyle,Trace,FigNum,ch(SubFigNum))
			else
				error('  >> DG_SetTraceStyle: illegal subplot number given')
			end
		end
end

%/////////////////////////////////////////////////////////////////////////////////////
function DoTraceStyle(TraceStyle,Trace,FigNum,SubFigNum)

	figure(FigNum);
	AllTraces = get(SubFigNum,'Children');
	AllTraces = flipud(AllTraces(:));
  if isstr(Trace)
		if (~iscell(TraceStyle))
			for ind_tr = 1:length(AllTraces)
				set(AllTraces(ind_tr),'LineStyle',TraceStyle);
			end;
		elseif (length(TraceStyle)==length(AllTraces))
			for ind_tr = 1:length(AllTraces)
				set(AllTraces(ind_tr),'LineStyle',TraceStyle{ind_tr});
			end;
		else
			error(' >> DG_SetTraceStyle : illegal trace specified');
		end
	elseif (Trace<=length(AllTraces))
		set(AllTraces(Trace),'LineStyle',TraceStyle);
	else
		error(' >> DG_SetTraceStyle : illegal trace specified');
	end
