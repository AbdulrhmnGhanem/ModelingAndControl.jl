function DG_SetSymbolSize(FontSize,Trace,FigNum,SubFigNum)
%//
%// Set the symbol size for a plot
%//
%// Parameters 
%//		- FontSize		: New font size for symbols on the trace
%//										if Trace = '*', an array of sizes can be given, one for each trace
%//		- Trace				: index of the trace to be changed ('*' = all traces)
%//		- FigNum			: Pointer to the figure as returned by gcf
%//										if empty -> use the current figure
%//		- SubFigNum		: Index of the subplot as is used by subplot
%//										if empty  -> apply to current plot on figure
%//										if '*'    -> apply to all plots on figure

switch nargin
	case 1,
		error('  >> DG_SetSymbolSize: At least a symbol and a trace are needed')

	case 2,
		FigNum			= gcf;
		SubFigNum		= gca;
		DoSymbolSize(FontSize,Trace,FigNum,SubFigNum)

	case 3,
		SubFigNum		= gca;
		DoSymbolSize(FontSize,Trace,FigNum,SubFigNum)

	case 4,
		ch = get(FigNum,'children');
		ch = flipud(ch(:));
		figure(FigNum)
		if isstr(SubFigNum)
			for ind_ch = 1:length(ch)
				DoSymbolSize(FontSize,Trace,FigNum,ch(ind_ch))
			end
		else
			if (SubFigNum <= length(ch))
				DoSymbolSize(FontSize,Trace,FigNum,ch(SubFigNum));
			else
				error(' >> DG_SetSymbolSize : illegal number of subplots');
			end
		end
end

%///////////////////////////////////////////////////////////////////////
function DoSymbolSize(FontSize,Trace,FigNum,SubFigNum)

	figure(FigNum);
	AllTraces = get(SubFigNum,'Children');
	AllTraces = flipud(AllTraces(:));
  if isstr(Trace)
		if (length(FontSize)==1)
			for ind_tr = 1:length(AllTraces)
				set(AllTraces(ind_tr),'MarkerSize',FontSize);
			end;
		elseif (length(FontSize)==length(AllTraces))
			for ind_tr = 1:length(AllTraces)
				set(AllTraces(ind_tr),'MarkerSize',FontSize(ind_tr));
			end;
		else
			error(' >> DG_SetSymbolSize : illegal trace specified');
		end
	elseif (Trace<=length(AllTraces))
		set(AllTraces(Trace),'MarkerSize',FontSize);
	else
		error(' >> DG_SetSymbolSize : illegal trace specified');
	end
