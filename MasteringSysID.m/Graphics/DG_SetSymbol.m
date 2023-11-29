function DG_SetSymbol(Symbol,Trace,FigNum,SubFigNum)
%//
%// Set the symbol size for a plot
%//
%// Parameters 
%//		- Symbol			: [ + | o | * | . | x | square | diamond | v | ^ | > | < | pentagram | hexagram | {none} ]
%//										if Trace = '*', a list of symbols can be given to apply to each trace successively
%//		- Trace				: index of the trace to be changed ('*' = all traces)
%//		- FigNum			: Pointer to the figure as returned by gcf
%//										if empty -> use the current figure
%//		- SubFigNum		: Index of the subplot as is used by subplot
%//										if empty  -> apply to current plot on figure
%//										if '*'    -> apply to all plots on figure

switch nargin
	case 1,
		error('  >> DG_SetSymbol: At least a symbol and a trace are needed')

	case 2,
		FigNum			= gcf;
		SubFigNum		= gca;
		DoSymbol(Symbol,Trace,FigNum,SubFigNum)

	case 3,
		SubFigNum		= gca;
		DoSymbol(Symbol,Trace,FigNum,SubFigNum)

	case 4,
		ch = get(FigNum,'children');
		ch = flipud(ch(:));
		figure(FigNum)
		if isstr(SubFigNum)
			for ind_ch = 1:length(ch)
				DoSymbol(Symbol,Trace,FigNum,ch(ind_ch));
			end
		else
			if (SubFigNum < length(ch))
				DoSymbol(Symbol,Trace,FigNum,ch(SubFigNum));
			else
				error(' >> DG_SetLabelValue : illegal number of subplots');
			end
		end
end

function DoSymbol(Symbol,Trace,FigNum,SubFigNum)

	figure(FigNum);
	AllTraces = get(SubFigNum,'Children');
	AllTraces = flipud(AllTraces(:));
  if isstr(Trace)
		if (length(Symbol)==1)
			for ind_tr = 1:length(AllTraces)
				set(AllTraces(ind_tr),'Marker',Symbol);
			end;
		elseif (length(Symbol)==length(AllTraces))
			for ind_tr = 1:length(AllTraces)
				set(AllTraces(ind_tr),'Marker',Symbol(ind_tr));
			end;
		else
			error(' >> DG_SetSymbol : illegal trace specified');
		end
	elseif (Trace<=length(AllTraces))
		set(AllTraces(Trace),'Marker',Symbol);
	else
		error(' >> DG_SetSymbol : illegal trace specified');
	end
