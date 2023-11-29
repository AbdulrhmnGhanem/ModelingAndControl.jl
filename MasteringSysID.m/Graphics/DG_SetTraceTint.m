function DG_SetTraceTint(Tint,Trace,FigNum,SubFigNum)
%//
%// Set the trace tint for a plot (0% to 100%)
%// Assumption: at least one color index has a color value of 1 for renormalisation
%//							or all the color indexesa are 0
%//
%// Parameters 
%//		- Color				: New tint in %
%//		- Trace				: index of the trace to be changed ('*' = all traces)
%//		- FigNum			: Pointer to the figure as returned by gcf
%//										if empty -> use the current figure
%//		- SubFigNum		: Index of the subplot as is used by subplot
%//										if empty  -> apply to current plot on figure
%//										if '*'    -> apply to all plots on figure

if any((Tint <= 0)|(Tint>100))
	error(' >> DG_SetTraceTint: Illegal tint value');
end


switch nargin
	case 1,
		error('  >> DG_SetTraceTint: At least a symbol and a trace are needed')

	case 2,
		FigNum			= gcf;
		SubFigNum		= gca;
		DoTraceTint(Tint,Trace,FigNum,SubFigNum)

	case 3,
		SubFigNum		= gca;
		DoTraceTint(Tint,Trace,FigNum,SubFigNum)

	case 4,
		ch = get(FigNum,'children');
		ch = flipud(ch(:));
		figure(FigNum)
		if isstr(SubFigNum)
			for ind_ch = 1:length(ch)
				DoTraceTint(Tint,Trace,FigNum,ch(ind_ch))
			end
		else
			if (SubFigNum <= length(ch))
				DoTraceTint(Tint,Trace,FigNum,ch(SubFigNum))
			else
				error(' >> DG_SetTraceTint: illegal number of subplots')
			end
		end
end

%///////////////////////////////////////////////////////////////////////
function DoTraceTint(Tint,Trace,FigNum,SubFigNum)

	figure(FigNum);
	AllTraces = get(SubFigNum,'Children');
	AllTraces = flipud(AllTraces(:));
  if isstr(Trace)
		if (length(Tint)==1)
			for ind_tr = 1:length(AllTraces)
				%// Read the old color and calculate the new tint
				ColorVal	= get(AllTraces(ind_tr),'Color');
				LUserData = get(AllTraces(ind_tr),'UserData');
				if isempty(LUserData)
					set(AllTraces(ind_tr),'UserData',ColorVal);
				else
					ColorVal = LUserData;
				end
				ColorVal = ChangeTint(ColorVal,Tint);
				%// Set the new 'color'
				set(AllTraces(ind_tr),'Color',ColorVal);
			end;
		elseif (length(Tint)==length(AllTraces))
			for ind_tr = 1:length(AllTraces)
				ColorVal	= get(AllTraces(ind_tr),'Color');
				LUserData = get(AllTraces(ind_tr),'UserData');
				if isempty(LUserData)
					set(AllTraces(ind_tr),'UserData',ColorVal);
				else
					ColorVal = LUserData;
				end
				%// Read the old color and calculate the new tint
				ColorVal = ChangeTint(ColorVal,Tint(ind_tr));
				%// Set the new 'color'
				set(AllTraces(ind_tr),'Color',ColorVal);
			end;
		else
			error(' >> DG_SetTraceTint : illegal trace specified');
		end
	elseif (Trace<=length(AllTraces))
		%// Read the old color and calculate the new tint
		ColorVal	= get(AllTraces(Trace),'Color');
		LUserData = get(AllTraces(Trace),'UserData');
		if isempty(LUserData)
			set(AllTraces(Trace),'UserData',ColorVal);
		else
			ColorVal = LUserData;
		end
		%// Set the new 'color'
		ColorVal	= ChangeTint(ColorVal,Tint);
		set(AllTraces(Trace),'Color',ColorVal);
	else
		error(' >> DG_SetTraceTint : illegal trace specified');
	end
 
%// Color definitions
function TintedColor = ChangeTint(TraceColor,Tint);
%// Renormalise the color and then take the tint requested

	%// Renormalise the color
	NormalisationFactor = max(TraceColor);
	%// If the color is black, take white s a reference
	if (NormalisationFactor == 0)
		TintedColor = [1 1 1]*(100-Tint)/100;
	else
		%//TintedColor = TraceColor/NormalisationFactor*Tint/100;
		TintedColor = TraceColor + (100-Tint)/100*[1 1 1];
		TintedColor(TintedColor>1) = 1;
	end
