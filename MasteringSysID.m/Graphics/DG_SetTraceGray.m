function DG_SetTraceGray(lColor,Trace,FigNum)
%//
%// Set the trace width for a plot
%//
%// Parameters 
%//		- LineWidth		: New width of the trace
%//		- Trace			: index of the trace to be changed ('*' = all traces)
%//		- lColor		: color array or scalar for grayscale
%//		- FigNum		: if empty -> use the current plot
%//						  if figure number -> apply to all plots on figure
 
if (nargin == 1)
	Traces = get(get(gcf,'CurrentAxes'),'Children');
	for ind_tr = 1:length(Traces)
		if length(lColor) == 1
			set(Traces(ind_tr),'Color',lColor*[1 1 1]);
		elseif length(lColor) == 3
			set(Traces(ind_tr),'Color',lColor);
		end
	end;
elseif (nargin == 2)
	Traces = get(gca,'Children');
	if isstr(Trace)
		for ind_tr = 1:length(Traces)
			if length(lColor) == 1
				set(Traces(ind_tr),'Color',lColor*[1 1 1]);
			elseif length(lColor) == 3
				set(Traces(ind_tr),'Color',lColor);
			end
		end;
	elseif (Trace<=length(Traces))
		if length(lColor) == 1
			set(Traces(Trace),'Color',lColor*[1 1 1]);
		elseif length(lColor) == 3
			set(Traces(Trace),'Color',lColor);
		end
	end
elseif (nargin == 3)
	ch = get(gcf,'children');
	for ind_ch = 1:length(ch)
		figure(FigNum)
		Traces = get(ch(ind_ch),'Children');
		if isstr(Trace)
			for ind_tr = 1:length(Traces)
				if length(lColor) == 1
					set(Traces(ind_tr),'Color',lColor*[1 1 1]);
				elseif length(lColor) == 3
					set(Traces(ind_tr),'Color',lColor);
				end
			end;
		elseif (Trace<=length(Traces))
			if length(lColor) == 1
				set(Traces(Trace),'Color',lColor*[1 1 1]);
			elseif length(lColor) == 3
				set(Traces(Trace),'Color',lColor);
			end
		end
	end
end