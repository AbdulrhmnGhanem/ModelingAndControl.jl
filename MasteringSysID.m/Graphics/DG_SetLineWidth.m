function DG_SetLineWidth(LineWidth,FigNum,SubFigNum)
%//
%// Set the font size for a plot
%//
%// Parameters 
%//		- LineWidth		: New width of the line
%//		- FigNum			: Pointer to the figure as returned by gcf
%//										if empty -> use the current figure
%//		- SubFigNum		: Index of the subplot as is used by subplot
%//										if empty  -> apply to current plot on figure
%//										if '*'    -> apply to all plots on figure

switch nargin
	case 1,
		FigNum		= gcf;
		SubFigNum = gca;
		DoLineWidth(LineWidth,FigNum,SubFigNum);

	case 2,
		SubFigNum = gca;
		DoLineWidth(LineWidth,FigNum,SubFigNum);

	case 3,
		AxisNum = get(gcf,'children');
		AxisNum = flipud(AxisNum(:));
		if isstr(SubFigNum)
			for ind_ch = 1:length(AxisNum)
				DoLineWidth(LineWidth,FigNum,AxisNum(ind_ch));
			end
		else
			if (SubFigNum <= length(length(AxisNum))) 
				DoLineWidth(LineWidth,FigNum,AxisNum(SubFigNum));
			else
				error(' >> DG_SetLineWidth : illegal number of subplots');
			end
		end
end

%//////////////////////////////////////////////////////////////////////////////////////
function DoLineWidth(LineWidth,FigNum,SubFigNum)

	ch = get(FigNum,'children');
	figure(FigNum)
	theAxis = find(ch==SubFigNum);
	if (length(theAxis)==1)
		set(ch(theAxis),'LineWidth',LineWidth)
	else
		error(' >> DG_SetLineWidth : illegal subplot specified');
	end
