function DG_SetFontSize(FontSize,FigNum,SubFigNum)
%//
%// Set the font size for a plot
%//
%// Parameters 
%//		- FontSize		: New size of the font
%//		- FigNum			: Pointer to the figure as treturned by gcf
%//										if empty -> use the current figure
%//		- SubFigNum		: Index the subplot as is used by subplot
%//										if empty  -> apply to current plot on figure
%//										if '*'    -> apply to all plots on figure

if nargin == 1
	FigNum		= gcf;
	SubFigNum = gca;
	DoFontSize(FontSize,FigNum,SubFigNum);
elseif (nargin == 2)
	SubFigNum = gca;
	DoFontSize(FontSize,FigNum,SubFigNum);
else
	AxisNum = get(gcf,'children');
	AxisNum = flipud(AxisNum(:));
	if isstr(SubFigNum)
		for ind_ch = 1:length(AxisNum)
			DoFontSize(FontSize,FigNum,AxisNum(ind_ch));
		end
	else
		if (SubFigNum)<=length(AxisNum)
			DoFontSize(FontSize,FigNum,AxisNum(SubFigNum));
		else
			error('  >> DG_SetFontName: illegal subplot');
		end
	end
end

%//////////////////////////////////////////////////////////////////////////////////////
function DoFontSize(FontSize,FigNum,SubFigNum)

	ch = get(FigNum,'children');
	figure(FigNum)
	theAxis = find(ch==SubFigNum);
	if (length(theAxis)==1)
		set(ch(theAxis),'FontSize',FontSize);
	else
		error(' >> DG_SetFontSize : illegal subplot specified');
	end
