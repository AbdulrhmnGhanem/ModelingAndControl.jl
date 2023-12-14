function DG_SetFontName(FontName,FigNum,SubFigNum)
%//
%// Set the font for a plot
%//
%// Parameters 
%//		- FontName		: New font to be used !!! WARNING !!! Font name has to be exact!
%//		- FigNum			: Pointer to the figure as returned by gcf
%//										if empty -> use the current figure and the current (sub)plot
%//		- SubFigNum		: Index of the subplot to be used, as numbered in the subplot command
%//										if empty  -> apply to current plot on figure
%//										if '*'    -> apply to all plots on figure

if nargin == 0
    FontName = 'HelveticaLTStd-Roman';
	FigNum		= gcf;
	SubFigNum = gca;
	DoFontName(FontName,FigNum,SubFigNum);

elseif nargin == 1
	FigNum		= gcf;
	SubFigNum = gca;
	DoFontName(FontName,FigNum,SubFigNum);
elseif (nargin == 2)
	SubFigNum = gca;
	DoFontName(FontName,FigNum,SubFigNum);
else
	AxisNum = get(gcf,'children');
	AxisNum = flipud(AxisNum(:));
	if isstr(SubFigNum)
		for ind_ch = 1:length(AxisNum)
			DoFontName(FontName,FigNum,AxisNum(ind_ch));
		end
	else
		if (SubFigNum)<=length(AxisNum)
			DoFontName(FontName,FigNum,AxisNum(SubFigNum));
		else
			error('  >> DG_SetFontName: illegal subplot');
		end
	end
end

%//////////////////////////////////////////////////////////////////////////////////////
function DoFontName(FontName,FigNum,SubFigNum)

	ch = get(FigNum,'children');
	figure(FigNum)
	theAxis = find(ch==SubFigNum);
	if (length(theAxis)==1)
		set(ch(theAxis),'FontName',FontName);
	else
		error(' >> DG_SetFontName : illegal subplot specified');
	end
