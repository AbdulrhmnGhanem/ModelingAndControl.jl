function DG_SetFont(FontName,FigNum)
%//
%// Set the font size for a plot
%//
%// Parameters 
%//		- FontName			: New Font Name (Default = HelveticaLTStd-Roman)
%//		- FigNum			: if empty -> use the current plot
%//										if figure number -> apply to all plots on figure

switch nargin
    case 0,
        FontName = 'HelveticaLTStd-Roman';
        set(get(gcf,'CurrentAxes'),'FontName',FontName);
    case 1,
        set(get(gcf,'CurrentAxes'),'FontName',FontName);
    otherwise
        ch = get(gcf,'children');
        for ind_ch = 1:length(ch)
            figure(FigNum)
            set(ch(ind_ch),'FontName',FontName);
        end
end