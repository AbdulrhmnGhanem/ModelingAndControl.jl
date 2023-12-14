function DG_SetLabelValue(Axis,Value,FigNum, SubFigNum)
%//
%// Set the font size for a plot
%//
%// Parameters 
%//		- Axis				: the axis to be used 'X','Y','Z'
%//		- Value				: 'auto' or array of label values (numbers in increasing order, no strings,) 
%//                           Empty = auto
%//		- FigNum			: Pointer to the figure as returned by gcf
%//										if empty -> use the current figure
%//		- SubFigNum		: Index of the subplot as is used by subplot
%//										if empty  -> apply to current plot on figure
%//										if '*'    -> apply to all plots on figure
%//
%//     Example: DG_SetLabelValue('Y', [-120, -80, -40], 1, 2)

Axis = upper(Axis);
if (Axis(1)>='X')&(Axis(1)<='Z')
	
	switch nargin
		case 1,
			Value			= 'auto';
			FigNum		= gcf;
			SubFigNum = gca;
			figure(FigNum);
			SetLabelForOne(SubFigNum,Axis,Value);
			set()

		case 2,
			FigNum		= gcf;
			SubFigNum = gca;
			figure(FigNum);
			SetLabelForOne(SubFigNum,Axis,Value)

		case 3,
			SubFigNum = gca;
			figure(FigNum)
			SetLabelForOne(SubFigNum,Axis,Value)

		case 4,
			ch = get(gcf,'children');
			ch = flipud(ch(:));
			figure(FigNum)
			if ischar(SubFigNum)
				for ind_ch = 1:length(ch)
					SetLabelForOne(ch(ind_ch),Axis,Value);
				end
			else
				if (SubFigNum < length(ch))
					SetLabelForOne(ch(SubFigNum),Axis,Value);
				else
					error(' >> DG_SetLabelValue : illegal number of subplots');
				end
			end
	end
else
	error([' >> Invalid axis in DG_SetLabelValue',Axis])
end

%/////////////////////////////////////////////////////////////////////////////////////////////////
function SetLabelForOne(cFig,Axis,Value)

%// Built up axis specific labels
LabelMode = [Axis(1),'TickMode'];
LabelVal	= [Axis(1),'Tick'];


%// Check if valid Value
if isstr(Value)
	%// This is either an automatic spec or an error
	Value = upper(Value);
	if strcmp(Value,'AUTO')
		set(cFig,LabelMode,'auto');
	elseif strcmp(Value,'MAN')
		set(cFig,LabelMode,'man');
	else
		error([' >> Invalid value in DG_SetLabelValue',Value]);
	end

else
	%// This is a manual spec with a set of label values specified
	%// If the label value has been fixed, unlock it
	set(cFig,[Axis(1),'TickLabelMode'],'auto')

	set(cFig,LabelMode,'man');
	set(cFig,LabelVal,Value);
	
	LUserData = get(cFig,'UserData');
	switch Axis
		case 'X',
			LUserData.XLabelValue = Value;
		case 'Y',
			LUserData.YLabelValue = Value;
		case 'Z',
			LUserData.ZLabelValue = Value;
	end		
	set(cFig,'UserData',LUserData);

end
