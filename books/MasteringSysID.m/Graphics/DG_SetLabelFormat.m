function DG_SetLabelFormat(Axis,Value,FigNum,SubFigNum)
%//
%// Set the format of the axis labels for the specified axis
%//
%// Parameters 
%//		- Axis				: the axis to be used 'X','Y','Z'
%//		- Value				: 'auto', 'man' to freeze the currrent situation
%//										or format specifier (starts with %, example '%4.1e').
%//										if empty -> use the 'auto' setting
%//		- FigNum			: Pointer to the figure as returned by gcf
%//										if empty -> use the current figure
%//		- SubFigNum		: Index of the subplot as is used by subplot
%//										if empty  -> apply to current plot on figure
%//										if '*'    -> apply to all plots on figure

Axis = upper(Axis);
if (Axis(1)>='X')&(Axis(1)<='Z')
	
	switch nargin
		case 1,
			Value			= 'auto';
			FigNum		= gcf;
			SubFigNum = gca;
			figure(FigNum)
			SetFormatForOne(SubFigNum,Axis,Value);

		case 2,
			FigNum		= gcf;
			SubFigNum = gca;
			figure(FigNum)
			SetFormatForOne(SubFigNum,Axis,Value);

		case 3,
			SubFigNum = gca;
			figure(FigNum)
			SetFormatForOne(SubFigNum,Axis,Value);

		case 4,
			ch = get(FigNum,'children');
			ch = flipud(ch(:));
			if isstr(SubFigNum)
				for ind_ch = 1:length(ch)
					SetFormatForOne(ch(ind_ch),Axis,Value);
				end
			else
				if (SubFigNum <= length(ch))
					SetFormatForOne(ch(SubFigNum),Axis,Value);
				else
					error(' >> DG_SetLabelFormat : illegal number of subplots');
				end
			end
	end
else
	error([' >> Invalid axis in DG_SetLabelFormat',Axis])
end

%/////////////////////////////////////////////////////////////////////////////////////////////////
function SetFormatForOne(cFig,Axis,Value)

%// Built up axis specific labels
LabelMode     = [Axis(1),'TickLabelMode'];
LabelVal	    = [Axis(1),'TickLabel'];

%// Check if valid Value
if isstr(Value)
	%// This is either an automatic spec or an error
	if strcmpi(Value,'AUTO')
		set(cFig,LabelMode,'auto');
	elseif strcmpi(Value,'MAN')
		set(cFig,LabelMode,'man');
	elseif (Value(1)=='%')
		%// This is a manual spec with a set of label values specified
		set(cFig,LabelMode,'man');
		labelValue	= get(cFig,LabelVal);
		labels			= cell(size(labelValue,1),1);
		for ind = 1:length(labels)
            if iscell(labelValue)
                labels{ind} = sprintf(Value,str2num(labelValue{ind}));
            else
                labels{ind} = sprintf(Value,str2num(labelValue(ind,:)));
            end
		end
		set(cFig,LabelVal,labels);
	else
		error([' >> Invalid value in DG_SetLabelFormat',Value]);
	end
else
		error([' >> Invalid value type in DG_SetLabelFormat, no string found']);
end

LUserData = get(cFig,'UserData');
switch Axis
	case 'X',
		LUserData.XFormatValue = Value;
	case 'Y',
		LUserData.YFormatValue = Value;
	case 'Z',
		LUserData.ZFormatValue = Value;
end		
set(cFig,'UserData',LUserData);
