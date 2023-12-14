function  DG_Init4PDF(Figure, Width, Height)
%//
%// Prepare the figure at the final size for the processing as a pdf file
%//
%// INPUT 
%//  - Figure	: current figure handle (as obtained by gcf)
%//  - Width	: width of the figure in cm (Default = 7 cm)
%//  - Height	: Height of the figure in cm (Default = Width /1.618)


if (nargin <= 1), Width = 7; end;
if (nargin <= 2), Height = Width/1.618; end;

AllAxes = get(Figure,'children');
figure(Figure)
for ind_Ax = 1:length(AllAxes)
	set(AllAxes(ind_Ax),'ActivePositionProperty','Position')
	set(AllAxes(ind_Ax),'Units','centimeters');
	Position(:,ind_Ax) = get(AllAxes(ind_Ax),'Position');
	%// Scale first
	CurrPosition(1) = Position(1,ind_Ax);
	CurrPosition(2) = Position(2,ind_Ax);
	CurrPosition(3) = Width;
	CurrPosition(4) = Height;
	set(AllAxes(ind_Ax),'Position',CurrPosition);
	OuterPosition(:,ind_Ax) = get(AllAxes(ind_Ax),'OuterPosition');
	Position(:,ind_Ax) = get(AllAxes(ind_Ax),'Position');
end

%////////
%// Determine the position of the columns
ColPos = sort(Position(1,:));
ColPos = ColPos(diff([ColPos(1)-1,ColPos])~=0);
%// Determine the position of the rows
RowPos = sort(Position(2,:));
RowPos = RowPos(diff([RowPos(1)-1,RowPos])~=0);
%// Determine the offset between the plot rectangle and the bound of the plot
DeltaHor = OuterPosition(1,:)-Position(1,:);
DeltaVer = OuterPosition(2,:)-Position(2,:);

%// Get the spacing of the columns
for ColInd = 1:length(ColPos)
	%// Determine the maximal horizontal offset in a column
	SelAllInCol							= (Position(1,:) == ColPos(ColInd));
	DeltaHorFront(ColInd)		= max(Position(1,SelAllInCol)-OuterPosition(1,SelAllInCol));
	DeltaHorBack(ColInd)		= max((OuterPosition(1,SelAllInCol)+OuterPosition(3,SelAllInCol))-...
																(Position(1,SelAllInCol)+Position(3,SelAllInCol)));
end	
	
%// Get the spacing of the rows
for RowInd = 1:length(RowPos)
	%// Determine the maximal vertical offset in a row
	SelAllInRow							= (Position(2,:) == RowPos(RowInd));
	DeltaHorTop(RowInd)			= max(Position(2,SelAllInRow)-OuterPosition(2,SelAllInRow));
	DeltaHorBottom(RowInd)	= max((OuterPosition(2,SelAllInRow)+OuterPosition(4,SelAllInRow))-...
																(Position(2,SelAllInRow)+Position(4,SelAllInRow)));
end
	
for ColInd = 1:length(ColPos)
	SelAllInCol							= (Position(1,:) == ColPos(ColInd));
	%// Origin = endpoint of previous column + maximal horizontal offset in the column
	if (ColInd == 1)
		HorOriginofColumn = DeltaHorFront(ColInd);
	else
		%// Retrieve new position of previous graph
		HorOriginofColumn = sum(DeltaHorFront(1:ColInd))+sum(DeltaHorBack(1:ColInd-1))+sum(Position(3,1:ColInd-1)*1.05);
	end
	for RowInd = 1:length(RowPos)
		%// Select the element to move
		SelAllInRow = (Position(2,:) == RowPos(RowInd));
		Sel					= (SelAllInCol& SelAllInRow);
		if (RowInd == 1)
			VertOriginofRow = DeltaHorTop(RowInd);
		else
			%// Retrieve new position of previous graph
			VertOriginofRow = sum(DeltaHorTop(1:RowInd))+sum(DeltaHorBottom(1:RowInd-1))+sum(Position(4,1:RowInd-1)+Position(3,1:RowInd-1)*0.05);
		end
		%// Move it
		RCurrPosition(1) = HorOriginofColumn;
		RCurrPosition(2) = VertOriginofRow;
		RCurrPosition(3) = Position(3,1);
		RCurrPosition(4) = Position(4,1);
		set(AllAxes(Sel),'Position',RCurrPosition);
	end
end

%// Redraw the labels
for ind_Ax = 1:length(AllAxes)
	if strncmpi(get(AllAxes(ind_Ax),'YTickLabelMode'),'man',3)
		set(AllAxes(ind_Ax),'YTickLabelMode','auto');
		drawnow
		set(AllAxes(ind_Ax),'YTickLabelMode','man');
		UserData = get(AllAxes(ind_Ax),'UserData');
		labelValue	= get(AllAxes(ind_Ax),'YTickLabel');
		if ~isempty(UserData)
			if isfield(UserData,'YLabelValue')
				labelValue = UserData.YLabelValue;
			end
			if isfield(UserData, 'YFormatValue')
				labels			= cell(size(labelValue,1),1);
				for ind = 1:length(labels)
					labels{ind} = sprintf(UserData.YFormatValue,str2num(labelValue(ind,:)));
				end
				set(AllAxes(ind_Ax),'YTickLabel',labels);
			end
		end
	end
	
	if strncmpi(get(AllAxes(ind_Ax),'XTickLabelMode'),'man',3)
		set(AllAxes(ind_Ax),'XTickLabelMode','auto');
		set(AllAxes(ind_Ax),'XTickLabelMode','man');
		UserData = get(AllAxes(ind_Ax),'UserData');
		labelValue	= get(AllAxes(ind_Ax),'XTickLabel');
		if ~isempty(UserData)
			if isfield(UserData,'XLabelValue')
				labelValue = UserData.XLabelValue;
			end
			if isfield(UserData, 'XFormatValue')
				labels			= cell(size(labelValue,1),1);
				for ind = 1:length(labels)
					labels{ind} = sprintf(UserData.XFormatValue,str2num(labelValue(ind,:)));
				end
				set(AllAxes(ind_Ax),'XTickLabel',labels);
			end
		end
	end
	
	if strncmpi(get(AllAxes(ind_Ax),'ZTickLabelMode'),'man',3)
		set(AllAxes(ind_Ax),'ZTickLabelMode','auto');
		set(AllAxes(ind_Ax),'ZTickLabelMode','man');
		UserData = get(AllAxes(ind_Ax),'UserData');
		labelValue	= get(AllAxes(ind_Ax),'ZTickLabel');
		if ~isempty(UserData)
			if isfield(UserData,'ZLabelValue')
				labelValue = UserData.ZLabelValue;
			end
			if isfield(UserData, 'ZFormatValue')
				labels			= cell(size(labelValue,1),1);
				for ind = 1:length(labels)
					labels{ind} = sprintf(UserData.ZFormatValue,str2num(labelValue(ind,:)));
				end
				set(AllAxes(ind_Ax),'ZTickLabel',labels);
			end
		end
	end
end




%// Resize the window
%// This is a user thing

disp(' Please adapt the window to show ALL parts of the plots');
figure(Figure)
keyboard