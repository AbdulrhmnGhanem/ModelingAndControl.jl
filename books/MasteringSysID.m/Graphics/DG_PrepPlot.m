function DG_PrepPlot(Freq,Data,Type,Color,units,titleText,options);
%//
%// Plot the specified complex column vector(s) record versus frequency, in publication form
%//
%// ARGUMENTS:
%//		- Freq	:		Vector of frequencies
%//		- Type	:   'A'mplitude in dB,					=> (complex) y data are converted to dB and plotted
%//               'P'hase in degree,					=> complex y data are converted to angle in degree and plotted
%//               'M'agnitude (linear scale)	=> (complex) data are taken magnitude of and plotted
%//               'R'eal number plot					=> real numbers are plot on a linear scale
%//								'DScale' dBm scale (LSNA)		=> Scale is the number of time domain points in the acquired signal
%//
%//		- Color	:   'B' for black and white plot, 
%//               'C' for full line color change, 
%//               'G' for full line greyscale , 
%//               'L' for linestyle black
%//								'S' for symbol black
%//                 for extended: specify 1 style / trace {'TR1','TR2',}
%//                     each trace string consits of 3 cells:
%//											'C=COLOR'
%//												'BLA' = 'black','BLU' = 'blue','RED' = 'red','GRE' = 'green','MAG' = 'magenta','YEL' = 'yellow'
%//												default = color as by column number
%//											'S=Symbol'
%//												[ + | o | * | . | x | square | diamond | v | ^ | > | < | pentagram | hexagram | {none} ]
%//												default = no symbol
%//											'L=linetype'
%//												[ {-} | -- | : | -. | none ]
%//												default = full line
%//											'T=tint%'
%//												tint of the current trace. Default = 100%
%//         
%//		- units :   2 strings,	
%//                 + first for X,
%//                     leave empty for automatic adaptation with Hz as a unit,			=> ''
%//											specify units without brackets for automatic scaling				=> 'Length'
%//											bracketed units for direct print														=> 'Freq [Hz]'
%//									+ second for Y axis		:
%//		- titletext	:	Title of the plot
%//		- options		: special things to be asked for, cell array of strings
%//										'XMIN=value' : minimum value for the X-axis
%//										'XMAX=value' : maximum value for the X-axis
%//										'YMIN=value' : minimum value for the Y-axis
%//										'YMAX=value' : maximum value for the Y-axis


%// Work on the current figure on top
lFigureId = gcf;

if (nargin < 7), options = {};end;
if (nargin < 5), units{1} = ''; end;

%// Empty string means Hz
if (isempty(units{1}))
	units{1}	= 'Freq [Hz]';
end

%// Scale the frequency vector
if strfind(units{1},'[')
	decFreq			= max(floor(max(log10(abs(Freq(Freq~=0))))/3)*3);
	Freq				= Freq./(10.^decFreq);
	units{1}		= getUnitString(units{1}, decFreq);
else
	decFreq = 0;
end

if nargin < 6
	titleText = '';
end

decData = 0;

%// plot amplitude or phase
switch upper(Type(1))

	case 'D'
		
		lScale = str2num(Type(2:end));
		if isempty(units{2})
			units{2} = 'P [dBm]';
		end
		ScaledData = 10.0*log10(abs((2*Data/lScale).^2/50e-3));
		if strfind(units{2},'[')
			decData = max(floor(max(log10(abs(ScaledData)))/3)*3);
			ScaledData = ScaledData./(10.^decData);
		end
		if ((max(max(abs(ScaledData)))<1)||(max(max(abs(ScaledData)))>1000)),error('DG_PrepPlot:ranger error');end;
		plot(Freq, ScaledData)
		if (decData~=0)
			units{2}		= getUnitString(units{2}, decData);
		end

	case 'A'

		if isempty(units{2})
			units{2} = 'Mag [dB]';
		end
		ScaledData = 20.0*log10(abs(Data));
		if strfind(units{2},'[')
			decData = max(floor(max(log10(abs(ScaledData)))/3)*3);
			ScaledData = ScaledData./(10.^decData);
		end
		if ((max(max(abs(ScaledData)))<1)||(max(max(abs(ScaledData)))>1000)),error('DG_PrepPlot:ranger error');end;
		plot(Freq, ScaledData)
		if (decData~=0)
			units{2}		= getUnitString(units{2}, decData);
		end

	case 'M'

		if isempty(units{2})
			units{2} = 'Mag [V]';
		end
		ScaledData = abs(Data);

		if strfind(units{2},'[')
			%// If there is more than 1 value, take the maximum
			decData = max(floor(max(log10(ScaledData))/3)*3);
			ScaledData = ScaledData./(10.^decData);
		end
		if ((max(max(abs(ScaledData)))<1)||(max(max(abs(ScaledData)))>1000)),error('DG_PrepPlot:ranger error');end;
		plot(Freq, ScaledData)
		if (decData~=0)
			units{2}		= getUnitString(units{2}, decData);
		end
		


	case 'P'

		if isempty(units{2})
			units{2} = '\Phi [o]';
		end
		ScaledData = 180.0/pi*angle(Data);
		if strfind(units{2},'[')
			decData = max(floor(max(log10(abs(ScaledData)))/3)*3);
			ScaledData = ScaledData./(10.^decData);
		end
		if ((max(max(abs(ScaledData)))<1)||(max(max(abs(ScaledData)))>1000)),error('DG_PrepPlot:ranger error');end;
		plot(Freq, ScaledData)
		if (decData~=0)
			units{2}		= getUnitString(units{2}, decData);
		end

	case 'R'

		if isempty(units{2})
			units{2} = 'Voltage [V]';
		end
		ScaledData = Data;
		if strfind(units{2},'[')
			%// If there is more than 1 value, take the maximum
			decData = max(floor(max(log10(abs(ScaledData)))/3)*3);
			ScaledData = ScaledData./(10.^decData);
		end
		if ((max(max(abs(ScaledData)))<1)||(max(max(abs(ScaledData)))>1000)),error('DG_PrepPlot:ranger error');end;
		plot(Freq, ScaledData)
		if (decData~=0)
			units{2}		= getUnitString(units{2}, decData);
		end


end

if isempty(options)
	%// Axis and labels are selected automatically
	Limits = axis;
	%// Adapt the frequency range to meet strict limits
	Limits(1) = min(Freq);
	Limits(2) = max(Freq);
	minData		= Limits(1);
	maxData		= Limits(2);
	meanData	= (maxData+minData)/2;

	if (abs(maxData-minData)*100>abs(meanData))
		%// calculate the rounded size of 2 3 4 slices
		Delta			= (maxData-minData);
		Slice			= Delta ./[2;3;4];
		rndSlice	= round(Slice.*10.^(-floor(log10(Slice)))).*10.^(floor(log10(Slice)));
		%// Replace the slices on the axis
		minsliced = floor(minData./rndSlice).*rndSlice;
		maxsliced = ceil(maxData./rndSlice).*rndSlice;
		%// Actual number of slices
		ActualSlices = (maxsliced-minsliced)./rndSlice;
		%// Differences in slice units
		minDelta	= (minsliced-minData)./rndSlice;
		maxDelta	= (maxsliced-maxData)./rndSlice;
		%// Count the trials with minimal number of 1/2 empty minimum cells
		fullBounds	= (abs(minDelta)<0.5)+(abs(maxDelta)<0.5);
		sel					= find(fullBounds==max(fullBounds));
		%// If there is more than 1 combination left, sort on minimal freespace
		[dist,distSel] = sort(abs(minDelta(sel))+abs(maxDelta(sel)));
		%// Select the contribution
		sel				= sel(distSel(1));
		rndSlice	= rndSlice(sel);
		minsliced = minsliced(sel);
		maxsliced = maxsliced(sel);
		DG_SetLabelValue('X',minsliced:rndSlice:maxsliced,lFigureId);
		Limits(1) = minsliced;
		Limits(2) = maxsliced;
	else
		%// Fix the labels 
		DG_SetLabelFormat('X','auto',lFigureId)
		DG_SetLabelValue('X','auto',lFigureId)
	end

	minData		= min(min(ScaledData));
	maxData		= max(max(ScaledData));
	meanData	= (maxData+minData)/2;
	if (abs(maxData-minData)*100>abs(meanData))
		%// calculate the rounded size of 2 3 4 slices
		Delta			= (maxData-minData);
		Slice			= Delta ./[2;3;4];
		rndSlice	= round(Slice.*10.^(-floor(log10(Slice)))).*10.^(floor(log10(Slice)));
		%// Replace the slices on the axis
		minsliced = floor(minData./rndSlice).*rndSlice;
		maxsliced = ceil(maxData./rndSlice).*rndSlice;
		%// Actual number of slices
		ActualSlices = (maxsliced-minsliced)./rndSlice;
		%// Differences in slice units
		minDelta	= (minsliced-minData)./rndSlice;
		maxDelta	= (maxsliced-maxData)./rndSlice;
		%// Count the trials with minimal number of 1/2 empty minimum cells
		fullBounds	= (abs(minDelta)<0.5)+(abs(maxDelta)<0.5);
		sel					= find(fullBounds==max(fullBounds));
		%// If there is more than 1 combination left, sort on minimal freespace
		[dist,distSel] = sort(abs(minDelta(sel))+abs(maxDelta(sel)));
		%// Select the contribution
		sel				= sel(distSel(1));
		rndSlice	= rndSlice(sel);
		minsliced = minsliced(sel);
		maxsliced = maxsliced(sel);
		DG_SetLabelValue('Y',minsliced:rndSlice:maxsliced,lFigureId);
		Limits(3) = minsliced;
		Limits(4) = maxsliced;
	else
		%// Fix the labels 
		DG_SetLabelFormat('Y','auto',lFigureId)
		DG_SetLabelValue('Y','auto',lFigureId)
	end
else
	%// Some axis or label limits are set manually
	for optind = 1:length(options)
		lopt = options{optind};
		switch upper(lopt{1})
			case 'XMIN'
				Limits(1) = lopt{2}/(10^decFreq);
			case 'XMAX'
				Limits(2) = lopt{2}/(10^decFreq);
			case 'YMIN'
				Limits(3) = lopt{2}/(10^decData);
			case 'YMAX'
				Limits(4) = lopt{2}/(10^decData);
		end
	end
end

axis(Limits)

DG_SetFontName('HelveticaLTStd-Roman');
DG_SetFontSize(14)

ylabel(units{2})
if ~isempty(titleText),title(titleText);end

%// If there is only one trace, use a simple line
if iscell(Color)
	%// An extended call is used, read in the manual settings
	for lTraceInd = 1:size(Data,2)
		lAttrCell = Color{lTraceInd};
		for lAttrInd = 1:length(lAttrCell)
			lAttrVal = lAttrCell{lAttrInd};
			switch upper(lAttrVal(1))
				case 'C'
					DG_SetTraceColor(lAttrVal(3:5),lTraceInd,lFigureId);
				case 'L'
					DG_SetTraceStyle(lAttrVal(3:end),lTraceInd,lFigureId);
				case 'S'
					DG_SetSymbol(lAttrVal(3),lTraceInd,lFigureId);
				case 'T'
					DG_SetTraceTint(str2num(lAttrVal(3:end)),lTraceInd,lFigureId);
      end           
		end
	end
else
	%// An predefined call is used
	if size(Data,2)==1
		%// Color mode with one line => red trace
		switch upper(Color(1))
			case {'C'}
				DG_SetTraceColor('RED','*',lFigureId)			
			case {'B','G','L'}
				%// All other modes with one line => black trace
				DG_SetTraceColor('BLACK','*',lFigureId);	
			case {'S'}
				DG_SetTraceColor('BLACK','*',lFigureId);
				DG_SetSymbol('+','*',lFigureId);
		end
	else 
		%// There is more than 1 trace in the data
		switch upper(Color(1))
			
			case {'C'}
				%// A color plot is asked for, color for each trace
				for lTraceInd = 1:size(Data,2)
					switch (lTraceInd)
						case 1
							DG_SetTraceColor('RED',lTraceInd,lFigureId);
						case 2
							DG_SetTraceColor('GREEN',lTraceInd,lFigureId);
						case 3
							DG_SetTraceColor('BLUE',lTraceInd,lFigureId);
						case 4
							DG_SetTraceColor('CYAN',lTraceInd,lFigureId);
					end
				end
				
			case {'G'}
				%// A greyscale plot is asked for
				for lTraceInd = 1:size(Data,2)
					DG_SetTraceColor('BLACK',lTraceInd,lFigureId);
					switch (lTraceInd)
						case 1
							DG_SetTraceTint(100,lTraceInd,lFigureId);
						case 2
							DG_SetTraceTint(50,lTraceInd,lFigureId);
						case 3
							DG_SetTraceTint(30,lTraceInd,lFigureId);
						case 4
							DG_SetTraceTint(10,lTraceInd,lFigureId);
					end
				end
				
			case {'L'}
				for lTraceInd = 1:size(Data,2)
					DG_SetTraceColor('BLACK',lTraceInd,lFigureId);
					switch (lTraceInd)
						case 1
							DG_SetTraceStyle('-',lTraceInd,lFigureId);
						case 2
							DG_SetTraceStyle(':',lTraceInd,lFigureId);
						case 3
							DG_SetTraceStyle('--',lTraceInd,lFigureId);
						case 4
							DG_SetTraceStyle('-.',lTraceInd,lFigureId);
					end
				end

			case {'S'}
				for lTraceInd = 1:size(Data,2)
					DG_SetTraceColor('BLACK',lTraceInd,lFigureId);
					switch (lTraceInd)
						case 1
							DG_SetSymbol('+',lTraceInd,lFigureId);
						case 2
							DG_SetSymbol('o',lTraceInd,lFigureId);
						case 3
							DG_SetSymbol('x',lTraceInd,lFigureId);
						case 4
							DG_SetSymbol('*',lTraceInd,lFigureId);
					end
				end
				
			case {'B'}
				DG_SetTraceColor('BLACK','*',lFigureId);
				
		end
	end
end
%// Set the font size for the labels
xlabel(units{1})

%// Set the line width
DG_SetTraceWidth(2.0,'*',lFigureId)
DG_SetLineWidth(1.0,lFigureId)



%///////////////////////////////////////////////////////////////////////////////////////////////////////////
function unit = getUnitString(unit, decFactor)

switch decFactor
	case -3
		mul	= 'm';
	case 0
		mul	= '';
	case 3
		mul	= 'k';
	case 6
		mul	= 'M';
	case 9
		mul	= 'G';
	otherwise
		mul = ['10^',num2str(decFactor)];
end

unit = strrep(unit, '[',['[',mul]);
