function DG_PlotSij(Freq,Sij,Type,Color);
%//
%// Plot the specified S-parameter record versus frequency, in publication form
%//
%// ARGUMENTS:
%//		- Type	: 'A'mplitude in dB or 'P'hase in degree
%//		- Color	: 'B'W for black and white, 'C'OLOR for full , 'G' for greyscale 

DG_PrepPlot(Freq,Sij,Type,Color)