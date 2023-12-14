%// Initialise the size of the figure
t = linspace(0,1,100);
plot(t,sin(2*pi*t),'k-',t,-sin(4*pi*t),'k--',t,cos(4*pi*t),'k-.')

DG_SetFontName('Helvetica')
DG_SetFontSize(14)
DG_SetLabelFormat('Y','%6.2f')
title('bla')
xlabel('blix')
ylabel('bly')
DG_SetTraceWidth(1.5,'*')
DG_SetTraceWidth(1.0,'*')
DG_SetTraceTint(30,1)

%// Initialise the size of the figure to get a fixed size pdf
DG_Init4PDF(gcf);
%// Save the file under pdf formaat
DG_MakePDF('test', gcf);

%// Initialise the size of the figure
t = linspace(0,1,100);
subplot(321);plot(t,sin(2*pi*t),'k-',t,-sin(4*pi*t),'k--',t,cos(4*pi*t),'k-.')
subplot(322);plot(t,sin(3*pi*t),'k-',t,-sin(6*pi*t),'k--',t,cos(6*pi*t),'k-.')
subplot(323);plot(t,sin(4*pi*t),'k-',t,-sin(8*pi*t),'k--',t,cos(8*pi*t),'k-.')
subplot(324);plot(t,sin(5*pi*t),'k-',t,-sin(10*pi*t),'k--',t,cos(10*pi*t),'k-.')
subplot(325);plot(t,sin(4*pi*t),'k-',t,-sin(8*pi*t),'k--',t,cos(8*pi*t),'k-.')
subplot(326);plot(t,sin(5*pi*t),'k-',t,-sin(10*pi*t),'k--',t,cos(10*pi*t),'k-.')


subplot(321);
DG_SetFontSize(14,gcf,1)
DG_SetFontName('HelveticaLTStd-Roman',gcf,1)
title('Plot1');
ylabel('kak_2')
xlabel('pis^2')
DG_SetLabelFormat('Y','%6.2f',gcf,1)
DG_SetLineWidth(2,gcf,1)
DG_SetSymbol('xo+','*',gcf,1)
DG_SetSymbolSize([12,10,8],'*',gcf,1)
DG_SetTraceStyle({'-','--',':'},'*',gcf,1)
DG_SetTraceTint([90,70,50],'*',gcf,1)
DG_SetTraceWidth([4,3,2],'*',gcf,1)

subplot(322);
DG_SetFontSize(14,gcf,2)
DG_SetFontName('HelveticaLTStd-Roman',gcf,2)
title('Plot2');
ylabel('kak_2')
xlabel('pis^2')
DG_SetTraceColor({'blue','red','green'},'*',gcf,2)

subplot(323)
DG_SetFontSize(14,gcf,3)
DG_SetFontName('HelveticaLTStd-Roman',gcf,3)
title('Plot3');
ylabel('kak_2')
xlabel('pis^2')

subplot(324)
DG_SetFontSize(14,gcf,4)
DG_SetFontName('HelveticaLTStd-Roman',gcf,4)
title('Plot4');
ylabel('kak_2')
xlabel('pis^2')

subplot(325)
DG_SetFontSize(14,gcf,5)
DG_SetFontName('HelveticaLTStd-Roman',gcf,5)
title('Plot5');
ylabel('kak_2')
xlabel('pis^2')

subplot(326)
DG_SetFontSize(14,gcf,6)
DG_SetFontName('HelveticaLTStd-Roman',gcf,6)
title('Plot6');
ylabel('kak_2')
xlabel('pis^2')


DG_Init4PDF(gcf,5);

%// Save the file under pdf formaat
DG_MakePDF('test2', gcf);
