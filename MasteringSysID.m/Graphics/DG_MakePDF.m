function  DG_MakePDF(FileName, Figure);
%//
%// Make a pdf file of the specified figure window under foxed circumstances

set(Figure, 'PaperPositionMode', 'auto');
set(Figure, 'PaperUnits', 'centimeters');

OuterPosition = get(get(Figure,'CurrentAxes'),'OuterPosition');

% if isempty(findstr(FileName,'.ps'))
%   FileName = [pwd, '.ps'];
% end

aa=dir('*.pdf');
if any(strcmpi(FileName,{aa.name}))
 delete(FileName);
end

if isempty(findstr(FileName,'\'))
  FileName = [pwd, '\', FileName];
end


%// Rescale to get 1.6 ratio
%// Left bottom width height
%//set(Figure, 'PaperPosition', OuterPosition);
print(Figure,'-dpdf',FileName)
disp(FileName)
eval(['!"C:\Program Files (x86)\Adobe\Acrobat 6.0\Acrobat\acrobat.exe" "',FileName,'" &']);
keyboard