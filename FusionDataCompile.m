function FusionDataCompile
%I want this function to detect all .xls files in a folder, open them,
%separate out data and write it to another .xls file.

%There should be a fusion sheet and a protein sheet.

%File names of the source xls sheet should be in the first row. Time on the
%left most column


%ABOUT THIS SCRIPT
%This script is to be used with data recieved from
%MiniStackFusionDataProcessV2. The data is in xls sheets with A being time
%in frames corresponding to the fusion event, B is the fusion data
%intensity, C is the protein data intensity.

%AW 1/12/22

%Gather Data About Files
myFolder = pwd;
filePattern = fullfile(myFolder, '*.xls');
xlsFiles   = dir(filePattern);
numberOfFiles = size(xlsFiles);
numberOfFiles = numberOfFiles(1,1);


%Initial xlsx setup
xlsxFileWriteName = 'CompiledFusionData.xlsx';
writematrix('Time Fusion', xlsxFileWriteName,'Sheet',1,'Range','A1');
writematrix('Time Protein', xlsxFileWriteName,'Sheet',2,'Range','A1');


timeArray = (-25:474);
%DO SOME FRAME TO TIME CONVERSION HERE


timeArray = transpose(timeArray);
writematrix(timeArray, xlsxFileWriteName,'Sheet',1,'Range','A2');
writematrix(timeArray, xlsxFileWriteName,'Sheet',2,'Range','A2');


%WRITE TIME TO COLUMN A MAKE SURE TO OFFSET THE REST OF THE DATA SO IT
%DOESNT OVER WRITE A.

for k = 1:length(xlsFiles)
  baseFileName = xlsFiles(k).name;
  %fullFileName = fullfile(myFolder, baseFileName);
  %fprintf('Now reading %s\n', fullFileName);
  
  %split this and write it into an excel or xls sheet
  tempData = xlsread(baseFileName);
  tempFusionData = tempData(:,2);
  tempProteinData = tempData(:,3);
  
  
  %Get int to xlsx cells
  k + 1
  
  activeXlsxCell = NumberToXLSColumnConverter(k + 1)
  activeXlsxCellName = append(activeXlsxCell,'1');
  activeXlsxCellData = append(activeXlsxCell,'2');
  
  
  %APPEND FILE NAME ABOVE DATA
  %Just write it to the first row in the column... too much work to append
  %a string to an array
  writematrix(baseFileName, xlsxFileWriteName,'Sheet',1,'Range',activeXlsxCellName);
  writematrix(baseFileName, xlsxFileWriteName,'Sheet',2,'Range',activeXlsxCellName);

  %WRITE DATA TO EXCEL SHEET 1 FUSION DATA SHEET 2 PROTEIN
  writematrix(tempFusionData, xlsxFileWriteName,'Sheet',1,'Range',activeXlsxCellData);
  writematrix(tempProteinData, xlsxFileWriteName,'Sheet',2,'Range',activeXlsxCellData);
  
  
end  
end



%ADDITIONAL FUNCTIONS

%NumberToXLSColumnConverter
function out = NumberToXLSColumnConverter(int)
%max int 26 then start appending another letter
alphabetString = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

%I want to build an array based on the int number. 
%NOTE max xls columns is 256...

%TripleLetters
if int >= 703
    
    disp('Over 702 files, please put some into a new folder')
    return;
else

    %Double Letters
    if (int > 26) && (int < 703) 
        dividedBy26 = floor(int/26); 
        remainder = int - 26*dividedBy26;
    
        %When getting to Z something odd happens
        if remainder == 0
            remainder = 26;
            dividedBy26 = dividedBy26 - 1;        
        end
    
        intArray = [dividedBy26 remainder];
    
    else
        %SingleLetters
        %number is going to be a single Letter
        intArray = [int];
    end
end

%OUTPUT
letters = alphabetString(intArray);
out = letters;
end