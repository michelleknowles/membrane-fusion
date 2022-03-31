function stkwrite(movi,file,pth)
if nargin<3
   [file,pth] = uiputfile('*.stk','Where to put the Stack');
end

if length(file)<=4
   file = [file,'.stk'];
end
if ~strcmp(file(end-3),'.')
   file = [file,'.stk'];
end

if ~(isa(movi,'uint8')|isa(movi,'uint16'))
   disp('Input Must Be 8 or 16 bit');
   return
end
nframes = size(movi,3);
Width = size(movi,2);
Height = size(movi,1);
if isa(movi,'uint8')
   BitDepth = 8;
end
if isa(movi,'uint16')
   BitDepth = 16;
end
Compression = 1;%Compression = none;
PhotometricInterpretation = 1;
%calculate RowsPerStrip
%8000 ~= Width*RowsPerStrip*BitDepth/8
RowsPerStrip = round(8000*8/Width/BitDepth);
%unless RowsPerStrip is Larger Than Height;
RowsPerStrip = min(RowsPerStrip,Height);
stripsPerImage = floor((Height+RowsPerStrip -1)/RowsPerStrip);
StripByteCounts = BitDepth*Width*RowsPerStrip/8*ones(stripsPerImage,1);
StripByteCounts(end) = BitDepth*Width*Height/8 - ...
   sum(StripByteCounts(1:end-1));
StripOffsets = [8;cumsum(StripByteCounts)+8];
StripOffsets = StripOffsets(1:end-1);
if RowsPerStrip == Height
   StripByteCounts = BitDepth*Width*RowsPerStrip/8;
   StripOffsets = 8;
end
XResolution = 72;
YResolution = 72;
ResolutionUnit = 2;
Predictor = 1;
UIC2Tag = ones(1,6*nframes);
UIC3Tag = ones(1,2*nframes);
numberOfFields = 15;

fid = fopen([pth,file],'w+');
fwrite(fid,[0 0 0 0 0 0 0 0],'uint8');
h = waitbar(0,'Saving Stack');
for f = 1:nframes
   switch BitDepth
   case 8
      fwrite(fid,movi(:,:,f)','uint8');
   case 16
      fwrite(fid,movi(:,:,f)','uint16');
   end
   waitbar(f/nframes,h)
end
close(h)
pIFD = ftell(fid);

%fill in the file header
fseek(fid,0,'bof');
fwrite(fid,[73,73,42,0],'uint8');
fwrite(fid,pIFD,'uint32');

%go to the end of the movi and add the IFD
fseek(fid,pIFD,'bof');

fwrite(fid,numberOfFields,'uint16');

%fill in the fields
%field 254 newsubfiletype
fwrite(fid,254,'uint16');
%type long
fwrite(fid,4,'uint16');
%number 1
fwrite(fid,1,'uint32');
%value 0
fwrite(fid,0,'uint32');

%field 256 Width
fwrite(fid,256,'uint16');
%type long
fwrite(fid,4,'uint16');
%number 1
fwrite(fid,1,'uint32');
%value Height
fwrite(fid,Width,'uint32');

%field 257 Height
fwrite(fid,257,'uint16');
%type long
fwrite(fid,4,'uint16');
%number 1
fwrite(fid,1,'uint32');
%value Height
fwrite(fid,Height,'uint32');

%field 258 BitDepth
fwrite(fid,258,'uint16');
%type short
fwrite(fid,3,'uint16');
%number 1
fwrite(fid,1,'uint32');
%value Bitdepth
fwrite(fid,BitDepth,'uint32');

%field 259 Compression
fwrite(fid,259,'uint16');
%type short
fwrite(fid,3,'uint16');
%number 1
fwrite(fid,1,'uint32');
%value Compression
fwrite(fid,Compression,'uint32');

%field 262 PhotometricInterpretation
fwrite(fid,262,'uint16');
%type short
fwrite(fid,3,'uint16');
%number 1
fwrite(fid,1,'uint32');
%value PhotometricInterpretation
fwrite(fid,PhotometricInterpretation,'uint32');

%field 270 ? may need to write an image discription
%see metamorph, onee null terminated string for each plane
%in the stack

%field 273 StripOffsets
fwrite(fid,273,'uint16');
%type long
fwrite(fid,4,'uint16');
%number
if length(StripOffsets)>1
fwrite(fid,length(StripOffsets),'uint32');
%Write the Strip Offsets After the IFD
%Add the pointer to the strip offsets once they are written
%Use the location pStripOffsets as the place to write the pointer
ppStripOffsets = ftell(fid);
fwrite(fid,0,'uint32'); %skip forward 4bytes
else
   fwrite(fid,1,'uint32');
   fwrite(fid,StripOffsets,'uint32');
end


%field 278 RowsPerStrip
fwrite(fid,278,'uint16');
%type long
fwrite(fid,4,'uint16');
%number 1
fwrite(fid,1,'uint32');
%value RowsPerStrip
fwrite(fid,RowsPerStrip,'uint32');

%field 279 StripByteCounts
fwrite(fid,279,'uint16');
%type long
fwrite(fid,4,'uint16');
%number
if length(StripByteCounts)>1
fwrite(fid,length(StripByteCounts),'uint32');
ppStripByteCounts = ftell(fid);
fwrite(fid,0,'uint32'); %skip forward 4bytes
%see use of this variable above at pStripOffsets
else
   fwrite(fid,1,'uint32');
   fwrite(fid,StripByteCounts,'uint32');
end

%field 282 XResolution
fwrite(fid,282,'uint16');
%type Rational
fwrite(fid,5,'uint16');
%number 1
fwrite(fid,1,'uint32');
ppXResolution = ftell(fid);
fwrite(fid,0,'uint32'); %skip forward 4bytes

%field 283 YResoulution
fwrite(fid,283,'uint16');
%type Rational
fwrite(fid,5,'uint16');
%number 1
fwrite(fid,1,'uint32');
ppYResolution = ftell(fid);
fwrite(fid,0,'uint32'); %skip forward 4bytes

%field 296 Resolution uint
fwrite(fid,296,'uint16');
%type short
fwrite(fid,3,'uint16');
%number 1
fwrite(fid,1,'uint32');
%value Resolution uint
fwrite(fid,ResolutionUnit,'uint32');

%field 305 Software Could be include if needed
%not sure Weather to leave this field blank
%Indicate matlab's creation or MetaMorphs
%currently this tag is ignored

%field 306 Date ignored as well

%field 317 Predictor
fwrite(fid,317,'uint16');
%type short
fwrite(fid,3,'uint16');
%number 1
fwrite(fid,1,'uint32');
%value = Predictor
fwrite(fid,Predictor,'uint32');

%field 33628 UIC1Tag Ignored for now

%field 33629 UIC2Tag
fwrite(fid,33629,'uint16');
%type Rational
fwrite(fid,5,'uint16');
%number nframes
fwrite(fid,nframes,'uint32');
ppUIC2Tag = ftell(fid);
fwrite(fid,0,'uint32'); %skip forward 4bytes

%field 33630 UIC3Tag
fwrite(fid,33630,'uint16');
%type Rational
fwrite(fid,5,'uint16');
%number nframes
fwrite(fid,nframes,'uint32');
ppUIC3Tag = ftell(fid);
fwrite(fid,0,'uint32'); %skip forward 4bytes

%Write Four bytes of zeros
fwrite(fid,0,'uint32');

%now Write the tables
if length(StripOffsets)>1
   pStripOffsets = ftell(fid);
   fwrite(fid,StripOffsets,'uint32');
end
if length(StripByteCounts)>1
   pStripByteCounts = ftell(fid);
   fwrite(fid,StripByteCounts,'uint32');
end

pXResolution = ftell(fid);
fwrite(fid,[XResolution,1],'uint32');
pYResolution = ftell(fid);
fwrite(fid,[YResolution,1],'uint32');
pUIC2Tag = ftell(fid);
fwrite(fid,UIC2Tag,'uint32');
pUIC3Tag = ftell(fid);
fwrite(fid,UIC3Tag,'uint32');

if length(StripOffsets)>1
   fseek(fid,ppStripOffsets,'bof');
   fwrite(fid,pStripOffsets,'uint32');
end

if length(StripByteCounts)>1
   fseek(fid,ppStripByteCounts,'bof');
   fwrite(fid,pStripByteCounts,'uint32');
end

fseek(fid,ppXResolution,'bof');
fwrite(fid,pXResolution,'uint32');

fseek(fid,ppYResolution,'bof');
fwrite(fid,pYResolution,'uint32');

fseek(fid,ppUIC2Tag,'bof');
fwrite(fid,pUIC2Tag,'uint32');

fseek(fid,ppUIC3Tag,'bof');
fwrite(fid,pUIC3Tag,'uint32');

fclose(fid);
