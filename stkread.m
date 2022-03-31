function out = stkread(file,pth,align)
if nargin<1
   [file,pth] = uigetfile('*.*','Choose the stack');
end
if ~file,return,end
if length(file)<=4
   file = [file,'.stk'];
end
if ~strcmp(file(end-3),'.')
   file = [file,'.stk'];
end

%get the information from the stack needed

fid = fopen([pth,file]);
if fid == -1
   fprintf('Can''t open file')
   return
end
if fread(fid,4,'uint8')' ~= [73 73 42 0];
   fprintf('Not a stack')
   return
end
tagpos = fread(fid,1,'uint32');
fseek(fid,tagpos,'bof');
num_tags = fread(fid,1,'uint16');
ps = (0:num_tags-1)*12 + tagpos + 2;
tags = [];
for p = ps
   fseek(fid,p,'bof');
   tags(end+1) = fread(fid,1,'uint16');
end

fseek(fid,ps(tags==256)+8,'bof');
tag.Width = fread(fid,1,'uint32');
fseek(fid,ps(tags==257)+8,'bof');
tag.Height = fread(fid,1,'uint32');
fseek(fid,ps(tags==258)+8,'bof');
tag.BitDepth = fread(fid,1,'uint32');
fseek(fid,ps(tags==33629)+4,'bof');
nframes = fread(fid,1,'uint32');
%nframes = floor(tagpos*8/tag.Width/tag.Height/tag.BitDepth);
cdims = [tag.Width,tag.Height];
mdims = [tag.Height,tag.Width];
ImageCount = prod(cdims);
M = zeros(mdims);
switch tag.BitDepth;
case 8
   M=uint8(M);								
   temp = uint8(zeros(tag.Width,tag.Height));
case 16
   M = uint16(M);
   temp = uint16(zeros(tag.Width,tag.Height));
end
if nargin ~=3
   out=M(:,:,ones(1,nframes));
else
   out = M(:,:,ones(1,nframes+1));
end
fseek(fid,8,'bof');
h = waitbar(0,'Loading Movie');
for f = 1:nframes
   switch tag.BitDepth;
   case 8
      temp(:) = uint8(fread(fid,ImageCount,'uint8'));
   case 16
      temp(:) = uint16(fread(fid,ImageCount,'uint16'));
   end
   out(:,:,f) = temp';
   waitbar(f/nframes,h)
end
close(h)
fclose(fid);
if nargin ==3
   if size(out,1) == size(align,1) & size(out,2) == size(align,2)
      out(:,:,end) = align;
   else
      'movie and alignment image are not the same size'
      out(:,:,end) = [];
   end
end
