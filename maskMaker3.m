function [out1, out2]=maskMaker3(dc,da)

% PURPOSE: to be a more adaptable mask making program for C-A.

% INPUT: dc = diameter of circle in pixels, odd numbers 5 - 9
% da = diameter of annulus in pixels, odd numbers 7 - 13.

% OUTPUT: masks for circle and annulus to use with ministacks to determine
% C-A or dF/S.

% Adapted from maskMaker.m (AW) by BB on 7/21/21

Z = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
O = [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0];
T = [0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0];
F = [0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0];
F1 = [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0];
F2 = [0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0];

R1 = [0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0];
R2 = [0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0];
R3 = [0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];
R4 = [0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0];
R5 = [0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0];
R6 = [0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0];
R7 = [0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0];
R8 = [0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0];
R9 = [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0];
R10 = [0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0];
R11 = [0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0];
R12 = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];


if dc == 5;
    mask13 = [Z; Z; Z; Z; Z; Z; Z; Z; Z; Z; O; T; F; T; O; Z; Z; Z; Z; Z; Z; Z; Z; Z; Z;];
elseif dc == 7
    mask13 = [Z; Z; Z; Z; Z; Z; Z; Z; Z; O; T; F; F1; F; T; O; Z; Z; Z; Z; Z; Z; Z; Z; Z;];
elseif dc == 9
    mask13 = [Z; Z; Z; Z; Z; Z; Z; Z; O; T; F; F1; F2; F1; F; T; O; Z; Z; Z; Z; Z; Z; Z; Z;];
end

if da == 7;
    ring13 = [Z; Z; Z; Z; Z; Z; Z; Z; R1; R2; R3; R4; R4; R4; R3; R2; R1; Z; Z; Z; Z; Z; Z; Z; Z;];
elseif da == 9;
    ring13 = [Z; Z; Z; Z; Z; Z; Z; R1; R2; R3; R4; R5; R5; R5; R4; R3; R2; R1; Z; Z; Z; Z; Z; Z; Z;];
elseif da == 11;
    ring13 = [Z; Z; Z; Z; Z; Z; R1; R2; R3; R4; R5; R6; R6; R6; R5; R4; R3; R2; R1; Z; Z; Z; Z; Z; Z;];
elseif da == 13; 
    ring13 = [Z; Z; Z; Z; Z; R1; R2; R3; R4; R5; R6; R7; R7; R7; R6; R5; R4; R3; R2; R1; Z; Z; Z; Z; Z;];
elseif da == 15
    ring13 = [Z; Z; Z; Z; R1; R2; R3; R4; R5; R6; R7; R8; R8; R8; R7; R6; R5; R4; R3; R2; R1; Z; Z; Z; Z;];
elseif da == 17
    ring13 = [Z; Z; Z; R1; R2; R3; R4; R5; R6; R7; R8; R9; R9; R9; R8; R7; R6; R5; R4; R3; R2; R1; Z; Z; Z;];
elseif da == 19
    ring13 = [Z; Z; R1; R2; R3; R4; R5; R6; R7; R8; R9; R10; R10; R10; R9; R8; R7; R6; R5; R4; R3; R2; R1; Z; Z;];
elseif da == 21
    ring13 = [Z; R1; R2; R3; R4; R5; R6; R7; R8; R9; R10; R11; R11; R11; R10; R9; R8; R7; R6; R5; R4; R3; R2; R1; Z;];
end

maskConverted = uint16(mask13);
ringConverted = uint16(ring13);


out1 = maskConverted;
out2 = ringConverted;
end