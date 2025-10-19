%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving an Assignment, Project or Exam
% ======================================
%
% Write a script that populates the Q structs with
% your answers. This script is an example of that.
% 
% Once you are done ...
%   1) Run the solution script (this script).
%   2) Run the associated xxSubmit.p file to 
%      grade your answers and generate a MAT file.
%
% Always include a comment above any Matlab script.
% Matlab automatically uses it as a help screen.
%
% Try it. To see this comment in Matalb, enter:
%   "help apeExample"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A blank line terminates the help screen so
% These lines do not show up in the help screen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize

% Clean up workspace
clear all; clc;

% SN variable must contain Student Number
% This must be right or solution will not be graded
SN    = 12345678;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculations

% Q1: Scalars
% Do not round scalars
% Keep answers as accurate as possible
Q1.K1 = 123;            % integer
Q1.K2 = 1e2;            % integer with exponent
Q1.K3 = 1.23;           % floating point
Q1.K4 = 1.23e-6;        % floating point with exponent
Q1.K5 = 1+j*2;          % complex

% Q2: Vectors
% May be specified as row or column vectors
% Length must be right
Q2.Va = [1 2 3];                  % 1x3 row vec
Q2.Vb = [1 2 3]';                 % Use transpose operator for col vec
Q2.Vc = [1                        % Another way to make a 3x1 col vec
         2
         3];

% Q3: Matrices
% Both dimensions must be right
Q3.M1 = [1 2                      % 2x3 matrix
         3 4
         5 6];
Q3.M2 = [1 2;3 4];                % Use ';' to separate rows
Q3.M3 = eye(4)*1+j;               % Matrices can contain complex numbers
Q3.M4 = eye(4)*tf('s');           % or even the Laplace operator

% Q4: LTI Objects
% Create using tf(), zpk() or Laplace Operator tf('s')
% NEVER use chgTimeUnit() - leave in units of (s)
Q4.C  = tf(1, [2 3]);             % Use tf()
Q4.G  = zpk(-4, [-5 -6], 7);      % Use zpk()
Q4.H1 = 2*tf('s')/(3*tf('s')+4);  % Use s operator

% You may find it clearer to pre-define the Laplace operator
s     = tf('s');
Q4.H2 = 2*s/(3*s+4);              % Same as H1 but easier to read

% Q5: Text Message
% Enclose in single quotes.
Q5.X1 = 'Advice may be useful or pleasant, but never both.';

% Use sprintf() to include line breaks and variables.
% See the help screen on 'Formatting Text' to include special characters
% like the single-quote symbol.
pos = 100;
neg = -pos;
Q5.X2 = sprintf('\nIf you think %d and %d are basically the same number,\nlet''s apply that logic to your bank account and see if you\ndon''t change your mind.\n', pos, neg);
Q5.X2 = sprintf('\nIf you think %d and %d are basically the same number,\ntell that to the person who issues you''re patcheck.\nThey''ll be delighted to hear it.\n', pos, neg);

% Use the display() function to display formatted text.
display(Q5.X1);
display(Q5.X2);

% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
