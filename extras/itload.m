% File:   itload.m
% Brief:  Load an IT++ itfile content to Matlab/Octave workspace
% Author: Tony Ottosson
%
% $Date$
% $Revision$
%
% -------------------------------------------------------------------------
%
% IT++ - C++ library of mathematical, signal processing, speech processing,
%        and communications classes and functions
%
% Copyright (C) 1995-2005  (see AUTHORS file for a list of contributors)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
%
% -------------------------------------------------------------------------

function itload(fname)

endian = -1;

[fid, err_msg] = fopen(fname, 'rb');
if (fid == -1)
  fname = [fname '.it'];
  [fid, err_msg2] = fopen(fname, 'rb');
  if (fid == -1)
    error(err_msg);
  end
end

[d, n] = fread(fid, 5, 'char');
if (n ~= 5 | d(1:4) ~= [73 84 43 43]')
  error('Not an IT++ file!');
end

if (d(5) > 2)
  error('Unknown IT++ file version');
end

while (1)
  p = ftell(fid);
  [d, n] = fread(fid, 1, 'char');
  if (n ~= 1)
    break;
  end
  if (d ~= endian)
    fclose(fid);
    if (d == 0)
      [fid, err_msg] = fopen(fname, 'rb', 'ieee-le');
    else
      [fid, err_msg] = fopen(fname, 'rb', 'ieee-be');
    end
    endian = d;
  end
  fseek(fid, p+1, 'bof');
  [d1, n] = fread(fid, 3, 'uint32'); % Read header, data, and total block sizes
  name = fgetstr(fid);
  type = fgetstr(fid);

  fseek(fid, p+d1(1), 'bof');
  if (length(type) == 0)
% A deleted entry -> skip it

  % ------- bin -----------------------------
  elseif (strcmp(type, 'bin'))
    [d, n] = fread(fid, 1, 'char');
    assignin('caller', name, d);
  % ------- int16 (=short) -----------------------------
  elseif (strcmp(type, 'int16'))
    [d, n] = fread(fid, 1, 'int16');
    assignin('caller', name, d);
  % ------- int32 (=int) -----------------------------
  elseif (strcmp(type, 'int32'))
    [d, n] = fread(fid, 1, 'int32');
    assignin('caller', name, d);
  % ------- int64 (=long_long) -----------------------------
  elseif (strcmp(type, 'int64'))
    [d, n] = fread(fid, 1, 'int64');
    assignin('caller', name, d);
  % ------- float32 (=float) -----------------------------
  elseif (strcmp(type, 'float32'))
    [d, n] = fread(fid, 1, 'float32');
    assignin('caller', name, d);
  % ------- float64 (=double) -----------------------------
  elseif (strcmp(type, 'float64'))
    [d, n] = fread(fid, 1, 'float64');
    assignin('caller', name, d);
  % ------- float32_complex -----------------------------
  elseif (strcmp(type, 'float32_complex'))
    [d, n] = fread(fid, 2, 'float32');
    d = complex(d(1), d(2));
    assignin('caller', name, d);
  % ------- float64_complex -----------------------------
  elseif (strcmp(type, 'float64_complex'))
    [d, n] = fread(fid, 2, 'float64');
    d = complex(d(1), d(2));
    assignin('caller', name, d);

  % ------- bvec -----------------------------
  elseif (strcmp(type, 'bvec'))
    [size, n] = fread(fid, 1, 'uint32');
    [d, n] = fread(fid, size, 'char');
    assignin('caller', name, d);
  % ------- ivec -----------------------------
  elseif (strcmp(type, 'ivec'))
    [size, n] = fread(fid, 1, 'uint32');
    [d, n] = fread(fid, size, 'int32');
    assignin('caller', name, d);
  % ------- llvec -----------------------------
  elseif (strcmp(type, 'llvec'))
    [size, n] = fread(fid, 1, 'uint32');
    [d, n] = fread(fid, size, 'int64');
    assignin('caller', name, d);
  % ------- fvec -----------------------------
  elseif (strcmp(type, 'fvec'))
    [size, n] = fread(fid, 1, 'uint32');
    [d, n] = fread(fid, size, 'float32');
    assignin('caller', name, d);
  % ------- dvec -----------------------------
  elseif (strcmp(type, 'dvec'))
    [size, n] = fread(fid, 1, 'uint32');
    [d, n] = fread(fid, size, 'float64');
    assignin('caller', name, d);
  % ------- fcvec -----------------------------
  elseif (strcmp(type, 'fcvec'))
    [size, n] = fread(fid, 1, 'uint32');
    [d, n] = fread(fid, size*2, 'float32');
    d = complex(d(1:2:end), d(2:2:end));
    assignin('caller', name, d);
  % ------- dcvec -----------------------------
  elseif (strcmp(type, 'dcvec'))
    [size, n] = fread(fid, 1, 'uint32');
    [d, n] = fread(fid, size*2, 'float64');
    d = complex(d(1:2:end), d(2:2:end));
    assignin('caller', name, d);
  % ------- string -----------------------------
  elseif (strcmp(type, 'string'))
    [size, n] = fread(fid, 1, 'uint32');
    [d, n] = fread(fid, size, 'char');
    d = char(d);
    assignin('caller', name, d);
  % ------- bmat -----------------------------
  elseif (strcmp(type, 'bmat'))
    [size, n] = fread(fid, 2, 'uint32');
    [d, n] = fread(fid, size', 'char');
    assignin('caller', name, d);
  % ------- imat -----------------------------
  elseif (strcmp(type, 'imat'))
    [size, n] = fread(fid, 2, 'uint32');
    [d, n] = fread(fid, size', 'int32');
    assignin('caller', name, d);
  % ------- llmat -----------------------------
  elseif (strcmp(type, 'llmat'))
    [size, n] = fread(fid, 2, 'uint32');
    [d, n] = fread(fid, size', 'int64');
    assignin('caller', name, d);
  % ------- fmat -----------------------------
  elseif (strcmp(type, 'fmat'))
    [size, n] = fread(fid, 2, 'uint32');
    [d, n] = fread(fid, size', 'float32');
    assignin('caller', name, d);
  % ------- dmat -----------------------------
  elseif (strcmp(type, 'dmat'))
    [size, n] = fread(fid, 2, 'uint32');
    [d, n] = fread(fid, size', 'float64');
    assignin('caller', name, d);
  % ------- fcmat -----------------------------
  elseif (strcmp(type, 'fcmat'))
    [size, n] = fread(fid, 2, 'uint32');
    [d, n] = fread(fid, size(1)*size(2)*2, 'float32');
    d = reshape(complex(d(1:2:end), d(2:2:end)), size(1), size(2));
    assignin('caller', name, d);
  % ------- dcmat -----------------------------
  elseif (strcmp(type, 'dcmat'))
    [size, n] = fread(fid, 2, 'uint32');
    [d, n] = fread(fid, size(1)*size(2)*2, 'float64');
    d = reshape(complex(d(1:2:end), d(2:2:end)), size(1), size(2));
    assignin('caller', name, d);
  % ------- bArray -----------------------------
  elseif (strcmp(type, 'bArray'))
    [size, n] = fread(fid, 1, 'uint32');
    [d, n] = fread(fid, size, 'char');
    assignin('caller', name, d);
  % ------- iArray -----------------------------
  elseif (strcmp(type, 'iArray'))
    [size, n] = fread(fid, 1, 'uint32');
    [d, n] = fread(fid, size, 'int32');
    assignin('caller', name, d);
  % ------- llArray -----------------------------
  elseif (strcmp(type, 'llArray'))
    [size, n] = fread(fid, 1, 'uint32');
    [d, n] = fread(fid, size, 'int64');
    assignin('caller', name, d);
  % ------- fArray -----------------------------
  elseif (strcmp(type, 'fArray'))
    [size, n] = fread(fid, 1, 'uint32');
    [d, n] = fread(fid, size, 'float32');
    assignin('caller', name, d);
  % ------- dArray -----------------------------
  elseif (strcmp(type, 'dArray'))
    [size, n] = fread(fid, 1, 'uint32');
    [d, n] = fread(fid, size, 'float64');
    assignin('caller', name, d);
  % ------- fcArray -----------------------------
  elseif (strcmp(type, 'fcArray'))
    [size, n] = fread(fid, 1, 'uint32');
    [d, n] = fread(fid, size*2, 'float32');
    d = complex(d(1:2:end), d(2:2:end));
    assignin('caller', name, d);
  % ------- dcArray -----------------------------
  elseif (strcmp(type, 'dcArray'))
    [size, n] = fread(fid, 1, 'uint32');
    [d, n] = fread(fid, size*2, 'float64');
    d = complex(d(1:2:end), d(2:2:end));
    assignin('caller', name, d);


  % ------- else -----------------------------
  else
    warning(['Unknown type: ' type]);
  end
  fseek(fid, p+d1(3), 'bof');
end

fclose(fid);

  

function str = fgetstr(fid)
str = '';
while (1)
  [d, n] = fread(fid, 1, 'char');
  if (d == 0)
    break;
  end
  str = [str char(d)];
end
