% File:   itsave.m
% Brief:  Saves Matlab/Octave workspace variables to an IT++ itfile
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
%
% Usage: itsave("fname.it", x, y, z) saves variables x, y, and z to the 
%        IT++ file "fname.it"

function itsave(fname, varargin)

if (nargin > 1)
  vars = varargin;
else
  % vars = evalin('caller', 'who -variables')';
  error('syntax error: itsave(fname, variables ...)');
end

nargs = size(vars,2);

endianity = 1; % Always big endian for now...
file_version = 1; % The current file-version of it_file

if (endianity == 0)
  [fid, err_msg] = fopen(fname, 'wb', 'ieee-le');
else
  [fid, err_msg] = fopen(fname, 'wb', 'ieee-be');
end

if (fid == -1)
  error(err_msg);
end

% Write a file header consisting of "IT++" and a char containing the file
% version number 
fprintf(fid, 'IT++%c', file_version);

for ai=1:nargs
  if (exist('OCTAVE_VERSION')) % check for octave
      vname = deblank(argn(ai+1,:)); % octave way of getting parameter name
  else
      vname = inputname(ai+1); % matlab way of getting parameter name
  end
  v = vars{ai};

  if (isa(v, 'double')) % double precision floating point type
    is_scalar = (size(v,1)==1 & size(v,2)==1); % true if scalar (for future use)
    is_vector = (size(v,1)==1 | size(v,2)==1); % true if a vector (or a scalar)

    % Writes a char telling what binary format is used
    fwrite(fid, [endianity], 'char');

    if (isreal(v)) % Check if real values
      hdr_bytes = 1 + 3 * 4 + size(vname,2)+1 + 5;
      if (is_vector)
        data_bytes = 1 * 4 + 8 * prod(size(v));
      else % a matrix
        data_bytes = 2 * 4 + 8 * prod(size(v));
      end
      block_bytes = hdr_bytes + data_bytes;

      % Writes a header
      fwrite(fid, [hdr_bytes data_bytes block_bytes], 'uint32');
      % Writes variable name as string
      fprintf(fid, '%s%c', vname); fwrite(fid, 0, 'char');

      if (is_vector)
        fprintf(fid, 'dvec'); fwrite(fid, 0, 'char');
        fwrite(fid, size(v(:),1), 'uint32');

      else % a matrix
        fprintf(fid, 'dmat'); fwrite(fid, 0, 'char');
        fwrite(fid, size(v), 'uint32');
      end

      fwrite(fid, v, 'float64');

    else % complex values
      % 3*4 is for sized, 6 is for type (dcvec)
      hdr_bytes = 1 + 3 * 4 + size(vname,2)+1 + 6;  
      if (is_vector)
        data_bytes = 1 * 4 + 2 * 8 * prod(size(v));
      else % a matrix
        data_bytes = 2 * 4 + 2 * 8 * prod(size(v));
      end
      block_bytes = hdr_bytes + data_bytes;

      % Writes a header
      fwrite(fid, [hdr_bytes data_bytes block_bytes], 'uint32');
      % Writes variable name as string
      fprintf(fid, '%s', vname);  fwrite(fid, 0, 'char');

      if (is_vector)
        fprintf(fid, 'dcvec'); fwrite(fid, 0, 'char');
        fwrite(fid, size(v(:),1), 'uint32');
        for i=1:size(v(:),1)
          fwrite(fid, real(v(i)), 'float64');
          fwrite(fid, imag(v(i)), 'float64');
        end

      else % a matrix
        fprintf(fid, 'dcmat'); fwrite(fid, 0, 'char');
        fwrite(fid, size(v), 'uint32');
        for j=1:size(v,2)
          for i=1:size(v,1)
            fwrite(fid, real(v(i,j)), 'float64');
            fwrite(fid, imag(v(i,j)), 'float64');
          end
        end
      end

    end % real or complex

  else
    warning(['Variable ''' vname ''' is not a vector or matrix.  Not saved.']);
  end

end

fclose(fid);
