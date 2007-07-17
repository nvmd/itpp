% File:   itsave.m
% Brief:  Saves Matlab/Octave workspace variables to an IT++ itfile
% Author: Tony Ottosson and Adam Piatyszek
%
% Usage: itsave("fname.it", var1, [var2], ...)
%
% This function saves a set of Matlab/Octave workspace variables to an IT++
% file format. Currently, only vectors and 2-D matrices can be saved. The
% type of data saved is detected automatically and can be one of the
% following types: bvec, bmat, ivec, imat, vec, mat, cvec, cmat.
%
% -------------------------------------------------------------------------
%
% IT++ - C++ library of mathematical, signal processing, speech processing,
%        and communications classes and functions
%
% Copyright (C) 1995-2007  (see AUTHORS file for a list of contributors)
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

function itsave(fname, varargin)

if (nargin > 1)
  vars = varargin;
else
  error('Usage: itsave("fname.it", var1, [var2], ...)');
end

nargs = size(vars,2);

% The current file-version of it_file
file_version = 3;

[fid, err_msg] = fopen(fname, 'wb', 'ieee-le');
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

  is_scalar = all(size(v) == 1); % true if scalar (for future use)
  is_vector = (sum(size(v) > 1) <= 1); % true if a vector (or a scalar)

  if (~imag(v) && (v == round(v))) % binary or integer type
    if ((max(v) <= 1) && (min(v) >= 0)) % binary type
      % Calculate sizes
      if (is_vector)
	hdr_bytes = 3 * sizeof(uint64(0)) + size(vname,2)+1 ...
	    + sizeof('bvec')+1 +1;
	data_bytes = sizeof(uint64(0)) + sizeof(char(0)) * prod(size(v));
      else % a matrix
	hdr_bytes = 3 * sizeof(uint64(0)) + size(vname,2)+1 ...
	    + sizeof('bmat')+1 +1;
	data_bytes = 2 * sizeof(uint64(0)) + sizeof(char(0)) * prod(size(v));
      end
      block_bytes = hdr_bytes + data_bytes;

      % Write header sizes
      fwrite(fid, [hdr_bytes data_bytes block_bytes], 'uint64');
      % Write variable name as string
      fprintf(fid, '%s%c', vname); fwrite(fid, 0, 'char');

      % Write data type string, empty description string and data size
      if (is_vector)
	fprintf(fid, 'bvec'); fwrite(fid, 0, 'char'); fwrite(fid, 0, 'char');
	fwrite(fid, size(v(:),1), 'uint64');
      else % a matrix
	fprintf(fid, 'bmat'); fwrite(fid, 0, 'char'); fwrite(fid, 0, 'char');
	fwrite(fid, size(v), 'uint64');
      end
      % Write data
      fwrite(fid, v, 'char');

    else % integer type
      % Calculate sizes
      if (is_vector)
	hdr_bytes = 3 * sizeof(uint64(0)) + size(vname,2)+1 ...
	    + sizeof('ivec')+1 +1;
	data_bytes = sizeof(uint64(0)) + sizeof(int32(0)) * prod(size(v));
      else % a matrix
	hdr_bytes = 3 * sizeof(uint64(0)) + size(vname,2)+1 ...
	    + sizeof('imat')+1 +1;
	data_bytes = 2 * sizeof(uint64(0)) + sizeof(int32(0)) * prod(size(v));
      end
      block_bytes = hdr_bytes + data_bytes;

      % Write header sizes
      fwrite(fid, [hdr_bytes data_bytes block_bytes], 'uint64');
      % Write variable name as string
      fprintf(fid, '%s%c', vname); fwrite(fid, 0, 'char');

      % Write data type string, empty description string and data size
      if (is_vector)
	fprintf(fid, 'ivec'); fwrite(fid, 0, 'char'); fwrite(fid, 0, 'char');
	fwrite(fid, size(v(:),1), 'uint64');
      else % a matrix
	fprintf(fid, 'imat'); fwrite(fid, 0, 'char'); fwrite(fid, 0, 'char');
	fwrite(fid, size(v), 'uint64');
      end
      % Write data
      fwrite(fid, v, 'int32');
    end % binary or integer

  elseif (isa(v, 'double')) % double precision floating point type
    if (isreal(v)) % Check if real values
      % Calculate sizes
      if (is_vector)
	hdr_bytes = 3 * sizeof(uint64(0)) + size(vname,2)+1 ...
	    + sizeof('dvec')+1 + 1;
	data_bytes = sizeof(uint64(0)) + sizeof(double(0)) * prod(size(v));
      else % a matrix
	hdr_bytes = 3 * sizeof(uint64(0)) + size(vname,2)+1 ...
	    + sizeof('dmat')+1 + 1;
	data_bytes = 2 * sizeof(uint64(0)) + sizeof(double(0)) ...
	    * prod(size(v));
      end
      block_bytes = hdr_bytes + data_bytes;

      % Write a header
      fwrite(fid, [hdr_bytes data_bytes block_bytes], 'uint64');
      % Writes variable name as string
      fprintf(fid, '%s%c', vname); fwrite(fid, 0, 'char');

      if (is_vector)
	fprintf(fid, 'dvec'); fwrite(fid, 0, 'char'); fwrite(fid, 0, 'char');
	fwrite(fid, size(v(:),1), 'uint64');
      else % a matrix
	fprintf(fid, 'dmat'); fwrite(fid, 0, 'char'); fwrite(fid, 0, 'char');
	fwrite(fid, size(v), 'uint64');
      end
      fwrite(fid, v, 'float64');

    else % complex values
      % Calculate sizes
      if (is_vector)
	hdr_bytes = 3 * sizeof(uint64(0)) + size(vname,2)+1 ...
	    + sizeof('dcvec')+1 + 1;
	data_bytes = sizeof(uint64(0)) + 2 * sizeof(double(0)) ...
	    * prod(size(v));
      else % a matrix
	hdr_bytes = 3 * sizeof(uint64(0)) + size(vname,2)+1 ...
	    + sizeof('dcmat')+1 + 1;
	data_bytes = 2 * sizeof(uint64(0)) + 2 * sizeof(double(0)) ...
	    * prod(size(v));
      end
      block_bytes = hdr_bytes + data_bytes;

      % Writes header sizes
      fwrite(fid, [hdr_bytes data_bytes block_bytes], 'uint64');
      % Write variable name as string
      fprintf(fid, '%s', vname);  fwrite(fid, 0, 'char');

      if (is_vector)
	% Write data type string, empty description string and data size
	fprintf(fid, 'dcvec'); fwrite(fid, 0, 'char'); fwrite(fid, 0, 'char');
	fwrite(fid, size(v(:),1), 'uint64');
	% Write data
	for i=1:size(v(:),1)
	  fwrite(fid, real(v(i)), 'float64');
	  fwrite(fid, imag(v(i)), 'float64');
	end
      else % a matrix
        % Write data type string, empty description string and data size
	fprintf(fid, 'dcmat'); fwrite(fid, 0, 'char'); fwrite(fid, 0, 'char');
	fwrite(fid, size(v), 'uint64');
	% Write data
	for j=1:size(v,2)
	  for i=1:size(v,1)
	    fwrite(fid, real(v(i,j)), 'float64');
	    fwrite(fid, imag(v(i,j)), 'float64');
	  end
	end
      end
    end % real or complex
  else
      warning(['Variable ''' vname ''' is neither a vector nor matrix. Not saved.']);
  end

end

fclose(fid);
