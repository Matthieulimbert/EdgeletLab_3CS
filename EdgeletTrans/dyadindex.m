function ix = dyadindex(s,b,n)
% dyadindex -- returns the vector indices for a dyadic interval
%
%  Usage
%    ix = dyadindex(s,b,n)
%  Inputs
%    s     depth of splitting
%    b     dyadic interval index among 2^d possibilities at depth s, 
%    n     length of signal, has to be a dyadic number n=2^J
%  Outputs
%    ix     linear indices of all entries in that dyadic interval

nblocks = 2^s;
ix =  ( (b * (n/nblocks) + 1) : ((b+1)*n/nblocks) );

%
% Note: Code was adapted from WaveLab Version 802
%
% $RCSfile: dyadindex.m,v $
% $Date: 2006/05/01 17:44:15 $
% $Revision: 1.3 $
%
% Copyright (c) Hannes Helgason, California Institute of Technology, 2006
