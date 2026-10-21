function tf = contains(str, pat, varargin)
% Minimal Octave shim for MATLAB's contains (local testing only).
  if iscell(str)
    tf = cellfun(@(s) contains(s, pat, varargin{:}), str);
    return
  end
  if iscell(pat)
    tf = any(cellfun(@(p) ~isempty(strfind(str, p)), pat));
  else
    tf = ~isempty(strfind(str, pat));
  end
end
