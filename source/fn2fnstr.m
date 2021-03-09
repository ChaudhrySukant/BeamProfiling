function fns=fn2fnstr(fn);
%fns=fn2fnstr(fn) returns the filename as a TeX interpretable string
%i.e. replaces \ with \\, _ with \_

% created      : A. Wieser
% last modified: AW 2/23/2004


rpl = {'\','\\';...
       '_','\_'};          % replace assignments
       
fns = fn;
for i=1:size(rpl,1)
   fns = strrep(fns,rpl{i,1},rpl{i,2});
end
