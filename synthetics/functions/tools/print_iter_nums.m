function print_iter_nums(ii,nn,printInt)
% Prints iteration variable to stdout if it is a multiple of printInt.
% Starts new line every 10^th print.
% menandrin@gmail.com, 200110

% EXAMPLE
% nn = 5e4;
% for ii=1:nn
%     print_iter_nums(ii,nn,1e3)
% end

printValues   = linspace(printInt,1e3*printInt,1e3*printInt/printInt);
isNewLine     = rem(printValues,printInt*10)==0;
newLineValues = printValues(isNewLine);
newLineValues = [1, newLineValues, nn];
printValues   = [1, printValues  , nn];

if ismember(ii,newLineValues); fprintf(1,sprintf('\n%i -- ',nn)); end
if ismember(ii,printValues);   fprintf(1,sprintf('%i .. ',ii)); end
if ii==nn; fprintf(1,' done.\n'); end