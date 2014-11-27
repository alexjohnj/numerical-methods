% Reads data from a plaintext file containing two columns, x and f(x).
% The script then calculates the value f(x) for a user provided number
% by interpolating the data.

% Read File
fileName = input('Enter filename >> ', 's');
rawData = importdata(fileName);

% Get xint (the number we want to calculate an interpolated value of)
xint = NaN;

while true
    xint = input('Enter the value of x you want to interpolate in GPa >> ');
    if xint > max(rawData(:,1)) || xint < min(rawData(:,1))
       fprintf('%.2f does not lie in the data range\n', xint); 
    else
        break;
    end
end

% Check if f(xint) exists
[containsXint, index] = ismember(xint, rawData(:, 1));
if containsXint
  fprintf('f(%.1f) = %.2fÅ³\n', xint, rawData(index, 2));
  return;
end

% Find bracketing values
xLower = NaN;
xUpper = NaN;
fxLower = NaN;
fxUpper = NaN;

for ii = 2:length(rawData)
  if rawData(ii, 1) >= xint
    xLower = rawData(ii-1, 1);
    xUpper = rawData(ii, 1);

    fxLower = rawData(ii-1, 2);
    fxUpper = rawData(ii, 2);
    break;
  end
end

fxint = fxLower + (((fxUpper - fxLower) / (xUpper - xLower)) * (xint - xLower));
fprintf('f(%.1f) = %.2fÅ³\n', xint, fxint);
