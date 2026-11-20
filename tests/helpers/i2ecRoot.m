function root = i2ecRoot()
%I2ECROOT  Absolute path of the I2EC_2d repository root.
%
%   root = i2ecRoot()
%
% Resolved from this file's location (tests/helpers/), so the tests do not
% depend on the current folder.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));

end
