function clustInfo = importfile_clustinfo(infofile, dataLines)
%IMPORTFILE Import data from a text file
%  CLUSTERINFONEW = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  CLUSTERINFONEW = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  clustInfo = importfile("/home/jaejin/Parooa/Synapse/i/kiloSorted_DMR/MrCassius-190326/D2_AC_R1/KS2_7_AC/ClusterInfo/cluster_info_new.tsv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 24-Aug-2021 12:58:05

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 11);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["id", "Amplitude", "ContamPct", "KSLabel", "amp", "ch", "depth", "fr", "group", "n_spikes", "sh"];
opts.VariableTypes = ["double", "double", "double", "categorical", "double", "double", "double", "double", "categorical", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["KSLabel", "group"], "EmptyFieldRule", "auto");

% Import the data
clustInfo = readtable(infofile, opts);

end
