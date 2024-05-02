% Get the block name as a user input.
BlockName = input('Please enter the blockname to plot: ','s')
load(BlockName);

% Get the number of data points for neural data and TTL data.
num_datapoint1 = max(size(RawData.streams.Audt.data));
num_datapoint2 = max(size(RawData.streams.TTLv.data));

% Define time sequences to convert the unit of x axis.
neural_fs = RawData.streams.NRaw.fs;
ttl_fs = RawData.streams.TTLv.fs;
t1 = (1/neural_fs)*linspace(0,num_datapoint1-1,num_datapoint1);
t2 = (1/ttl_fs)*linspace(0,num_datapoint2-1,num_datapoint2);

% Plot Audt, NRaw, and TTLv from the RawData.
% figure
% plot(t1,RawData.streams.Audt.data); title "Audt"
% hold on
% figure
% plot(t1,RawData.streams.NRaw.data); title "Neural"
% figure
% plot(t2,RawData.streams.TTLv.data(1,:)); title "TTL"

% Plot only Audt and TTLv.
figure
plotyy(t1,RawData.streams.Audt.data,t2,RawData.streams.TTLv.data(1,:))