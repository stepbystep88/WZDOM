function [rangeInline, rangeCrossline] ...
    = bsGetWorkAreaRangeByTimeline(timeLine)
%% read the range of inline and crossline of a segy file
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
    rangeInline = [min(timeLine{1}(:, 1)), max(timeLine{1}(:, 1))];
    rangeCrossline = [min(timeLine{1}(:, 2)), max(timeLine{1}(:, 2))];
end