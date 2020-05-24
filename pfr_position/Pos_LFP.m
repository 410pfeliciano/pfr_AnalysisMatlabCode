position = [Pos1; Pos2; Pos3];

%%
[LFP, LFPt] = load_open_ephys_data_faster('100_CH8.continuous');

LFPt(end)- positiondata.time(end)

%%
sites = MazeDivision(positiondata.position, 5);