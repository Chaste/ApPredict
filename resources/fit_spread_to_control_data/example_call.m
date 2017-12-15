close all
clear all

load('AZ_hERG_cisapride.mat');

pic50s = AZ_cisapride(:,1);
hills = AZ_cisapride(:,2);

[pic50_spread] = fit_spread_parameters(pic50s)

[pic50_spread, hill_spread] = fit_spread_parameters(pic50s,hills)

