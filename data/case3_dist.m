function mpc = case3_dist
%CASE3_DIST  Power flow data for 3 bus radial distribution system
%   modified system of case4_dist, by removing bus 4
%   Please see CASEFORMAT for details on the case file format.

%   MATPOWER

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 1;

%% bus data
% bus_i  type  Pd  Qd  Gs  Bs  area  Vm  Va  baseKV  zone  Vmax  Vmin
mpc.bus = [
  1  1  0.2000  0.2000  0.0000  0.0000  1  1  0  12.5  1  1.05  0.95;
  2  1  0.2000  0.2000  0.0000  0.0000  1  1  0  12.5  1  1.05  0.95;
  3  3  0       0  0.0000  0.0000  1  1  0  12.5  1  1.05  0.95;
];

%% generator data
% bus  Pg  Qg  Qmax  Qmin Vg  mBase  status  Pmax  Pmin  Pc1  Pc2  Qc1min  Qc1max  Qc2min  Qc2max  ramp_agc  ramp_10  ramp_30  ramp_q  apf
mpc.gen = [
3  0.0000  0.0000  0.8  -0.8  1.0500  1  1   0.8  -0.8  0  0  0  0  0  0  0  0  0  0  0;
];

%% branch data
% fbus  tbus  r  x  b  rateA  rateB  rateC  ratio  angle  status  angmin  angmax
mpc.branch = [
	3	1	0.3922	0.2470	0	0.8	0	0	0	0	1	-360	360;
	1	2	0.4930	0.2511	0	0.4	0	0	0	0	1	-360	360;    
];

%% generator cost data
% not relevant but there for parsing
mpc.gencost = [
  2 0 0 2 0 0;
];
