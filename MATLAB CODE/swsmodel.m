function [F,M] = swsmodel(D,R,P,resample_rate)
% [F,M] = swsmodel(D,R,H)  Sine wave speech analysis
%       D is a speech example sampled at R samples per second.
%       Return a sinusoid model of up to 4 components with each sinusoid
%       defined by a row of F (frequencies in Hz) and M (linear magnitude). 
%       Each column of F and M corresponds to H samples at R.
%       Rows of F are sorted with lowest frequency first; 
%       sinusoids cannot cross.
%
%       Relies on lpcfit.m and lpca2frq.m to form LPC model and convert it 
%       into frequencies.
% 2001-03-12 dpwe@ee.columbia.edu $Header: $

% Amount of samples to group together
H = 128;

% Target sampling rate
MyR = resample_rate;

p = round(MyR/1000);
q = round(R/1000);

% Resample to 8 kHz, so LPC only picks main formants
D = resample(D, p, q);
size(D)

% Step size in units of my sampling rate (force to be even integer)
HH = 2*round(H/R * MyR / 2);

% Form 8th-order LPC model (3 or 4 pole pairs)
[lpca,g] = lpcfit(D,P*2,HH,2*HH);

% Convert poles to sorted freqs and magnitudes
% If only 3 nonzero freqs are found, 4th row will have mag/frq zero
[fa, ma] = lpca2frq(lpca,g);

% Convert frqs into Hz
F = fa'*MyR/(2*pi);

M = ma';

