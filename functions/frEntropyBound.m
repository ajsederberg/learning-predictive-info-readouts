function [H_per_spk, H] = frEntropyBound(fr)
% entropy per spike in bits given a firing rate (spikes per bin)
H_per_spk = -log2(fr) - (1 - fr).*log2(1 - fr)./fr;
H = H_per_spk.*fr;
