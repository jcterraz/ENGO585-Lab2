function [A] = compute_A(est_coords, range, target)

A = zeros(4,2);
A(1,1) = (target.one(1) - est_coords(1)) / range.one;

