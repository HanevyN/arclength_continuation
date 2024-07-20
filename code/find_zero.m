function f_0 = find_zero(x1,x2, y1,y2)
% linear interpolation of a zero between two points

    ratio = (0-y1) / (y2-y1); % interpolate between to find 0
    f_0 = x1 + ratio*(x2-x1);

end