function result = minmod_limiter(a, b, c)
    if (a > 0 && b > 0 && c > 0)
        result = min([a, b, c]);
    elseif (a < 0 && b < 0 && c < 0)
        result = max([a, b, c]);
    else
        result = 0;
    end
end