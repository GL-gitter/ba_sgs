function result = superbee_limiter(a, b)
    if a * b <= 0
        result = 0;
    else
        s = sign(a);
        c1 = min(2*abs(a), abs(b));
        c2 = min(abs(a), 2*abs(b));
        result = s * max(c1, c2);
    end
end
