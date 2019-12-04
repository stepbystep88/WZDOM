function possion = bsGetPossion(vp, vs)
    vp_vs = (vp ./ vs).^2;
    possion = (0.5*vp_vs - 1)./(vp_vs - 1);
end