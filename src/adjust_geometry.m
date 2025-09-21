function Pn = adjust_geometry(P0, DD, DH, TTH)
    Pn = P0(1,:);
    if size(P0,1) == 2
        Pn(1) = P0(1,1) + DD - P0(2,1);
        Pn(3) = P0(1,3) + (P0(2,3) - DH)/P0(2,4);
    end
    Pn(5) = TTH;
end