function mlist = add_blinking_fields(coords)
    mlist = struct();
    mlist.master.x = coords(:, 1);
    mlist.master.xc = mlist.master.x;
    mlist.blinking.newx = mlist.master.xc;
    mlist.master.y = coords(:, 2);
    mlist.master.yc = mlist.master.y;
    mlist.blinking.newy = mlist.master.yc;
    mlist.blinking.Lifetime = ones(size(mlist.master.x));
end
