function index=find_closest_member(A,x)

tmp = abs(A-x);
[vx, index] = min(tmp);

end