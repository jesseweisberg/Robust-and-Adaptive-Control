function t1 = cross(a,t)
%find the value of t where a crosses zero
n=prod(size(a));
j=0;
i=0;

while j < n,
  if a(n-j) > 0.,
    i=n-j;
    j=n+1;
  end;
  j=j+1;
end

% Index out of bounds check
if (i) < 1
    % Return rise time of 0 for error
    t1 = 0;
    return
end

pp=inv([t(i) 1.;t(i+1) 1.])*[a(i);a(i+1)];
t1=-pp(2)/pp(1);
return