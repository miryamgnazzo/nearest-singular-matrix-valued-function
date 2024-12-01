function error = discret_domain_double(N, d, A, f)
% d = degree of the polynomilal approximation
% N = number of points in the discretization procedure
% f = function to approximate

    pp = 2*pi*rand(1,N);
    
    ss =sqrt(rand(1,N));
    
    zr = cos(pp);
    zi = sin(pp);
    
    z = ss.*(zr + 1i*zi);
    
    points = exp((2*pi*1i).*((1:(d+1))./(d+1)));
    fk = evaluation_det(A, f, points);
    fk = fk.';
    
    V = fliplr(vander(points));
    c = V\fk;
    pz = polyval(flipud(c),z);
    
    fz = evaluation_det(A, f, z);
    error = max(abs(pz-fz));

kv = 2:2:8;
i = 1;

while (i <= length(kv))
    k = kv(i);
    N = N*k;

    pp = 2*pi*rand(1,N);
    
    ss =sqrt(rand(1,N));
    
    zr = cos(pp);
    zi = sin(pp);
    
    z = ss.*(zr + 1i*zi);
    
    points = exp((2*pi*1i).*((1:(d+1))./(d+1)));
    fk = evaluation_det(A, f, points);
    fk = fk.';
    
    V = fliplr(vander(points));
    c = V\fk;
    pz = polyval(flipud(c),z);
    
    fz = evaluation_det(A, f, z);
    error_new = max(abs(pz-fz));

    if (abs(error - error_new) <= 10^-12)
        i = 20;
    end

    error = error_new;

    i = i+1;
end
