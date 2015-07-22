function c = yr1(a,b,time_pts)

c = 0;

w1 = 0.50;
w2 = 0.25;
w3 = 0.25;

c = (w1 * functionREstrela(a,b)) + (w2 * functionA(a,b,time_pts)) + (w3 * functionM(a,b));

end

%Função M do artigo
function r = functionM(a,b)

[val max_a] = max(a);
[val max_b] = max(b);

[val min_a] = min(a);
[val min_b] = min(b);

if ( (max_a == max_b) && (min_a == min_b) )
    r = 1;
    return
end

if ( (max_a == max_b) || (min_a == min_b) )
    r = 0.5;
    return
end

if ( (max_a ~= max_b) && (min_a ~= min_b) )
    r = 0;
    return
end
end

function r = functionA(a,b,time_pts)

count = 0;

for i=1:(length(a)-1)
    if (functionL(a(i),a(i+1),time_pts(i),time_pts(i+1)) == functionL(b(i),b(i+1),time_pts(i),time_pts(i+1)))
        count = count + 1;
    end
end

r = count / (length(a) - 1);
end

function r = functionREstrela(a,b)

r = (corr(a,b) + 1) / 2;
end

function s = slope(x1,x2,t1,t2)

s = ((x2 - x1) / (t2 - t1));

end

function l = functionL(x1,x2,t1,t2)

l = sign(slope(x1,x2,t1,t2));

end