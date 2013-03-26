% Content: Spectral algorithm, created on 2013.3.23

function spectral(adj,k)
g = graph;
set_matrix(g,adj);
distxy(g);
xy = getxy(g);
Lap = diag(deg(g)) - adj;
[V,D] = eig(Lap);   
vec = V(:,1:k);
index = kmeans(vec,k);

for i = 1:k
    cx = 0;
    cy = 0;
    num = 0;
    maxR = 0;
    for j = 1:size(adj,1)
        if index(j) == i
            cx = cx + xy(j,1);
            cy = cy + xy(j,2);
            num = num + 1;
        end
    end
    cx = cx/num;
    cy = cy/num;
    
    for j = 1:size(adj,1)
        if index(j) == i
           r = sqrt((cx - xy(j,1))^2 + (cy - xy(j,2))^2);
            if r > maxR
                maxR = r;
            end
        end
    end
    
    plotcircle(cx,cy,maxR);
    hold on;
end

ndraw(g)
end