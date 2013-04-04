% File:         SACluster.m
% Description:  Implementation of SACluster algorithm
% Chen,Zhiwei, 2013.4.1


function SACluster(graph,l,c,sita,k)
w = ones(m,1);
temp = c*(1-c)*PA;
RA = 0;
while(l>0)
    RA = RA + temp;
    temp = temp*(1-c)*PA;
    l = l-1;
end

fb = 1 - exp(D.^2/(2*sita))
end
