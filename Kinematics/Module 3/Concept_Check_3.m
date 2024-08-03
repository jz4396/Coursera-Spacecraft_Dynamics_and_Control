%% Problem 2
R = eulerAd([3,2,1],[20,-10,120])
[n2v,n2e] = eig(R)
n2r = acos((trace(R)-1)/2)


%% Problem 3
R1=[1 0 0;0 0 1;0 -1 0]
R=R1*R1
[n3v,n3e] = eig(R)
n3r = acos((trace(R)-1)/2)
function R=rotated(dim,deg)
    orders = [1,2,3;3,1,2;2,3,1];
    r = [1 0 0; 0 cosd(deg) sind(deg); 0 -sind(deg) cosd(deg)];
    R=r(orders(dim,:),orders(dim,:));
end
function R=eulerAd(dims,rads)
    R=eye(3);
    for i = 1:length(dims)
        R=rotated(dims(i),rads(i))*R;
    end
end