# run as octave
# graddescent
u=-1:0.01:1;
v=-1:0.01:1; 
[X,Y]=meshgrid(u,v);
x=X;
y=Y;
z=4.*X.^2 - 4.*X.*Y + 2.*Y.^2;
r=surf(x,y,z);


