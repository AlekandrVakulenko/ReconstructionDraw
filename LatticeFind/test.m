
clc

Points(1).X = 1.5+0.1;
Points(1).Y = 1.5;
Points(2).X = -1.5-0.1;
Points(2).Y = 1.5;


[LatPar, Angle] = latticefind(Points, true)





















