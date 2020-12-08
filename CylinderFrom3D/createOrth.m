


% creates two orthogonal vectors to present

function [orth1 orth2] = createOrth(vect)
orth1 = [0 0 0];

[~, index] = max(abs(vect));


coord1 = mod(index,3)+1;
coord2 = mod(index+1,3)+1;

orth1(coord1) = 1;
orth1(coord2) = 0;

orth1(index) = -vect(coord1)./vect(index);
orth2 = cross(vect, orth1);

orth1 = orth1./norm(orth1);
orth2 = orth2./norm(orth2);

end




















