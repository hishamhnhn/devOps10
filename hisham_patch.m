%Program to calculate the parameters to design a rectangular patch antenna
%the user have to feed the values of frequency, dielectric constant, and
%height of the dielectric.
%program will calculate automatically the width and length of the patch
%and the width and length of the transition and transmission feed line.

f = input( 'input frequency f in Ghz: ');
Er = input ( 'input dielectric constant of the substrate ');
h = input( 'input height of substrate h in mm: ');
h=h/1000; %(in mm)
f=f*1e9; % turn frequency to HZ
c = 3e8; % speed of light
% calculating Width and Length of the Patch (Read Balis book Chapeter 14)
W = ( c / ( 2 * f ) ) * ( ( 2 / ( Er + 1 ) )^0.5); % Width of the patch
Er_eff = (Er+1)/2 + (( Er -1 )/2)*(1/(sqrt(1+(12*(h/W)))));  % Effective permitivity
L_eff = c/(2*f*sqrt(Er_eff)) % Effective length
a1 = ( Er_eff + 0.3 ) * ( ( W / h ) + 0.264 );
a2 = ( Er_eff - 0.258 ) * ( ( W / h ) + 0.8 );
delta_L = (0.412 * ( a1 / a2 )) * h;
L = L_eff - 2*delta_L; % Length of the patch
str=['width of the patch = ', num2str(W*1000), ' mm']
str=['length of the patch = ', num2str(L*1000), ' mm']
% Calculating the input impedance of the patch
Zo = 90 * Er^2*(L/W)^2/(Er-1);
Zl=90*((Er^2)/(Er-1))*(L/W)^2

% Calculating the strip transition line
Zt=sqrt(50*Zo);
% Zt=50;
a3=exp(Zt*sqrt(Er)/60); p=-4*h*a3; q=32*h^2;
Wt1=-(p/2) + sqrt((p/2)^2-q);
Wt2=-(p/2) - sqrt((p/2)^2-q); %width of the transition line
Er_t= (Er+1)/2 + (( Er -1 )/2)*(1/(sqrt(1+(12*(h/Wt2)))));
L_t=(c/f)/(4*sqrt(Er_t)); %length of transition line
str=['width of the transition line = ', num2str(Wt2*1000), ' mm']
str=['length of transition line = ', num2str(L_t*1000), ' mm']
str=['The impedance of the load is = ', num2str(Zo), ' Ohm']
str=['The impedasnce of transition line = ',num2str(Zt), ' Ohm']