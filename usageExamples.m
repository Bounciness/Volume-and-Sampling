
% Toy examples of sampling and volume of a cube
P = makeBody('cube',10);
X = genSamples(P,100,1000);
V_cube = Volume(P);

%a more complicated example, needs to be preprocessed/rounded
%P_complicated a predefined polytope, with fields
% A
% b A*x <= b
% A_eq
% b_eq A_eq*x = b_eq

%preprocess the polytope to get it in an acceptable format for the
%sampling/volume routines
[rP] = preprocess(P_complicated);

%now for samples or volume...
X = genSamples(rP,100,1000);
V = Volume(rP);