%chk=+ch4.chk
%nprocshared=4
#hf/cc-pvqz freq=(HPmodes,saveNM) geom=connectivity pop=full

Title

1 2
 C         -0.0000000004        0.0000000072       -0.0000000001
 H          0.7569725982        0.7569725837        0.3775945777
 H         -0.7569725797       -0.7569725942        0.3775945726
 H          0.7569725889       -0.7569725865       -0.3775945748
 H         -0.7569725891        0.7569725915       -0.3775945753

 1 2 1.0 3 1.0 4 1.0 5 1.0
 2 1 1.0
 3 1 1.0
 4 1 1.0
 5 1 1.0

