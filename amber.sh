tleap -f leap.in
sander -O -i min.in -o mini.out -p prmtop -c inpcrd -ref inpcrd -r mini.rst
ambpdb -p prmtop -c mini.rst > f024.aa.min.pdb
