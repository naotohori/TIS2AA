tleap -f leap.in
rm mini.* mdinfo
sander -i min.in -o mini.out -p prmtop -c inpcrd -r mini.rst
