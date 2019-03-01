Usage:

python do_it_all.py [-h] [-b MILLER_BOUNDS] [-l NLAYERS] [-v VACUUM]
                    [-o OUTDIR] [-f FILES]

optional arguments:
  -h, --help            show this help message and exit
  -b MILLER_BOUNDS, --miller_bounds MILLER_BOUNDS
                        Maximum index of the combinations of miller indices to
                        try
  -l NLAYERS, --nlayers NLAYERS
                        Number of layers for the output
  -v VACUUM, --vacuum VACUUM
                        Size of vacuum between slabs in angstroms
  -o OUTDIR, --outdir OUTDIR
                        Output directory
  -f FILES, --files FILES
                        Folder to look for the cif files. ex: ./*/*cif

Usage example:

srun -n 16 python do_it_all.py -b 1 -o "testout" -f "../*/*cif"