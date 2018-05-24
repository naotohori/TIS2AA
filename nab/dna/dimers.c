#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nabcode.h"
extern char NAB_rsbuf[];
static int mytaskid, numtasks;

static MOLECULE_T *m;

int main( argc, argv )
	int	argc;
	char	*argv[];
{
	nabout = stdout; /*default*/

	mytaskid=0; numtasks=1;
static STRING_T *__st0001__ = NULL;
static STRING_T *__st0002__ = NULL;
static STRING_T *__st0003__ = NULL;
m = fd_helix( STEMP( __st0001__, "abdna" ), STEMP( __st0002__, "aa" ), STEMP( __st0003__, "dna" ) );
putpdb( "aa_tt.pdb", m, NULL );

m = fd_helix( STEMP( __st0001__, "abdna" ), STEMP( __st0002__, "at" ), STEMP( __st0003__, "dna" ) );
putpdb( "at_at.pdb", m, NULL );

m = fd_helix( STEMP( __st0001__, "abdna" ), STEMP( __st0002__, "ac" ), STEMP( __st0003__, "dna" ) );
putpdb( "ac_gt.pdb", m, NULL );

m = fd_helix( STEMP( __st0001__, "abdna" ), STEMP( __st0002__, "ag" ), STEMP( __st0003__, "dna" ) );
putpdb( "ag_ct.pdb", m, NULL );

m = fd_helix( STEMP( __st0001__, "abdna" ), STEMP( __st0002__, "gg" ), STEMP( __st0003__, "dna" ) );
putpdb( "gg_cc.pdb", m, NULL );

m = fd_helix( STEMP( __st0001__, "abdna" ), STEMP( __st0002__, "ga" ), STEMP( __st0003__, "dna" ) );
putpdb( "ga_tc.pdb", m, NULL );

m = fd_helix( STEMP( __st0001__, "abdna" ), STEMP( __st0002__, "gc" ), STEMP( __st0003__, "dna" ) );
putpdb( "gc_gc.pdb", m, NULL );

m = fd_helix( STEMP( __st0001__, "abdna" ), STEMP( __st0002__, "ta" ), STEMP( __st0003__, "dna" ) );
putpdb( "ta_ta.pdb", m, NULL );

m = fd_helix( STEMP( __st0001__, "abdna" ), STEMP( __st0002__, "tg" ), STEMP( __st0003__, "dna" ) );
putpdb( "tg_ca.pdb", m, NULL );

m = fd_helix( STEMP( __st0001__, "abdna" ), STEMP( __st0002__, "cg" ), STEMP( __st0003__, "dna" ) );
putpdb( "cg_cg.pdb", m, NULL );


	exit( 0 );
}
