/*
 *            header file for  DYOP SIMULATOR
 *
 * -Cem Bozsahin
 *  Arizona State University
 *  1986
 */

/*	parameters		*/

#define MAXATR	16	/* max. # of attributes	 	*/
#define MAXPAGE 4000 	/* max. # of allocatable pages	*/
#define MAXVAL  017777777777  /* 2**31 -1			*/
#define ATRS    10	/* attrib. size for simulator  	*/
#define IMPLICIT -1	/* return value if partition is
			 * not in the partition table
			 */
#define TSTORED 001	/* tuple status bytes 		*/
#define TMOVED  002

/* -----------------------------------------------------*/

/* 	exogeneous variables	*/

/* domain information table	*/
struct d
{
	long range;	/* value range			*/
	int  size;	/* attrib. length in bytes	*/
} d [ MAXATR ];
int nattrib;		/* # of attrib. in the relation */
int pagesize;		/* page size for a partition(byte)*/
int b0;			/* page size in # of tuples	*/
int tsize;		/* tuple size (bytes)		*/
int stsize;		/* simulator tuple size : all attributes
			 * are hashed into range 2**31-1*/
int spsize;		/* simulator page size : b0 * stsize	*/
int seed;
int maxtupl;

/* -----------------------------------------------------*/

/* endogeneous variables	*/

/* partition table : serves as 1 st level directory	*/
struct p
{
	long pno;	/* partiton #			*/
	int ntup;	/* # of tuples in the part.	*/
	long paddr;	/* base address of part. on disk*/
	char splx;	/* split axis the partition has 
			 * been generated
			 */
	long plow;	/* lower bound of partition	*/
	long phigh;	/* upper bound of partition 	*/
} p [ MAXPAGE ];
int npart;		/* # of explicit partitions     */
int nsplit [ MAXATR ];  /* # of splits in domain i	*/
int tsplit;		/* total # of splits in whole data space*/
long tuple [ MAXATR ];  /* generated tuple values	*/
int ntupl;

/* -----------------------------------------------------*/

/* status variable	*/
char axis;		/* the previous split axis. Next will be
			 * axis+1 (mod nattrib) in case of cyclic split.
			 */

/* -----------------------------------------------------*/
char next;
FILE *ofp,*pfp,*sfp;

