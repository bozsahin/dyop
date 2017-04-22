#include <stdio.h>
#include <math.h>
#include "dyop.h"

/*
 *                      DYOP SIMULATOR
 *
 * -Cem Bozsahin
 *  Arizona State University
 *  1986
 */

main( argc, argv)
int  argc;
char *argv[];	/* arguments are parameter, dyop data and statistics
		 * file names .
		 */
{

	if (argc != 4){
		printf("dyop:usage is: fn parfile dfile statfile \n");
		exit(1);
	}
	if ((pfp=fopen(argv[1],"r")) == NULL ) {
		printf("dyop: can't open argument file ! \n");
		exit(1);
	}
        if ((ofp=fopen(argv[2],"w+")) == NULL ){
		printf("dyop: can't create data file ! \n");
		exit(1);
	}
	if ((sfp=fopen(argv[3],"w")) == NULL ) {
		printf("dyop: can't create summary file ! \n");
		exit(1);
	}
	fscanf(pfp,"%d %d %d %d",&seed,&maxtupl,&pagesize,&nattrib);
	
	init(argv[1]);

	/*
	 * main loop : simulate dyop for 'maxtupl' tuples
	 */
	next = 1;
	for (ntupl = 1; ntupl <= maxtupl;ntupl++){
		if(next)
			uniform();
		else
			ntupl--;
		next = 1;
		map();
	}
	stats();
}


/* 
 * set the parameters for the simulation : read in domain info
 * and find effective partition sizes
 */
init(ar1)
char *ar1;
{
	int i;

	for(i=0;i < nattrib;i++){
		fscanf(pfp,"%ld %d",&d[i].range,&d[i].size);
		tsize += d[i].size;
	}
	if (tsize > pagesize){
		printf("dyop: tuple size > page size . No way ! \n");
		exit(1);
	}
	b0 = pagesize/tsize;
	axis = nattrib-1;	/* 1st split is at 0th attrib.  */
	srandom(seed);		/* random # generator init.	*/
	stsize = nattrib*ATRS+1;/* simulator tuple size		*/
	spsize = stsize * b0;	/* simulator page size		*/
	/* 
	 * generate the very first partition ( # 0 )
	 */
	p[0].plow = 0;
	p[0].phigh = d[0].range;
 	for(i=1;i <= spsize ; i++)
		fprintf(ofp,"%c",' ');
	/* dump parameters	*/
	printf("\n\n dyop simulator \n\n");
	printf("argument file = %s\n",ar1);
	printf("# of tuples = %d pagesize = %d # of attr. = %d seed=%d\n",
		maxtupl,pagesize,nattrib,seed);
	printf("\n partition size in # of tuples = %d\n\n",b0);
	printf("attribute information :\n\n        range          size\n");
	for(i = 0; i < nattrib;i++)
		printf("%2d %10ld %10d\n",i,d[i].range,d[i].size);
}

/*
 * generate statistics	
 */
stats()
{
	int i;
        int waste;
	
	printf("\n\n\n dyop simulation statistics \n ------------------\n\n");
	printf("# of explicit partitions generated = %d\n",npart+1);
	printf("\n splitting information \n axis \t # of splits\n\n");
	for(i = 0; i < nattrib;i++)
		printf("%2d \t %3d\n",i,nsplit[i]);
	printf("\n average data space utilization = %8.3f %% \n",
		100*(double)maxtupl/((npart+1)*b0));
	waste = pagesize - ( b0*tsize);
	printf("\n waste of space due to truncation = %8.3f %% \n",
		100*(double)waste/pagesize);
	printf("\n total disk storage occupied by data = %d  bytes \n",
		(npart+1)*pagesize);
	for (i = 0;i <= npart; i++)
		fprintf(sfp,"%12ld %10d %12ld %3d %12ld %12ld \n",p[i].pno,
			p[i].ntup,p[i].paddr,p[i].splx,p[i].plow,p[i].phigh);
}

/*
 * uniform data generator
 */
uniform()
{
	int i;

	for(i = 0 ; i < nattrib;i++)
		tuple[i] = d[i].range * (double) random() / MAXVAL;
}

/* 
 * find explicit partition in the partition table
 * if not found return IMPLICIT
 * Binary search is impossible due to implicit partitions becoming
 * explicit after an overflow.
 */
fxp(px)
long px;
{
	int i;
 
	for(i = npart;i >= 0;i--)
		if( p[i].pno == px)
			return(i);
	return(IMPLICIT);
}

/* 
 * find the first implicit in 'm' to become explicit (i.e., implicit
 * partition with smallest # in 'm'
 */
long fmp(m)
long m;
{
	long mm;
	int splm;	/* Split # (in whole data space) partition m has 
			 * been generated.	
			 */

	if (tsplit == 1) return(1L);	/* #1 is in #0 at 1st split	*/
	splm = log10((double)(m+1))/log10(2.0);
	if (tsplit == splm) return (0L);/* No one embedded in.		*/
	mm = m;
	do
		mm += 01<<splm;		/* each part. is 2**splm apart. */
	while ((fxp(mm)) != IMPLICIT);
	return(mm);
}

/* 
 * find the explicit partition embedding m.
 * See formula (2) in reference
 */
embed(m,iptr)
long *m,*iptr;
{
	int ix,a;
	long xm,ixm,t;

	ix = IMPLICIT;
	xm = *m;
	ixm = xm;
	a = 0;
	if (xm > 0)
		a = log10((double) xm) / log10(2.0);
	while (xm >= 0 && (ix = fxp(xm)) == IMPLICIT && a >= 0){
		t = xm;
		if ( xm&01<<(a))
			xm -= 01<<a;
		if(t != xm) ixm = t;
		a--;
	}
	*iptr = ixm;
	*m = xm;
	return(ix);
}

/*
 * Called when the mapped explicit partition overflows.
 * Update the old partition and the new partition by moving
 * tuples to new partition. The newly generated partition is the implicit
 * partition with the smallest number in the overflowing partition.
 */
emerge(ixp,impp)
int ixp;	/* p table index of explicit part.	*/
long impp;	/* implicit p that will become explicit */
{
	int a,splno;
	char ax,x,tsta;
	long base1,base2;	/* partition offsets	*/
	long tbuf [ MAXATR ];	/* buffer for a tuple	*/
	long bit;
        int i,j;
        long low,high;

	/* find coordinates of emerging new partition	*/
	a = 0;
	if ( impp > 0)
		a = log10((double)impp)/log10(2.0);
	ax = a % nattrib;		/* split axis	*/
	splno = a/nattrib+1;		/* split #	*/
	/* find lower and upper bounds			*/
	x = ax;
	low = 0;
	high = d[ax].range;
	bit = impp & 01<<ax;
	for(i = 0; i < splno;i++){
		if (bit)
			low = high/2;
		else
			high /= 2;
		x += nattrib;
		bit = impp & 01<<x;
	}
        if (low == high){
		next=0;
		return;
	}
	/* Allocate a page for new partition		*/
	p[++npart].pno = impp;	/* make it explicit	*/
	p[npart].paddr = npart * b0 * stsize;
	fseek(ofp,0L,2);
	for(i=1; i <= spsize;i++)
		fprintf(ofp,"%c",' ');
        p[npart].plow = low;
	p[npart].phigh = high;
	p[npart].splx = ax;
	/* If the old & new partition were generated at the same axis,
 	 * update the upper bound of old partition (old one always in
	 * the first half of the partition domain).
	 */
	if (ax == p[ixp].splx)
		p[ixp].phigh = p[npart].plow;
	base1 = p[ixp].paddr;
	base2 = p[npart].paddr;
	/* copy the tuples into new partition on the basis
	 * of split axis values.
	 * Physical write of the tuple will occur in case of 
	 * "no overflow"  in map. 
	 * The only thing exchanged here is the existing tuples.
	 */
	a = p[ixp].ntup;
	for (i = 1; i <= a;){
		fseek(ofp,base1,0);
		fscanf(ofp,"%c",&tsta);
		if (tsta&TSTORED){
			for(j = 0; j < nattrib;j++)
				fscanf(ofp,"%10ld",&tbuf[j]);
			if(tbuf[ax] >= p[npart].plow &&
			   tbuf[ax] <= p[npart].phigh){
				fseek(ofp,base1,0);
				fprintf(ofp,"%c",TMOVED); /* moved tuple*/
				fseek(ofp,base2,0);
				fprintf(ofp,"%c",TSTORED);
				for(j=0;j<nattrib;j++)
					fprintf(ofp,"%10ld",tbuf[j]);
				base2 += stsize;
				p[npart].ntup++;
				p[ixp].ntup--;
			}
			i++;
		}
		base1 += stsize;
	}
}

/*
 * partition mapping routine for dyop : see the reference
 */
map()
{
	char c;
	long m;		/* part. # that may contain the tuple	*/
	long impp;
	int ix,a,i;
        long base;

	m=impp=0;
        m = hash();
	/* find the explicit partition # */
	ix = embed(&m,&impp);
	if (ix == IMPLICIT ){
	     printf("dyop:fatal error.Partition not found\n");
	     exit(1);
	}
	/* if explicit partition overflows, update old and new
	 * partition. If it doesn't, then enter the tuple into the part.
	 */
	if ( p[ix].ntup == b0){
		next = 0;			/*  Tuple not to be written
						 *  in this pass. Try again.
						 */
		if (m != impp)
			emerge(ix,impp);	/* make impl. explicit	*/
		else {				/* no one embedded in.  */
			axis = (axis+1)%nattrib;/* cyclic split.	*/
			nsplit[axis]++;		/* have to split.	*/
			tsplit++;		/* now that there is an
						 * implicit partition in it:
						 */
			/* Find which one to become explicit.		*/
			impp = fmp(m);
			emerge(ix,impp);	
		}
	}
	else {
		/*  no overflow. Find a free slot in p			*/
		base = p[ix].paddr;
		fseek(ofp,base,0);
		fscanf(ofp,"%c",&c);
		while ( c&TSTORED){
			base += stsize;
			fseek(ofp,base,0);
			fscanf(ofp,"%c",&c);
		}
		fseek(ofp,base,0);
		fprintf(ofp,"%c",TSTORED);	/* fill the slot	*/
		for(i=0;i<nattrib;i++)
			fprintf(ofp,"%10ld",tuple[i]);
		p[ix].ntup++;
	}
}

/*	Find the partition # ( implicit or explicit ) using formula 1
 * 	in reference. 
 * 	Essentially, this is Multidimensional Linear Dynamic Hashing.
 */
long hash ()
{
	long di [ MAXATR ];	/* Di values in reference.	*/
	long z,m;
	int i,j;

	m = 0;
	/* find di = Vi / |Di|*2**Li for all attributes		*/
	for (i = 0; i < nattrib; i++)
		di[i] = (tuple [i]<<nsplit[i]) /(double) d[i].range;
	/* find the partition # using formula (1) in reference	*/
	for(i=0;i<nattrib;i++){
		z = 0;
		for (j = 0;j < nsplit[i];j++)
			if(di[i]&(01<<j))
				z += 01<<(nattrib*(nsplit[i]-1-j));
		m += (01<<i)*z;
	}
        return(m);
}
