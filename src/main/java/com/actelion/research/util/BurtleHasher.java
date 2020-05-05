package com.actelion.research.util;

import java.util.Date;

/**
 * http://burtleburtle.net/bob/c/lookup3.c
 -------------------------------------------------------------------------------
	lookup3.c, by Bob Jenkins, May 2006, Public Domain.

	These are functions for producing 32-bit hashes for hash table lookup.
	hashword(), hashlittle(), hashbig(), mix(), and final() are externally 
	useful functions.  Routines to test the hash are included if SELF_TEST 
	is defined.  You can use this free for any purpose.  It has no warranty.

	You probably want to use hashlittle().  hashlittle() and hashbig()
	hash byte arrays.  hashlittle() is is faster than hashbig() on
	little-endian machines.  Intel and AMD are little-endian machines.

	If you want to find a hash of, say, exactly 7 integers, do
	  a = i1;  b = i2;  c = i3;
	  mix(a,b,c);
	  a += i4; b += i5; c += i6;
	  mix(a,b,c);
	  a += i7;
	  final(a,b,c);
	then use c as the hash value.  If you have a variable length array of
	4-byte integers to hash, use hashword().  If you have a byte array (like
	a character string), use hashlittle().  If you have several byte arrays, or
	a mix of things, see the comments above hashlittle().
	-------------------------------------------------------------------------------
	
 */
public class BurtleHasher {
	
	// The static constructor prevents many new calls
	private static BurtleHasherABC abcHashlittleInteger = new BurtleHasherABC(0,0,0); 
	
	public static int hashsize(long n) {
		int v = 1<<(n);
		return v;
	}

	/**
	 *
	 * @param n number of bits set in mask.
	 * @return
	 */
	public static int hashmask(int n)  {
		int v = (hashsize(n)-1);
		
		return v;
	}
	public static long rot(long x, long k) {
		long v = (((x)<<(k)) ^ ((x)>>>(32-(k))));
		return v;
	}

	/*
	-------------------------------------------------------------------------------
	mix -- mix 3 32-bit values reversibly.

	This is reversible, so any information in (a,b,c) before mix() is
	still in (a,b,c) after mix().

	If four pairs of (a,b,c) inputs are run through mix(), or through
	mix() in reverse, there are at least 32 bits of the output that
	are sometimes the same for one pair and different for another pair.
	This was tested for:
	* pairs that differed by one bit, by two bits, in any combination
	  of top bits of (a,b,c), or in any combination of bottom bits of
	  (a,b,c).
	* "differ" is defined as +, -, ^, or ~^.  For + and -, I transformed
	  the output delta to a Gray code (a^(a>>1)) so a string of 1's (as
	  is commonly produced by subtraction) look like a single 1-bit
	  difference.
	* the base values were pseudorandom, all zero but one bit set, or 
	  all zero plus a counter that starts at zero.

	Some k values for my "a-=c; a^=rot(c,k); c+=b;" arrangement that
	satisfy this are
	    4  6  8 16 19  4
	    9 15  3 18 27 15
	   14  9  3  7 17  3
	Well, "9 15 3 18 27 15" didn't quite get 32 bits diffing
	for "differ" defined as + with a one-bit base and a two-bit delta.  I
	used http://burtleburtle.net/bob/hash/avalanche.html to choose 
	the operations, constants, and arrangements of the variables.

	This does not achieve avalanche.  There are input bits of (a,b,c)
	that fail to affect some output bits of (a,b,c), especially of a.  The
	most thoroughly mixed value is c, but it doesn't really even achieve
	avalanche in c.

	This allows some parallelism.  Read-after-writes are good at doubling
	the number of bits affected, so the goal of mixing pulls in the opposite
	direction as the goal of parallelism.  I did what I could.  Rotates
	seem to cost as much as shifts on every machine I could lay my hands
	on, and rotates are much kinder to the top and bottom bits, so I used
	rotates.
	-------------------------------------------------------------------------------
	*/
	private static void mix(BurtleHasherABC abc) {
		long a = abc.a;
		long b = abc.b;
		long c = abc.c;

		a -= c;
		a ^= rot(c, 4);
		c += b;
		b -= a;
		b ^= rot(a, 6);
		a += c;
		c -= b;
		c ^= rot(b, 8);
		b += a;
		a -= c;
		a ^= rot(c, 16);
		c += b;
		b -= a;
		b ^= rot(a, 19);
		a += c;
		c -= b;
		c ^= rot(b, 4);
		b += a;

		abc.a = a;
		abc.b = b;
		abc.c = c;

	}
	
	/*
	--------------------------------------------------------------------
	mix -- mix 3 64-bit values reversibly.
	mix() takes 48 machine instructions, but only 24 cycles on a superscalar
	  machine (like Intel's new MMX architecture).  It requires 4 64-bit
	  registers for 4::2 parallelism.
	All 1-bit deltas, all 2-bit deltas, all deltas composed of top bits of
	  (a,b,c), and all deltas of bottom bits were tested.  All deltas were
	  tested both on random keys and on keys that were nearly all zero.
	  These deltas all cause every bit of c to change between 1/3 and 2/3
	  of the time (well, only 113/400 to 287/400 of the time for some
	  2-bit delta).  These deltas all cause at least 80 bits to change
	  among (a,b,c) when the mix is run either forward or backward (yes it
	  is reversible).
	This implies that a hash using mix64 has no funnels.  There may be
	  characteristics with 3-bit deltas or bigger, I didn't test for
	  those.
	--------------------------------------------------------------------
	*/
	public static void mix64(BurtleHasherABC abc) {
		
		abc.a = abc.a - abc.b;
		abc.a = abc.a - abc.c;
		abc.a = abc.a ^ (abc.c >> 43);
		abc.b = abc.b - abc.c;
		abc.b = abc.b - abc.a;
		abc.b = abc.b ^ (abc.a << 9);
		abc.c = abc.c - abc.a;
		abc.c = abc.c - abc.b;
		abc.c = abc.c ^ (abc.b >> 8);
		abc.a = abc.a - abc.b;
		abc.a = abc.a - abc.c;
		abc.a = abc.a ^ (abc.c >> 38);
		abc.b = abc.b - abc.c;
		abc.b = abc.b - abc.a;
		abc.b = abc.b ^ (abc.a << 23);
		abc.c = abc.c - abc.a;
		abc.c = abc.c - abc.b;
		abc.c = abc.c ^ (abc.b >> 5);
		abc.a = abc.a - abc.b;
		abc.a = abc.a - abc.c;
		abc.a = abc.a ^ (abc.c >> 35);
		abc.b = abc.b - abc.c;
		abc.b = abc.b - abc.a;
		abc.b = abc.b ^ (abc.a << 49);
		abc.c = abc.c - abc.a;
		abc.c = abc.c - abc.b;
		abc.c = abc.c ^ (abc.b >> 11);
		abc.a = abc.a - abc.b;
		abc.a = abc.a - abc.c;
		abc.a = abc.a ^ (abc.c >> 12);
		abc.b = abc.b - abc.c;
		abc.b = abc.b - abc.a;
		abc.b = abc.b ^ (abc.a << 18);
		abc.c = abc.c - abc.a;
		abc.c = abc.c - abc.b;
		abc.c = abc.c ^ (abc.b >> 22);
	}

	/*
	 * -------------------------------------------------------------------------------
	 * final -- final mixing of 3 32-bit values (a,b,c) into c
	 * 
	 * Pairs of (a,b,c) values differing in only a few bits will usually produce
	 * values of c that look totally different. This was tested for pairs that
	 * differed by one bit, by two bits, in any combination of top bits of
	 * (a,b,c), or in any combination of bottom bits of (a,b,c). "differ" is
	 * defined as +, -, ^, or ~^. For + and -, I transformed the output delta to
	 * a Gray code (a^(a>>1)) so a string of 1's (as is commonly produced by
	 * subtraction) look like a single 1-bit difference. the base values were
	 * pseudorandom, all zero but one bit set, or all zero plus a counter that
	 * starts at zero.
	 * 
	 * These constants passed: 14 11 25 16 4 14 24 12 14 25 16 4 14 24 and these
	 * came close: 4 8 15 26 3 22 24 10 8 15 26 3 22 24 11 8 15 26 3 22 24
	 * -------------------------------------------------------------------------------
	 */
	private static void finalMix(BurtleHasherABC abc) 
	{ 
		long a = abc.a;
		long b = abc.b;
		long c = abc.c;
		
		
	  c ^= b; 
	  c -= rot(b,14); 
	  a ^= c; 
	  a -= rot(c,11); 
	  b ^= a; 
	  b -= rot(a,25); 
	  c ^= b; 
	  c -= rot(b,16); 
	  a ^= c; 
	  a -= rot(c,4);  
	  b ^= a; 
	  b -= rot(a,14); 
	  c ^= b; 
	  c -= rot(b,24); 
	  
		abc.a = a;
		abc.b = b;
		abc.c = c;
	}

	/*
	--------------------------------------------------------------------
	 This works on all machines.  To be useful, it requires
	 -- that the key be an array of uint32_t's, and
	 -- that all your machines have the same endianness, and
	 -- that the length be the number of uint32_t's in the key

	 The function hashword() is identical to hashlittle() on little-endian
	 machines, and identical to hashbig() on big-endian machines,
	 except that the length has to be measured in uint32_ts rather than in
	 bytes.  hashlittle() is more complicated than hashword() only because
	 hashlittle() has to dance around fitting the key bytes into registers.
	--------------------------------------------------------------------
	*/
	public static int hashword(String w, long initval)
	{
		// uint32_t *k;                         /* the key, an array of uint32_t values */
		// size_t  length;                       /* the length of the key, in uint32_ts */
		// uint32_t  initval;               /* the previous hash, or an arbitrary value */

		
		int length = w.length();
		
		byte [] k = w.getBytes();
		long a,b,c;

	  /* Set up the internal state */
	  a = b = c = 0xdeadbeef + (((long)length)<<2) + initval;

	  BurtleHasherABC abc = new BurtleHasherABC(a,b,c);
	  
	  
	  /*------------------------------------------------- handle most of the key */
	  int cc = 0;
	  while (length > 3)
	  {
		  abc.a += k[cc+0];
		  abc.b += k[cc+1];
		  abc.c += k[cc+2];
	    mix(abc);
	    length -= 3;
	    cc += 3;
	  }

	  /*------------------------------------------- handle the last 3 uint32_t's */
	  switch(length)                     /* all the case statements fall through */
	  { 
	  case 3 : abc.c+=k[2];
	  case 2 : abc.b+=k[1];
	  case 1 : abc.a+=k[0];
	    finalMix(abc);
	  case 0:     /* case 0: nothing left to add */
	    break;
	  }
	  /*------------------------------------------------------ report the result */
	  return (int)abc.c;
	}


	/*
	-------------------------------------------------------------------------------
	hashlittle() -- hash a variable-length key into a 32-bit value
	  k       : the key (the unaligned variable-length array of bytes)
	  length  : the length of the key, counting by bytes
	  initval : can be any 4-byte value
	Returns a 32-bit value.  Every bit of the key affects every bit of
	the return value.  Two keys differing by one or two bits will have
	totally different hash values.

	The best hash table sizes are powers of 2.  There is no need to do
	mod a prime (mod is sooo slow!).  If you need less than 32 bits,
	use a bitmask.  For example, if you need only 10 bits, do
	  h = (h & hashmask(10));
	In which case, the hash table should have hashsize(10) elements.

	If you are hashing n strings (uint8_t **)k, do it like this:
	  for (i=0, h=0; i<n; ++i) h = hashlittle( k[i], len[i], h);

	By Bob Jenkins, 2006.  bob_jenkins@burtleburtle.net.  You may use this
	code any way you wish, private, educational, or commercial.  It's free.

	Use for hash table lookup, or anything where one collision in 2^^32 is
	acceptable.  Do NOT use for cryptographic purposes.
	-------------------------------------------------------------------------------
	*/

	public static int hashlittle(String key, long initval)
	{
	  long a,b,c;

	  /* Set up the internal state */
	  
	  int length = key.length();
	  
	  a = b = c = 0xdeadbeef + ((long)length) + initval;

	                     /* need to read the key one byte at a time */
	  byte [] k = key.getBytes();

	  BurtleHasherABC abc = new BurtleHasherABC(a,b,c); 
	  
	  int cc=0;
	    /*--------------- all but the last block: affect some 32 bits of (a,b,c) */
	    while (length > 12)
	    {
	    	abc.a += k[cc+0];
	    	abc.a += ((long)k[cc+1])<<8;
	    	abc.a += ((long)k[cc+2])<<16;
	    	abc.a += ((long)k[cc+3])<<24;
	    	abc.b += k[cc+4];
	    	abc.b += ((long)k[cc+5])<<8;
	    	abc.b += ((long)k[cc+6])<<16;
	    	abc.b += ((long)k[cc+7])<<24;
	    	abc.c += k[cc+8];
	    	abc.c += ((long)k[cc+9])<<8;
	    	abc.c += ((long)k[cc+10])<<16;
	    	abc.c += ((long)k[cc+11])<<24;
	      mix(abc);
	      length -= 12;
	      cc += 12;
	    }

	    /*-------------------------------- last block: affect all 32 bits of (c) */
	    switch(length)                   /* all the case statements fall through */
	    {
	    case 12: abc.c+=((long)k[cc + 11])<<24;
	    case 11: abc.c+=((long)k[cc + 10])<<16;
	    case 10: abc.c+=((long)k[cc + 9])<<8;
	    case 9 : abc.c+=k[8];
	    case 8 : abc.b+=((long)k[cc + 7])<<24;
	    case 7 : abc.b+=((long)k[cc + 6])<<16;
	    case 6 : abc.b+=((long)k[cc + 5])<<8;
	    case 5 : abc.b+=k[4];
	    case 4 : abc.a+=((long)k[cc + 3])<<24;
	    case 3 : abc.a+=((long)k[cc + 2])<<16;
	    case 2 : abc.a+=((long)k[cc + 1])<<8;
	    case 1 : abc.a+=k[cc + 0];
	             break;
	    case 0 : return (int)abc.c;
	    }
	  

	  finalMix(abc);
	  return (int)abc.c;
	}
	public static int hashlittle(byte [] k, long initval) {
		return hashlittle(k, initval, k.length);
	}

	/**
	 *
	 * @param k
	 * @param initval
	 * @param size the hash value will be calculated for k up to size fields.
	 * @return
	 */
	public static int hashlittle(byte [] k, long initval, int size) {
		
	  long a,b,c;

	  /* Set up the internal state */
	  
	  int length = size;
	  
	  a = b = c = 0xdeadbeef + ((long)length) + initval;

	  /* need to read the key one byte at a time */
	  BurtleHasherABC abc = new BurtleHasherABC(a,b,c); 
	  
	  int cc=0;
	    /*--------------- all but the last block: affect some 32 bits of (a,b,c) */
	    while (length > 12)
	    {
	    	abc.a += k[cc+0];
	    	abc.a += ((long)k[cc+1])<<8;
	    	abc.a += ((long)k[cc+2])<<16;
	    	abc.a += ((long)k[cc+3])<<24;
	    	abc.b += k[cc+4];
	    	abc.b += ((long)k[cc+5])<<8;
	    	abc.b += ((long)k[cc+6])<<16;
	    	abc.b += ((long)k[cc+7])<<24;
	    	abc.c += k[cc+8];
	    	abc.c += ((long)k[cc+9])<<8;
	    	abc.c += ((long)k[cc+10])<<16;
	    	abc.c += ((long)k[cc+11])<<24;
	      mix(abc);
	      length -= 12;
	      cc += 12;
	    }

	    /*-------------------------------- last block: affect all 32 bits of (c) */
	    switch(length)                   /* all the case statements fall through */
	    {
	    case 12: abc.c+=((long)k[cc + 11])<<24;
	    case 11: abc.c+=((long)k[cc + 10])<<16;
	    case 10: abc.c+=((long)k[cc + 9])<<8;
	    case 9 : abc.c+=k[8];
	    case 8 : abc.b+=((long)k[cc + 7])<<24;
	    case 7 : abc.b+=((long)k[cc + 6])<<16;
	    case 6 : abc.b+=((long)k[cc + 5])<<8;
	    case 5 : abc.b+=k[4];
	    case 4 : abc.a+=((long)k[cc + 3])<<24;
	    case 3 : abc.a+=((long)k[cc + 2])<<16;
	    case 2 : abc.a+=((long)k[cc + 1])<<8;
	    case 1 : abc.a+=k[cc + 0];
	             break;
	    case 0 : return (int)abc.c;
	    }
	  

	  finalMix(abc);
	  return (int)abc.c;
	}
	
	
	/**
	 * MvK 10.02.2010
	 * @param k
	 * @param initval
	 * @return
	 */
	public static int hashlittle(int [] k, long initval) {
		return hashlittle(k, initval, k.length);
	}
	/**
	 * 
	 * @param k
	 * @param initval
	 * @param size the hash value will be calculated for k up to size fields. 
	 * @return
	 */
	public static int hashlittle(int [] k, long initval, int size) {
	  long a,b,c;

	  /* Set up the internal state */
	  
	  int length = size;
	  
	  a = b = c = 0xdeadbeef + ((long)length) + initval;

	  /* need to read the key one byte at a time */
	  abcHashlittleInteger.a = a;
	  abcHashlittleInteger.b = b;
	  abcHashlittleInteger.c = c;
	   
	  
	  int cc=0;
	    /*--------------- all but the last block: affect some 32 bits of (a,b,c) */
	    while (length > 12) {
	    	abcHashlittleInteger.a += k[cc+0];
	    	abcHashlittleInteger.a += ((long)k[cc+1])<<8;
	    	abcHashlittleInteger.a += ((long)k[cc+2])<<16;
	    	abcHashlittleInteger.a += ((long)k[cc+3])<<24;
	    	abcHashlittleInteger.b += k[cc+4];
	    	abcHashlittleInteger.b += ((long)k[cc+5])<<8;
	    	abcHashlittleInteger.b += ((long)k[cc+6])<<16;
	    	abcHashlittleInteger.b += ((long)k[cc+7])<<24;
	    	abcHashlittleInteger.c += k[cc+8];
	    	abcHashlittleInteger.c += ((long)k[cc+9])<<8;
	    	abcHashlittleInteger.c += ((long)k[cc+10])<<16;
	    	abcHashlittleInteger.c += ((long)k[cc+11])<<24;
	    	
	    	mix(abcHashlittleInteger);
	    	length -= 12;
	    	cc += 12;
	    }

	    /*-------------------------------- last block: affect all 32 bits of (c) */
	    switch(length)                   /* all the case statements fall through */
	    {
	    case 12: abcHashlittleInteger.c+=((long)k[cc + 11])<<24;
	    case 11: abcHashlittleInteger.c+=((long)k[cc + 10])<<16;
	    case 10: abcHashlittleInteger.c+=((long)k[cc + 9])<<8;
	    case 9 : abcHashlittleInteger.c+=k[8];
	    case 8 : abcHashlittleInteger.b+=((long)k[cc + 7])<<24;
	    case 7 : abcHashlittleInteger.b+=((long)k[cc + 6])<<16;
	    case 6 : abcHashlittleInteger.b+=((long)k[cc + 5])<<8;
	    case 5 : abcHashlittleInteger.b+=k[4];
	    case 4 : abcHashlittleInteger.a+=((long)k[cc + 3])<<24;
	    case 3 : abcHashlittleInteger.a+=((long)k[cc + 2])<<16;
	    case 2 : abcHashlittleInteger.a+=((long)k[cc + 1])<<8;
	    case 1 : abcHashlittleInteger.a+=k[cc + 0];
	             break;
	    case 0 : return (int)abcHashlittleInteger.c;
	    }
	  

	    finalMix(abcHashlittleInteger);
	  
	    return (int)abcHashlittleInteger.c;
	}





	/* used for timings */
	public static void driver1()
	{
	  int size = 256;
	  

	  Date daStart = new Date();
	  String s = "";
	  for (int i=0; i<size; i++) 
		  s += 'x';
	  long h = 0;
	  for (int i=0; i<1; ++i) 
	  {
	    h = hashlittle(s,h);
	  }
	  Date daEnd = new Date();
	  
	  long delta = daEnd.getTime() - daStart.getTime();
	  System.out.println("time: " + delta);
	}

	/* check that every input bit changes every output bit half the time */
	private final static int HASHSTATE = 1;
	private final static long HASHLEN = 1;
	private final static long MAXPAIR = 60;
	private final static int MAXLEN = 70;
	private static void driver2()
	{
		/*
	  char [] qa = new char[MAXLEN+1]; 
	  char [] qb = new char[MAXLEN+2]; 
	  int qaCC = 0; 
	  int qbCC = 1; 
	  // *b = &qb[1];
	  long [] c = new long [HASHSTATE];
	  long [] d = new long [HASHSTATE];
	  int i=0, j=0, k, l, m=0, z=0;
	  long [] e = new long [HASHSTATE];
	  long [] f = new long [HASHSTATE];
	  long [] g = new long [HASHSTATE];
	  long [] h = new long [HASHSTATE];
	  long [] x = new long [HASHSTATE];
	  long [] y = new long [HASHSTATE];
	  long hlen = 0;

	  System.out.println("No more than %d trials should ever be needed \n" + MAXPAIR/2);
	  boolean bFin = false;
	  for (hlen=0; hlen < MAXLEN; ++hlen)
	  {
	    z=0;
	    for (i=0; i<hlen; ++i)  //----------------------- for each input byte
	    {
	      for (j=0; j<8; ++j)   //------------------------ for each input bit, 
	      {
		for (m=1; m<8; ++m) //------------ for serveral possible initvals
		{
		  for (l=0; l<HASHSTATE; ++l)
		    e[l]=f[l]=g[l]=h[l]=x[l]=y[l]=~((long)0);

	      	  //---- check that every output bit is affected by that input bit
		  for (k=0; k<MAXPAIR; k+=2)
		  { 
		    long finished=1;
		    // keys have one bit different
		    for (l=0; l<hlen+1; ++l) {a[l] = b[l] = (char)0;}
		    //* have a and b be two keys differing in only one bit 
		    a[i] ^= (k<<j);
		    a[i] ^= (k>>(8-j));
		     c[0] = hashlittle(a, hlen, m);
		    b[i] ^= ((k+1)<<j);
		    b[i] ^= ((k+1)>>(8-j));
		     d[0] = hashlittle(b, hlen, m);
		    // check every bit is 1, 0, set, and not set at least once 
		    for (l=0; l<HASHSTATE; ++l)
		    {
		      e[l] &= (c[l]^d[l]);
		      f[l] &= ~(c[l]^d[l]);
		      g[l] &= c[l];
		      h[l] &= ~c[l];
		      x[l] &= d[l];
		      y[l] &= ~d[l];
		      if (e[l]|f[l]|g[l]|h[l]|x[l]|y[l]) finished=0;
		    }
		    if (finished) break;
		  }
		  if (k>z) z=k;
		  if (k==MAXPAIR) 
		  {
		     printf("Some bit didn't change: ");
		     printf("%.8x %.8x %.8x %.8x %.8x %.8x  ",
		            e[0],f[0],g[0],h[0],x[0],y[0]);
		     printf("i %ld j %ld m %ld len %ld\n",i,j,m,hlen);
		  }
		  if (z==MAXPAIR) {
			  bFin = true;
			  break;
			  }
		  
		  
		  
		}
		if(bFin)
			break;
	      }
			if(bFin)
				break;
	    }
	    
	    if(bFin)
	    if (z < MAXPAIR)
	    {
	      printf("Mix success  %2ld bytes  %2ld initvals  ",i,m);
	      printf("required  %ld  trials\n",z/2);
	    }
	  }
	  printf("\n");
	  */
	}

	/* Check for reading beyond the end of the buffer and alignment problems */
	public static void driver3()
	{
	  // uint8_t q[] = "This is the time for all good men to come to the aid of their country...";
	  String q = "This is the time for all good men to come to the aid of their country...";
	  
	  // uint8_t qq[] = "xThis is the time for all good men to come to the aid of their country...";
	  String qq = "xThis is the time for all good men to come to the aid of their country...";
	  // uint8_t qqq[] = "xxThis is the time for all good men to come to the aid of their country...";
	  String qqq = "xxThis is the time for all good men to come to the aid of their country...";
	  // uint8_t qqqq[] = "xxxThis is the time for all good men to come to the aid of their country...";
	  String qqqq = "xxxThis is the time for all good men to come to the aid of their country...";
	  
	  long h,i,j,ref,x,y;
	  // uint8_t *p;

	  System.out.println("Endianness.  These lines should all be the same (for values filled in):\n");
	  System.out.println("\n");
	  // System.out.println(hashword(q, ((q.length())-1)/4, 13));
	  // System.out.println(hashword(q, ((q.length())-5)/4, 13));
	  // System.out.println(hashword(q, ((q.length())-9)/4, 13)); 
	  System.out.println(hashword(q, 13));

	  String p1 = "";
	  String p2 = "";
	  
	  
	  
	  p1 = q.substring(0, q.length()-1);
	  p2 = q.substring(0, q.length()-2);
	  System.out.println(hashlittle(p1, 13) + " " +  hashlittle(p2, 13));
	  p1 = q.substring(0, q.length()-3);
	  p2 = q.substring(0, q.length()-4);
	  System.out.println(hashlittle(p1, 13) + " " +  hashlittle(p2, 13));
	  p1 = q.substring(0, q.length()-5);
	  p2 = q.substring(0, q.length()-6);
	  System.out.println(hashlittle(p1, 13) + " " +  hashlittle(p2, 13));
	  p1 = q.substring(0, q.length()-7);
	  p2 = q.substring(0, q.length()-8);
	  System.out.println(hashlittle(p1, 13) + " " +  hashlittle(p2, 13));
	  p1 = q.substring(0, q.length()-9);
	  p2 = q.substring(0, q.length()-10);
	  System.out.println(hashlittle(p1, 13) + " " +  hashlittle(p2, 13));
	  p1 = q.substring(0, q.length()-11);
	  p2 = q.substring(0, q.length()-12);
	  System.out.println(hashlittle(p1, 13) + " " +  hashlittle(p2, 13));
	  /*
	  p = &qq[1];
	  printf("%.8x %.8x %.8x %.8x %.8x %.8x %.8x %.8x %.8x %.8x %.8x %.8x\n",
	         hashlittle(p, sizeof(q)-1, 13), hashlittle(p, sizeof(q)-2, 13),
	         hashlittle(p, sizeof(q)-3, 13), hashlittle(p, sizeof(q)-4, 13),
	         hashlittle(p, sizeof(q)-5, 13), hashlittle(p, sizeof(q)-6, 13),
	         hashlittle(p, sizeof(q)-7, 13), hashlittle(p, sizeof(q)-8, 13),
	         hashlittle(p, sizeof(q)-9, 13), hashlittle(p, sizeof(q)-10, 13),
	         hashlittle(p, sizeof(q)-11, 13), hashlittle(p, sizeof(q)-12, 13));
	  p = &qqq[2];
	  printf("%.8x %.8x %.8x %.8x %.8x %.8x %.8x %.8x %.8x %.8x %.8x %.8x\n",
	         hashlittle(p, sizeof(q)-1, 13), hashlittle(p, sizeof(q)-2, 13),
	         hashlittle(p, sizeof(q)-3, 13), hashlittle(p, sizeof(q)-4, 13),
	         hashlittle(p, sizeof(q)-5, 13), hashlittle(p, sizeof(q)-6, 13),
	         hashlittle(p, sizeof(q)-7, 13), hashlittle(p, sizeof(q)-8, 13),
	         hashlittle(p, sizeof(q)-9, 13), hashlittle(p, sizeof(q)-10, 13),
	         hashlittle(p, sizeof(q)-11, 13), hashlittle(p, sizeof(q)-12, 13));
	  p = &qqqq[3];
	  printf("%.8x %.8x %.8x %.8x %.8x %.8x %.8x %.8x %.8x %.8x %.8x %.8x\n",
	         hashlittle(p, sizeof(q)-1, 13), hashlittle(p, sizeof(q)-2, 13),
	         hashlittle(p, sizeof(q)-3, 13), hashlittle(p, sizeof(q)-4, 13),
	         hashlittle(p, sizeof(q)-5, 13), hashlittle(p, sizeof(q)-6, 13),
	         hashlittle(p, sizeof(q)-7, 13), hashlittle(p, sizeof(q)-8, 13),
	         hashlittle(p, sizeof(q)-9, 13), hashlittle(p, sizeof(q)-10, 13),
	         hashlittle(p, sizeof(q)-11, 13), hashlittle(p, sizeof(q)-12, 13));
	        */
	  System.out.println("\n");
/*	  
	  for (h=0, b=buf+1; h<8; ++h, ++b)
	  {
	    for (i=0; i<MAXLEN; ++i)
	    {
	      len = i;
	      for (j=0; j<i; ++j) *(b+j)=0;

	      // these should all be equal 
	      ref = hashlittle(b, len, (uint32_t)1);
	      *(b+i)=(uint8_t)~0;
	      *(b-1)=(uint8_t)~0;
	      x = hashlittle(b, len, (uint32_t)1);
	      y = hashlittle(b, len, (uint32_t)1);
	      if ((ref != x) || (ref != y)) 
	      {
		printf("alignment error: %.8x %.8x %.8x %ld %ld\n",ref,x,y,h,i);
	      }
	    }
	  }
*/
	}

	/* check for problems with nulls */
/*	
	 void driver4()
	{
	  uint8_t buf[1];
	  uint32_t h,i,state[HASHSTATE];


	  buf[0] = ~0;
	  for (i=0; i<HASHSTATE; ++i) state[i] = 1;
	  printf("These should all be different\n");
	  for (i=0, h=0; i<8; ++i)
	  {
	    h = hashlittle(buf, 0, h);
	    printf("%2ld  0-byte strings, hash is  %.8x\n", i, h);
	  }
	}
*/

	public static void main(String [] args) 
	{
		
		String id1 = "DasIstderJanua";
		
		int hash = hashlittle(id1, 13);
		
		System.out.println("hash1: " + hash);

		hash = (hash&hashmask(10));
		
		System.out.println("hash2: " + hash);
		
		
	  // driver1();   /* test that the key is hashed: used for timings */
//		 driver2();   /* test that whole key is hashed thoroughly */
//		 driver3();   /* test that nothing but the key is hashed */
	  // driver4();   /* test hashing multiple buffers (all buffers are null) */
	 
	}

	


}
