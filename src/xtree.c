// Gabriel Al-Ghalith. Efficient characterization of orthogonal metagenomic
//  taxonomy and pathway coverage with CrossTree. 2018.
#define VER "CrossTree v2.00c 'ADAMANT WEAKSAUCE' by Gabe"
#define VNO 5
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <stdint.h>
#ifndef _OPENMP
	#define omp_get_max_threads() 1
	#define omp_get_thread_num() 0
	#include <time.h>
	#define omp_get_wtime() ((double)clock()/CLOCKS_PER_SEC)
	#define omp_set_num_threads(x) x
#endif
#include <sys/mman.h>
#include <sys/stat.h>
#ifndef kmer_t
	#define kmer_t uint32_t
#endif
#ifndef PFXLEN
	#define PFXLEN 13
#endif
#define rix_t uint32_t
#define AMBIG 4
#include <math.h>
#include <zlib.h>

void * huge_malloc(size_t n) {
	void *ptr = 0;
	posix_memalign(&ptr, 1 << 21, n);
	madvise(ptr, n, MADV_HUGEPAGE);
	return ptr;
}
void * huge_calloc(size_t n) {
	void *ptr = huge_malloc(n);
	memset(ptr,0,n);
	return ptr;
}
const uint8_t CONV[32]  = {4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4};
const uint8_t RCONV[32] = {4,3,4,2,4,4,4,1,4,4,4,4,4,4,4,4,4,4,4,4,0,0,4,4,4,4,4,4,4,4,4,4};
const uint8_t RCONS[32] = {'N','T','N','G','N','N','N','C','N','N','N','N','N','N','N','N',
						   'N','N','N','N','A','A','N','N','N','N','N','N','N','N','N','N'};
uint8_t FLTF1_1[5] = {1,0,0,0,0}, FLTF2_1[5] = {0,1,1,1,0}, FLTF2_2[5] = {0,1,0,0,0};
uint8_t FLTR1_1[5] = {0,0,0,1,0}, FLTR2_1[5] = {1,1,1,0,0}, FLTR2_2[5] = {0,0,1,0,0};


#pragma pack(1)
typedef struct {
	kmer_t sfx;
	rix_t rix;
} KPod;

//typedef struct {uint32_t h1, h2;} hpair_t;

int strcmp_ptr(const void* a, const void* b)
	{ return strcmp(*(char**)a,*(char**)b); }

typedef union {uint32_t a[4]; __uint128_t n; uint8_t b[16];} MasterBin_t;
int binCmp(const void *a, const void *b) {
	MasterBin_t *A = *(MasterBin_t **)a, *B = *(MasterBin_t **)b;
	uint64_t valA = *(uint64_t *)(A->b+4), valB = *(uint64_t *)(B->b+4);
	return valA < valB ? -1 : valB < valA;
}

static inline kmer_t dna_to_num(char *s, int nl, int *err, kmer_t kshift) {
	kmer_t nib = 0, k = kshift;
	for (uint32_t i = 0; i < nl; ++i, k-=2) {
		kmer_t c = CONV[31 & s[i]];
		if (c == AMBIG) {*err = 1; return i;}
		nib |= c << k;
	}
	return nib;
}


// As a convenience function, let's also define the opposite of the above, num_to_dna
static inline char * num_to_dna(kmer_t num, int nl, char *s) {
	for (int i = nl-1; i >= 0; --i, num >>= 2)
		s[i] = "ACGT"[num & 3];
	return s;
}

// Key string delimited by null
static inline uint64_t binsearch_str(char **Strings, char *key, uint64_t N) {
	uint64_t lo = 0, hi = N;
	while (lo < hi) {
		uint64_t mid = lo + ((hi-lo) >> 1);
		int cmp = strcmp(key,Strings[mid]);
		if (cmp > 0) lo = mid+1;
		else if (cmp < 0) hi = mid;
		else return mid;
	}
	return -1;
}
// Key string delimited by tab, null, newline...
static inline uint64_t binsearch_str_d(char **Strings, char *key, uint64_t N) {
	uint64_t lo = 0, hi = N;
	while (lo < hi) {
		uint64_t mid = lo + ((hi-lo) >> 1);
		char *a = key, *b = Strings[mid];
		while (*b && *a==*b) ++a, ++b;
		if (!*b && (!*a || *a == '\n' || *a=='\t')) return mid;
		if (*a < *b) hi = mid;
		else lo = mid+1;
	}
	return -1;
}

// returns first possible match as if null existed at position key_len, matched as far as possible
// A helper function will translate position via Map_h* minus ranks missing. 
static inline uint64_t binsearch_str_L(char **Strings, char *key, uint64_t N, uint32_t key_len) {
	uint64_t lo = 0, hi = N;
	while (lo < hi) {
		uint64_t mid = lo + ((hi - lo) >> 1);
		int cmp = 0;
		//strcmp(key,Strings[mid]);
		uint32_t i = 0; char *ref = Strings[mid];
		for (; i < key_len && ref[i]; ++i) 
			if (key[i] != ref[i]) break;
		if (i == key_len /* && key[i] */ && ref[i]) cmp = -1;
		else cmp = key[i] - ref[i];
		//printf("Compared key:\n%s\nagainst ref:\n%s\nuntil pos %u: %d [lo=%lu, mid=%lu, hi=%lu]\n",
		//	key,ref,key_len,cmp,lo,mid,hi);
		//if (!cmp) return mid;
		if (cmp > 0) lo = mid + 1;
		else if (cmp < 0) hi = mid;
		else return mid;
	}
	return lo; 
}
static inline int cmpfunc(const void *a, const void *b) {
	return *(uint64_t*)a < *(uint64_t*)b ? -1 :
		*(uint64_t*)b < *(uint64_t*)a;
}
static inline int u16cmp(const void *a, const void *b) {
	return *(uint16_t*)a < *(uint16_t*)b ? -1 : 
		*(uint16_t*)b < *(uint16_t*)a;
}
// sort comparator like above but for kmer_t
static inline int kcmp(const void *a, const void *b) {
	kmer_t k1 = *(kmer_t *)a, k2 = *(kmer_t *)b;
	return k1 < k2 ? -1 : k2 < k1;
}


static inline int kpackcmp(const void *a, const void *b) {
	KPod *k1 = (KPod *)a, *k2 = (KPod *)b;
	if (k1->sfx < k2->sfx) return -1;
	if (k2->sfx < k1->sfx) return 1;
	if (k1->rix < k2->rix) return -1;
	if (k2->rix < k1->rix) return 1;
	return 0;
}

static inline uint64_t LBS_k(KPod *KP, uint64_t N, kmer_t k) {
	uint64_t L = 0, R = N;
	while (L < R) {
		uint64_t m = (L + R) >> 1;
		if (KP[m].sfx < k) L = m+1;
		else R = m;
	}
	return (L < N && KP[L].sfx == k) ? L : -1;
}

// Mirror of the above binary search but just for kmer_t 
static inline uint64_t skLBS(kmer_t *A, uint64_t N, kmer_t k) {
	uint64_t L = 0, R = N;
	while (L < R) {
		uint64_t m = (L + R) >> 1;
		if (A[m] < k) L = m+1;
		else R = m;
	}
	return (L < N && A[L] == k) ? L : -1;
}

// Parses a fasta or fastq-formatted file, storing data in QBucket and HBucket, and returning the number of queries read
static inline uint64_t get_queries(gzFile in, uint8_t **QBucket, char **HBucket, uint8_t *head, uint8_t *line, int qChunk, uint64_t szmax) {
	uint64_t nq = 0; uint8_t *QB_ptr = *QBucket; char *H_ptr = *HBucket;
	uint8_t *eol; int len;
	while (nq < qChunk && QB_ptr - *QBucket <= szmax) {
		if (!gzgets(in,head,(1 << 20)-1)) break;
		eol = strchr(head,'\n');
		*eol = 0; 
		len = eol-head+1;
		memcpy(H_ptr,head+1,len-1);
		HBucket[nq] = H_ptr;
		H_ptr += len-1;
		
		if (!gzgets(in,line,(1 << 28)-1)) break;
		eol = strchr(line,'\n');
		*eol = 0; 
		len = eol-line+1;
		memcpy(QB_ptr,line,len);
		QBucket[nq] = QB_ptr;
		QB_ptr += len;
		if (*head == '@' && (!gzgets(in,line,1 << 28) || 
			!gzgets(in,line,1 << 28))) break;
		++nq;
	}
	return nq;
}

#define USAGE1 "USAGE: xtree {BUILD,ALIGN} [options]\n  "
#define USAGE2 USAGE1 "Options for both BUILD and ALIGN, with args: {seqs,log-out,threads,db}\n"
#define USAGE3 USAGE2 "BUILD Options\n  With args: {map,comp,k,db-out} <arg>\n"
#define USAGE4 USAGE3 "ALIGN Options\n  With args: {confidence,perq-out,ref-out,tax-out,cov-out,orthog-out}\n"
#define USAGE USAGE4 "  Without args: {redistribute,shallow-lca,copymem,doforage,half-forage,no-adamantium}\n"
int main(int argc, char *argv[]) {
	puts(VER);
	char *dbPath = 0, *seqPath = 0, *ixPath = 0, *covPath = 0;
	char *perqPath = 0, *taxPath = 0, *orthogPath = 0, *refPath = 0, *logPath = 0;
	int doBuild = 0, threads = omp_get_max_threads(), comp = -1, kchoice = 0;
	int doFullLCA = 1, doRedist = 0, doFastRedist = 0, doCopyMem = 0, doForage = 0, 
		doHalfForage = 0, adamantium = 1;
	double conf = 0.33; // Reasonable default for compression lv 2
	uint32_t nUniqMatches = 0;
	
	// Robust parser
	for (int a = 1; a < argc; ++a) {
		if (!strcmp(argv[a],"BUILD")) doBuild = 1; 
		else if (!strcmp(argv[a],"--map")) ixPath = argv[++a];
		else if (!strcmp(argv[a],"--comp")) comp = atoi(argv[++a]); //I
		else if (!strcmp(argv[a],"--k")) kchoice = atoi(argv[++a]); //I
		
		else if (!strcmp(argv[a],"ALIGN")) doBuild = 0; 
		else if (!strcmp(argv[a],"--confidence")) {
			double ctemp = atof(argv[++a]);
			if (ctemp <= 1) conf = ctemp, printf("Setting confprop %f\n",conf);
			else nUniqMatches = ctemp, printf("Setting min uniq ref matches to %u\n",nUniqMatches);
		}
		else if (!strcmp(argv[a],"--perq-out")) perqPath = argv[++a];
		else if (!strcmp(argv[a],"--ref-out")) refPath = argv[++a];
		else if (!strcmp(argv[a],"--tax-out")) taxPath = argv[++a];
		else if (!strcmp(argv[a],"--cov-out")) covPath = argv[++a];
		else if (!strcmp(argv[a],"--orthog-out")) orthogPath = argv[++a];
		else if (!strcmp(argv[a],"--redistribute")) doRedist = 1; //NA
		else if (!strcmp(argv[a],"--fast-redistribute")) doRedist = 1, doFastRedist = 1; //NA
		else if (!strcmp(argv[a],"--shallow-lca")) doFullLCA = 0; //NA
		else if (!strcmp(argv[a],"--copymem")) doCopyMem = 1; //NA
		else if (!strcmp(argv[a],"--doforage")) doForage = 1; //NA
		else if (!strcmp(argv[a],"--half-forage")) doHalfForage = 1; //NA
		else if (!strcmp(argv[a],"--adamantium")) adamantium = 1; //NA
		else if (!strcmp(argv[a],"--no-adamantium")) adamantium = 0; //NA
	
		
		// Options for both BUILD and ALIGN
		else if (!strcmp(argv[a],"--seqs")) seqPath = argv[++a]; 
		else if (!strcmp(argv[a],"--log-out")) logPath = argv[++a];
		else if (!strcmp(argv[a],"--threads")) threads = atoi(argv[++a]);
		else if (!strcmp(argv[a],"--db") || !strcmp(argv[a],"--db-out")) 
			dbPath = argv[++a];
		
		else {printf("Unrecognized option: %s\n",argv[a]); exit(1);}
	}
	threads = threads > 4096 ? 4096 : threads;
	comp=comp>2?2:comp; 
	omp_set_num_threads(threads); 
	printf("Using %d thread(s)\n",threads);
	if (argc < 4) {puts(USAGE); exit(1);}

	uint8_t *FltF1 = FLTF1_1, *FltF2 = FLTF2_1;
	uint8_t *FltR1 = FLTR1_1, *FltR2 = FLTR2_1;
	
	if (doBuild) {
		uint32_t PL = PFXLEN, SL = sizeof(kmer_t)*4;
		uint64_t K = PL+SL;
		if (comp == -1) comp = 0; // Default compression level
		if (comp > 2) {printf("Bad compression level! [%d] Range 0-2\n",comp); exit(1);}
		printf("Setting compression level to %d\n",comp);
		if (comp == 2) FltF2 = FLTF2_2, FltR2 = FLTR2_2;
		
		if (kchoice) K = kchoice;
		SL = K - PL;
		if (K < PL || !SL || SL > sizeof(kmer_t)*4) {printf("Bad K! [%lu]\n",K); exit(1);}
		printf("Building DB with K=%lu [PL %d, SL %d]\n",K,PL,SL);
		
		uint32_t kpre_shf = PL*2-2;
		kmer_t kpst_shf = SL*2-2;
		uint32_t pre_bshf = 32-(PL*2), pre_bshf_2 = pre_bshf+2;
		kmer_t pst_bshf = sizeof(kmer_t)*8-(SL*2), pst_bshf_2 = pst_bshf+2;

		FILE *in = fopen(seqPath,"rb");
		if (!in) {printf("ERROR: bad FASTA input: %s\n",seqPath); exit(2);}
		struct stat sb; int fno = fileno(in); fstat(fno,&sb);
		uint64_t fsz = sb.st_size;
		char *Raw = mmap(0, fsz+16, PROT_READ, MAP_SHARED, fno, 0);
		madvise(Raw,fsz,MADV_SEQUENTIAL);
		
		// We need: num records, num valid k-mers (each!), header locs
		if (Raw[0] != '>') {puts("Uh oh. Input FASTA looks fishy."); exit(2);}
		
		// Set up the tank for input -- up to 1 billion refs allowed
		uint64_t nbins = 1<<(2*PL),
			*Offsets = malloc(((uint64_t)1<<30)*sizeof(*Offsets)),
			*Nibs = huge_calloc((nbins+1)*sizeof(*Nibs));
		if (!Offsets || !Nibs) {puts("ERROR:OOM Offsets"); exit(3);}
		uint32_t ns = 0;
		//uint32_t mask32 = ((uint64_t)1 << (2*PL)) - 1;
		uint32_t shifty = 32 - 2*PL; // shift to zero out extra letters in the prefix left over due to the type being larger than the prefix
		kmer_t shifty2 = sizeof(kmer_t)*8 - 2*SL; // shift to zero out extra letters in the suffix left over due to the type being larger than the suffix
		uint32_t rc_rshift_pfx = 2*(PL-1); // shift to put letter in first position of prefix
		kmer_t rc_rshift_sfx = 2*(SL-1); // shift to put letter in first position of suffix
		double wtime = omp_get_wtime();
		uint64_t n_prefices = 0;
		/// KMER PREFIX COUNTING LOOP. Read a fasta record, decompose in to k-mers, and increment the prefix counts
		#pragma omp parallel for reduction(+:n_prefices) 
		for (uint64_t z = 0; z < fsz; ++z) {
			// The trick is to ignore all new lines that begin with the header, so each iteration starts at a sequence and ends before the next header.
			if (Raw[z] > 64 && Raw[z-1] == '\n') { // the ascii code 64 is the '@' symbol and 62 is '>', so looking above this means we're not at a header
				uint32_t ix;
				#pragma omp atomic capture
				ix = ns++; 
				Offsets[ix] = z; // Store where each sequence starts
				char *seq = Raw + z; // seq is a pointer to the beginning of the sequence
				int64_t y = 0; while (seq[y] && seq[y] != '\n') ++y; // y is the last index in the sequence
				// Loop over the sequence. This logic is further expanded to handle suffices in the "KMER BIN FILLING LOOP" later
				int64_t lastAmbig = -1;
				uint32_t pfx = 0, pfxRC = 0; // The prefix of the current k-mer. We can omit the sfx because we're only tallying prefixes here!
				for (int64_t j = 0; j < y; ++j) { // Start at the prefix length, and go to the end of the sequence
					uint32_t c = 0;  // c is the character we're looking at
					if (j >= SL) c = CONV[31 & seq[j-SL]]; // Get binary conversion of the letter, lag by suffix len (we don't want to overlap it)
					if (c==AMBIG) lastAmbig = j-SL, c = 0; // Mark last ambiguous position
					pfx = pfx << 2 | c; // Shift the prefix left by 2 bits and tack the new character onto the end
					pfx = pfx << shifty >> shifty; // Remove leftover bits
					
					// Now get the next character for the reverse complement's prefix
					c = RCONV[31 & seq[j]]; // starts on the other end of the kmer, so no need to lag
					c = c==AMBIG ? 0 : c; // Ambiguity already marked, so just cast to 0 ('A') to allow it to keep going
					pfxRC = pfxRC >> 2 | (uint32_t)c << rc_rshift_pfx; // Shift the prefix right by 2 bits and tack the new character onto the beginning

					// Also scan the current letter for ambigs so as to maintain parity with full k-mer scan later
					if (CONV[31 & seq[j]]==AMBIG) lastAmbig = lastAmbig > j ? lastAmbig : j; // Mark last ambiguous position
				
					// DEBUG: Print the prefix
					/* char pfxstr[PL+1]; pfxstr[PL] = 0;
					char pfxstrRC[PL+1]; pfxstrRC[PL] = 0;
					//for (int i = PL-1; i >= 0; --i) pfxstr[i] = "ACGT"[3 & (pfx >> 2*i)];
					num_to_dna(pfx,PL,pfxstr);
					num_to_dna(pfxRC,PL,pfxstrRC);
					printf("Ref %u, pos %u, pfx %s; pfxRC %s\n",ix,j-PL,pfxstr,pfxstrRC); */

					// Take the minimum of the two prefixes (normal and reverse complement)
					uint32_t final_pfx = pfx; 
					if (pfxRC < pfx) final_pfx = pfxRC;
					
					// If the read head is past the last ambiguous position, our current k-mer is complete, so look it up
					if (j >= lastAmbig + (int64_t)K) { 	
						if (!comp || j > K && FltF1[CONV[31 & seq[j-K]]] && FltF2[CONV[31 & seq[j-K-1]]] || 
								j < y - 1 && FltR1[CONV[31 & seq[j+1]]] && FltR2[CONV[31 & seq[j+2]]] ) { 
							#pragma omp atomic
							++Nibs[final_pfx]; // Increment the count for this prefix
							n_prefices++; // 
						} 
					} 
				}
			}
		}
		printf("There were %u records here of %lu prefices (%f s)\n",ns,n_prefices, omp_get_wtime()-wtime);
		Offsets = realloc(Offsets,ns*sizeof(*Offsets)); // Shrink the Offsets array to the correct size
		qsort(Offsets,ns,sizeof(*Offsets),cmpfunc);
		
		// Create the data structures
		uint64_t totWords = 0;
		for (uint64_t i = 0; i < nbins; ++i)
			totWords += Nibs[i];
		printf("In total, we need a structure that is %lu large.\n",totWords);
		printf("Now allocating %f GB of RAM...\n",(double)totWords*sizeof(KPod)/1073741824);
		wtime = omp_get_wtime();
		
		uint32_t *Lengths = malloc(ns*sizeof(*Lengths)); 
		KPod *KGrid = huge_malloc(totWords*sizeof(*KGrid));
		uint64_t *KIx = huge_calloc((nbins+1)*sizeof(*KIx));
		if (!KIx) {puts("OOM:KIx"); exit(3);}
		for (uint64_t i = 1; i < nbins; ++i) 
			KIx[i] = KIx[i-1] + Nibs[i-1];
		for (uint64_t i = 0; i < nbins; ++i) Nibs[i] = KIx[i];
		
		// Remember: PL is the prefix length, SL is the suffix length, and K is the k-mer length (PL+SL). 
		//uint32_t shifty = 32 - 2*PL; // shift to zero out extra letters in the prefix left over due to the type being larger than the prefix
		//kmer_t shifty2 = sizeof(kmer_t)*8 - 2*SL; // shift to zero out extra letters in the suffix left over due to the type being larger than the suffix
		//uint32_t rc_rshift_pfx = 2*(PL-1); // shift to put letter in first position of prefix
		//kmer_t rc_rshift_sfx = 2*(SL-1); // shift to put letter in first position of suffix

		/// KMER BIN FILLING LOOP. Read a fasta record, decompose in to k-mers, and add min(kmer, RC of kmer) to the KGrid.
		n_prefices = 0;
		#pragma omp parallel for schedule(dynamic) reduction(+:n_prefices)
		for (uint32_t i = 0; i < ns; ++i) {
			// Offsets stores indices into the whole Raw file for where each sequence starts
			char *seq = Raw + Offsets[i];
			uint64_t y = 0; while (seq[y] && seq[y] != '\n') ++y; // y is the last position in the sequence (inclusive!)
			int64_t lastAmbig = -1, tallyMade = 0, tallyPassed = 0;
			uint32_t pfx = 0; kmer_t sfx = 0; // stores binarized prefix and suffix of the k-mer
			uint32_t pfxRC = 0; kmer_t sfxRC = 0; // stores binarized prefix and suffix of the reverse complement of the k-mer
			for (int64_t j = 0; j < y; ++j) { // Start at the prefix length, and go to the end of the sequence
				uint32_t c = 0;  // c is the character we're looking at
				if (j >= SL) c = CONV[31 & seq[j-SL]]; // Get binary conversion of the letter, lag by suffix len (we don't want to overlap it)
					if (c==AMBIG) lastAmbig = j-SL, c = 0; // Mark last ambiguous position
				pfx = pfx << 2 | c; // Shift the prefix left by 2 bits and tack the new character onto the end
				pfx = pfx << shifty >> shifty; // remove the leading bits that are not part of the prefix
					
				// Now get the first character of the suffix, which is the current character
				kmer_t k = CONV[31 & seq[j]]; // the read head is always aligned with the end of the k-mer (also the end of the suffix)
				if (k==AMBIG) lastAmbig = lastAmbig > j ? lastAmbig : j, k = 0; // Mark last ambiguous position
				sfx = sfx << 2 | k; // Shift the suffix left by 2 bits and tack the new character onto the end
				sfx = sfx << shifty2 >> shifty2; // remove the leading bits that are not part of the suffix
				
				// Now get the next character for the reverse complement's prefix and suffix
				c = RCONV[31 & seq[j]]; // the read head is always aligned with the end of the k-mer (also the end of the suffix)
				c = c==AMBIG ? 0 : c; // Ambiguity already marked, so just cast to 0 ('A') to allow it to keep going
				pfxRC = pfxRC >> 2 | (uint32_t)c << rc_rshift_pfx; // Shift the prefix right by 2 bits and tack the new character onto the beginning
				if (j >= PL) k = RCONV[31 & seq[j-PL]]; // the read head is always aligned with the end of the k-mer (also the end of the suffix)
				k = k==AMBIG ? 0 : k; // Ambiguity already marked, so just cast to 0 ('A') to allow it to keep going
				sfxRC = sfxRC >> 2 | (kmer_t)k << rc_rshift_sfx; // Shift the suffix right by 2 bits and tack the new character onto the beginning

				// DEBUG: Print the prefix and suffix
				char pfxstr[PL+1]; pfxstr[PL] = 0;
				char pfxstrRC[PL+1]; pfxstrRC[PL] = 0;
				char sfxstr[SL+1]; sfxstr[SL] = 0;
				char sfxstrRC[SL+1]; sfxstrRC[SL] = 0;
				//for (int i = PL-1; i >= 0; --i) pfxstr[i] = "ACGT"[3 & (pfx >> 2*i)];
				/* num_to_dna(pfx,PL,pfxstr);
				num_to_dna(pfxRC,PL,pfxstrRC);
				num_to_dna(sfx,SL,sfxstr);
				num_to_dna(sfxRC,SL,sfxstrRC);
				char hasAmbig = j < lastAmbig + (int64_t)K ? 'Y' : 'N';
				printf("Ref %u, pos %d (%c). pfx sfx ; pfxRC sfxRC: %s %s ; %s %s\n",i,(int)j-K+1,hasAmbig,pfxstr,sfxstr,pfxstrRC,sfxstrRC); */

				// If the read head is past the last ambiguous position, our current k-mer is complete, so look it up
				if (j >= lastAmbig + (int64_t)K) { 
					// Determine which pairs (pfx,sfx or pfxRC,sfxRC) to keep by comparing the two. We simply keep the smaller of the two.
					uint32_t final_pfx = pfx; kmer_t final_sfx = sfx;
					if (pfxRC < pfx || pfxRC == pfx && sfxRC < sfx) final_pfx = pfxRC, final_sfx = sfxRC;
					// Now we have the final k-mer, so we can look it up in the KGrid
					// our choice of filter is simple: Are the characters before the prefix A then not A (or level 2: C)? If so, keep.
					if (!comp || j > K && FltF1[CONV[31 & seq[j-K]]] && FltF2[CONV[31 & seq[j-K-1]]] || 
								j < y - 1 && FltR1[CONV[31 & seq[j+1]]] && FltR2[CONV[31 & seq[j+2]]] ) { 
						if (comp) ++tallyPassed;
						++n_prefices;
						uint64_t pod_ix;
						#pragma omp atomic capture
						pod_ix = KIx[final_pfx]++;
						
						KGrid[pod_ix] = (KPod){final_sfx,i};
					} 
					++tallyMade;
					
				} 
				Lengths[i] = tallyMade + K - 1; 
			}
			//if (comp) printf("Query %d: %d kmers made, %d passed filter\n",i,tallyMade,tallyPassed);
		}


		printf("Time to catalogue all %lu prefices: %f\n",n_prefices, omp_get_wtime()-wtime);
		uint64_t numK = 0;
		#pragma omp parallel for schedule(dynamic,1024) reduction(+:numK)
		for (uint64_t i = 0; i < nbins; ++i) {
			if (KIx[i] <= Nibs[i]) continue;
			KPod *start = KGrid + Nibs[i];
			size_t num = KIx[i]-Nibs[i];
			qsort(start,num,sizeof(*start),kpackcmp);
			//char pfxstr[PL+1]; pfxstr[PL] = 0;
			//num_to_dna(i,PL,pfxstr); printf("Prefix %lu (%s): %lu kmers\n",i,pfxstr, num);
			numK += num;
		}
		printf("There were %lu k-mers.\n",numK);
		
		printf("Some stats on distributions!\n");
		// New: check for duplicates and multi-mapping k-mers...
		// Now including a step to reverse complement the k-mers in duplicate assessment
		uint64_t n_dupe = 0, n_multi = 0, n_total = 0;
		#pragma omp parallel for schedule(dynamic,1024) reduction(+:n_dupe,n_multi,n_total)
		for (uint64_t i = 0; i < nbins; ++i) {
			if (KIx[i] <= Nibs[i]) continue;
			n_total += KIx[i]-Nibs[i];
			// For later, we will convert this nib (which is a k-mer prefix!) into part of a k-mer which we can reverse complement later
			//char kmer_str[64] = {0}, rc_str[64] = {0};
			//num_to_dna(i,PL,kmer_str); // precompute the prefix shared by all k-mers in this bin
			//for (int j = 0; j < PL; ++j) rc_str[j+SL] = RCONS[31 & kmer_str[PL-1-j]]; // reverse complement the prefix
			//printf("* Prefix: %s; range to search for dupes: %lu <- %lu\n",kmer_str, Nibs[i]+1, KIx[i]);

			// First, check for exact duplicates
			for (uint64_t j = Nibs[i]+1; j < KIx[i]; ++j) {
				if (KGrid[j].sfx==KGrid[j-1].sfx) {
					++n_multi;
					if (KGrid[j].rix==KGrid[j-1].rix) ++n_dupe;
				}
			}
		}
		printf("Exact duplicates: %lu; across refs: %lu (%lu total)\n",n_dupe,n_multi, n_total);
		
		// Now write the results
		/* File structure V3:
		   1. Version byte [1] | compression level// no more rix_t size
		   2. Size of prefix
		   3. Size of suffix [1]
		   4. Size of kmer_t [1]
		   5. Num refs [4]
		   6. Num kmers [8]
		   7. All prefix indices [1 << (2 * #2) x 8]
		   8. All kmer data [(#2 + #4) * #6]
		   9. Number of adamantine kmer targets [8] // new
		   10. All adamantine kmer targets [(#2 + #4) * #9] // new
		   11. All reference lengths [#5 x 4] // new (v0.96)
		   12. Size of string data [8]
		   13. All strings [#9]
		   
		   14. Number of Ix1 in map [4] [0 means skip rest of file]
		   15. String size for Ix1 [8]
		   16. Ix1 strings dump [#12]
		   17. Number of Ix2 in map [4] [can be 0/skipped if no h2 map]
		   18. String size of Ix2 [8]
		   19. Ix2 strings dump [#15] 
		   //20. hpair_t dump [num ref by 8] // old
		   21. HPairs[0] dump [num ref by 4]
		   22. HPairs[1] dump [num ref by 4]
		*/
		
		uint64_t fileSz = 0, stringSz = 0;
		#pragma omp parallel for reduction(+:stringSz)
		for (uint32_t i = 0; i < ns; ++i) {
			uint64_t x = Offsets[i];
			uint64_t y = x; while (Raw[y] != '>') --y;
			stringSz += x-y - 1; // -1 for the '>' we're on, -1 '\n', but +1 '\0'
		}
		fileSz = 24 + sizeof(*Nibs)*(nbins+1) + sizeof(*KGrid)*numK + sizeof(*Lengths)*ns + stringSz;
		printf("Initial file size = %lu\n",fileSz+4);
		
		FILE *db = fopen(dbPath,"wb");
		if (!db) {puts("I/O error: invalid db output file"); exit(1);}
		setvbuf(db, 0, _IOFBF, 1<<22);
		wtime = omp_get_wtime();
		//uint8_t vno_derived = VNO + (comp); // comp can go up to 3, so the max value of vno_derived is VNO (currently 2) + 3 = 5.
		fputc((VNO << 4) | comp,db);  // 1: to extract the version number, do (vno_derived >> 4) & 0xF. to extract the compression level, do vno_derived & 0xF.
		fputc(PL,db); // 2
		fputc(SL,db); // 3
		fputc(sizeof(kmer_t),db); // 4
		size_t wrote = 4;
		wrote += fwrite(&ns,sizeof(ns),1,db); // 5
		wrote += fwrite(&numK,sizeof(numK),1,db); // 6
		uint64_t tally = 0;
		for (uint64_t i = 0; i < nbins+1; ++i) { // 7
			wrote += fwrite(&tally,sizeof(*Nibs),1,db);
			tally += KIx[i]-Nibs[i];
		}
		printf("First write: %f s.\n",omp_get_wtime()-wtime);
		wtime = omp_get_wtime();
		for (uint64_t i = 0; i < nbins; ++i) { // 8
			if (KIx[i] <= Nibs[i]) continue;
			uint64_t num = KIx[i]-Nibs[i];
			fwrite(KGrid + Nibs[i], sizeof(*KGrid), num,db);
		}
		printf("Second write: %f s.\n",omp_get_wtime()-wtime);
		
		// To determine the number of adamantine kmer targets, we need to go through the k-mers and count them
		
		
		// Compute (and potentially write) genome statistics, including adamantine k-mer targets
		uint64_t *Adamant_Prefix = 0; kmer_t *Adamant_Suffix = 0;
		uint64_t n_weak = 0;
		FILE *out = fopen(logPath,"wb"); // will be 0 if no log file specified
		if (out || adamantium) {
			if (adamantium) {
				Adamant_Prefix = huge_calloc((nbins+2)*sizeof(*Adamant_Prefix));
				Adamant_Prefix += 1; // Shift the pointer to the right by 1 so that we can use the 0th index as a sentinel
				// Prepopulate the two amantine k-mer bin prefix arrays to exactly mirror Nibs (and grow like KIx)
				for (uint64_t i = 0; i < nbins+1; ++i) Adamant_Prefix[i] = Nibs[i]; 
				Adamant_Suffix = malloc(numK*sizeof(*Adamant_Suffix)); // worst case is all k-mers are adamantine
			}
			uint32_t *TotK_m = huge_calloc((uint64_t)ns*sizeof(*TotK_m));
			uint32_t *TotUniq_m = huge_calloc((uint64_t)ns*sizeof(*TotUniq_m));
			uint32_t *TotUniq_dd_m = huge_calloc((uint64_t)ns*sizeof(*TotUniq_dd_m));
			uint32_t *TotWeak_m = huge_calloc((uint64_t)ns*sizeof(*TotWeak_m));
			#pragma omp parallel
			{
				int tid = omp_get_thread_num();
				uint32_t *TotK = TotK_m, *TotUniq = TotUniq_m, *TotWeak = TotWeak_m,
					*TotUniq_dd = TotUniq_dd_m;
				char kmer_str[128] = {0}, weakmer_str[128] = {0}; //, rc_str[128] = {0};
				
				if (tid) 
					TotWeak = huge_calloc((uint64_t)ns*sizeof(*TotWeak)),
					TotK = huge_calloc((uint64_t)ns*sizeof(*TotK)),
					TotUniq = huge_calloc((uint64_t)ns*sizeof(*TotUniq)),
					TotUniq_dd = huge_calloc((uint64_t)ns*sizeof(*TotUniq_dd));
				
				#pragma omp for schedule(dynamic,1024) reduction(+:n_weak)
				for (uint64_t i = 0; i < nbins; ++i) {
					if (KIx[i] <= Nibs[i]) continue; // empty
					uint32_t ambig = 0;
					uint64_t end = KIx[i], nd;
					kmer_t thisK = KGrid[Nibs[i]].sfx+1;
					
					for (uint64_t j = Nibs[i]; j < end; j += nd) {
						rix_t rix = KGrid[j].rix;
						
						// If new k-mer, check if ambig, store max value
						if (KGrid[j].sfx != thisK) {
							thisK = KGrid[j].sfx; 
							ambig = 0; 
							for (uint64_t k = j+1; k < end && KGrid[k].sfx == thisK; ++k)
								ambig |= KGrid[k].rix ^ rix; 
						}

						// Find number of in-ref copies
						nd = 1;
						for (uint64_t k = j+1; k < end && 
						KGrid[k].sfx == thisK && KGrid[k].rix == rix; ++k) ++nd;
						
						// Increment the appropriate variables
						if (!ambig) TotUniq[rix] += nd, ++TotUniq_dd[rix];
						if (!ambig && adamantium) {
							//printf("Kmer %s found in ref %u %u times\n",kmer_str,rix,nd);
							// Let's stress test the uniqueness of these supposedly "unambiguous" k-mers. We do this
							// by introducing every possible character mutation to the k-mer and seeing if it's still unique.
							// We can start by going through every position in the prefix and flipping the letter to each of the other 3.
							// Let's start with just the prefix, keeping the suffix constant. Let's directly modify prefix_mut.
							// Let's ignore the suffix entirely for now. Once we do all its mutations, we move onto suffix.
							int isWeak = 0; // is this k-mer weakly unique? (i.e. unique in this ref, but not strongly unique)
							uint32_t prefix_mut = i; kmer_t suffix_mut = thisK;
							for (int p = 0; p < PL; ++p) {
								//uint32_t prefix_mut = i; // the prefix of the k-mer, reset to the bin index
								for (uint32_t c = 0; c < 4; ++c) {
									prefix_mut = i; // the prefix of the k-mer, reset to the bin index
									if (c == (prefix_mut >> (2*(PL-1-p)) & 3)) continue; 
									prefix_mut = prefix_mut & ~(3 << (2*(PL-1-p))) | c << (2*(PL-1-p)); // set the p-th letter to c
									// Now that the prefix has mutated, we can look up the k-mer in the KGrid
									uint64_t ix_mut = Nibs[prefix_mut]; // the index of the first k-mer with this prefix
									uint64_t end_mut = KIx[prefix_mut]; // the index of the first k-mer with a different prefix

									// for debugging and printing, show old and new prefix
									/* num_to_dna(i,PL,kmer_str);
									num_to_dna(prefix_mut,PL,weakmer_str);
									num_to_dna(thisK,SL,weakmer_str+PL);
									kmer_str[PL] = 0; //weakmer_str[PL] = 0;
									printf("Prefix %s mutated to %s (range of matches %lu; %s)\n",kmer_str,weakmer_str,end_mut-ix_mut, weakmer_str);
 									*/
									// binary search for the current (original) suffix in the range of k-mers with the mutated prefix
									//static inline uint64_t LBS_k(KPod *KP, uint64_t N, kmer_t k)
									if (end_mut <= ix_mut) continue; // if there are no k-mers with this prefix, then it's unique
									uint64_t ix = LBS_k(KGrid+ix_mut, end_mut-ix_mut, thisK);

									if (ix != (uint64_t)-1) {
										//printf("--> Match found at %lu\n",ix);
										//for (uint64_t k = j+1; k < end && 
										//	KGrid[k].sfx == thisK && KGrid[k].rix == rix; ++k) ++nd;
										int difRef = 0; 
										for (uint64_t k = ix_mut + ix; k < end_mut && KGrid[k].sfx == thisK; ++k) 
											if (KGrid[k].rix != rix) { difRef = 1; break; }
										if (difRef) 
										{
											//num_to_dna(prefix_mut,PL,weakmer_str); num_to_dna(thisK,SL,weakmer_str+PL);
											isWeak = 1; break;
										}
									}
								}
								if (isWeak) break;
							}
							if (!isWeak) { // check the suffix now. We can directly modify suffix_mut. (keep old prefix, i)
								for (int p = 0; p < SL; ++p) {
									for (kmer_t c = 0; c < 4; ++c) {
										suffix_mut = thisK; // the original suffix of the k-mer
										if (c == (suffix_mut >> (2*(SL-1-p)) & 3)) continue; 
										suffix_mut = suffix_mut & ~(3 << (2*(SL-1-p))) | (kmer_t)c << (2*(SL-1-p)); // set the p-th letter to c
										// Now that the suffix has mutated, we can look up the k-mer in the KGrid
										uint64_t ix_mut = Nibs[i]; // the index of the first k-mer with this prefix
										uint64_t end_mut = KIx[i]; // the index of the first k-mer with a different prefix

										// for debugging and printing, show old and new suffix
										/* num_to_dna(thisK,SL,kmer_str+PL);
										num_to_dna(suffix_mut,SL,weakmer_str+PL);
										printf("Suffix %s mutated to %s (range of matches %lu; %s)\n",kmer_str+PL,weakmer_str+PL,end_mut-ix_mut,weakmer_str);
 										*/
										// binary search for the current (original) suffix in the range of k-mers with the mutated prefix
										//if (ix_mut == end_mut) continue; // if there are no k-mers with this prefix, then it's unique
										uint64_t ix = LBS_k(KGrid+ix_mut, end_mut-ix_mut, suffix_mut);
										// Test: see whether any reference found along the path is actually different from rix
										
										if (ix != (uint64_t)-1) { // if we found the k-mer (in a different ref?), then it's not unique
											//printf("--> Match found at %lu\n",ix);
											int difRef = 0; 
											for (uint64_t k = ix_mut+ix; k < end_mut && KGrid[k].sfx == suffix_mut; ++k) 
												if (KGrid[k].rix != rix) { difRef = 1; break; }
											if (difRef) 
											{
												//num_to_dna(i,PL,weakmer_str); num_to_dna(suffix_mut,SL,weakmer_str+PL);
												isWeak = 1; break;
											}
										}
									}
									if (isWeak) break;
								}
							}	
							//num_to_dna(i,PL,kmer_str); // the prefix shared by all k-mers in this bin
							//num_to_dna(thisK,SL,kmer_str+PL);  // append the suffix to the prefix
							//if (!isWeak) printf("%s\tSTRONG\t%u\n",kmer_str,rix);
							//else printf("%s\tWEAK\t%u\t%s\n",kmer_str,rix,weakmer_str);
							
							if (isWeak) {
								++TotWeak[rix];
								uint64_t adix;
								#pragma omp atomic capture
								adix = Adamant_Prefix[i]++;
								Adamant_Suffix[adix] = thisK;
								++n_weak; 
							}
						}
						TotK[rix] += nd;
					}
				}
				#pragma omp critical
				if (tid) for (uint32_t i = 0; i < ns; ++i) 
						TotK_m[i] += TotK[i], TotUniq_m[i] += TotUniq[i], TotUniq_dd_m[i] += TotUniq_dd[i], TotWeak_m[i] += TotWeak[i];
			}
			if (adamantium) {
				printf("There were %lu pre-adamantine k-mers.\n",n_weak);
				// Now go through and sort the suffixes for each prefix in Adamant_Prefix
				uint64_t unique_adamant = 0; 
				#pragma omp parallel for schedule(dynamic,1024) reduction(+:unique_adamant)
				for (uint64_t i = 0; i < nbins; ++i) {
					if (Adamant_Prefix[i] == Nibs[i]) continue; // empty bin
					uint64_t num = Adamant_Prefix[i]-Nibs[i];
					qsort(Adamant_Suffix+Nibs[i],num,sizeof(*Adamant_Suffix),kcmp);

					// Now we need to remove duplicates from the suffixes. We can do this by going through the sorted suffixes and
					// keeping track of the last suffix we saw. If the current suffix is the same as the last one, then we can
					// remove it. Otherwise, we keep it and update the last suffix.
					uint64_t last = Nibs[i];
					for (uint64_t j = Nibs[i]+1; j < Adamant_Prefix[i]; ++j) {
						if (Adamant_Suffix[j] == Adamant_Suffix[last]) continue;
						Adamant_Suffix[++last] = Adamant_Suffix[j];
					}
					unique_adamant += last - Nibs[i] + 1; // update the total number of unique suffixes
					Adamant_Prefix[i] = last+1; // update the prefix to point to the end of the unique suffixes
				}
				printf("There were %lu unique adamantine k-mers.\n",unique_adamant);

				// Now output them to the file using the new sorted deduplicated suffixes
				fwrite(&unique_adamant,sizeof(unique_adamant),1,db); // 9
				uint64_t counter = 0;
				for (uint64_t i = 0; i < nbins; ++i) { // 10
					uint64_t num = Adamant_Prefix[i]-Nibs[i];
					Adamant_Prefix[i] = Adamant_Prefix[i-1] + num; 
					if (num == 0) continue;
					counter += fwrite(Adamant_Suffix+Nibs[i],sizeof(*Adamant_Suffix),num,db); // 10
					// Rebuild from left to right as the starting index of the current bin
					
				}
				printf("Wrote %lu adamantine k-mers of size %lu\n",counter, sizeof(*Adamant_Suffix));
				Adamant_Prefix[nbins] = Adamant_Prefix[nbins-1]; // the last bin is the same as the second-to-last bin
				Adamant_Prefix -= 1; // Shift the pointer back to the left by 1
				fwrite(Adamant_Prefix,sizeof(*Adamant_Prefix),nbins+1,db); // 10 (part 2)
			} 
			if (out) {
				if (adamantium) fprintf(out,"Reference\tUsableLen\tTotalKmers\tUniqKmers\tUniqKmers_SC\tWeakKmers\n");
				else
				fprintf(out,"Reference\tUsableLen\tTotalKmers\tUniqKmers\tUniqKmers_SC\n");
				for (uint64_t i = 0; i < ns; ++i) {
					uint64_t x = Offsets[i];
					uint64_t y = x; while (Raw[y] != '>') --y;
					fwrite(Raw+y+1,1,x-y-2,out); // write the ref name
					if (adamantium) fprintf(out,"\t%u\t%u\t%u\t%u\t%u\n",Lengths[i], TotK_m[i], TotUniq_m[i], TotUniq_dd_m[i], TotWeak_m[i]);
					else fprintf(out,"\t%u\t%u\t%u\t%u\n",Lengths[i], TotK_m[i], TotUniq_m[i], TotUniq_dd_m[i]);
				}
			}

		}
		else if (!out) printf("No log file specified and no adamantium; won't produce tally\n");
		if (!adamantium) fwrite(&n_weak,sizeof(n_weak),1,db); // 9

		// Resume writing with #11
		wtime = omp_get_wtime();
		fwrite(Lengths,sizeof(*Lengths),ns,db); // 11 // new
		fwrite(&stringSz,sizeof(stringSz),1,db); // 12 
		for (uint32_t i = 0; i < ns; ++i) { // 13
			uint64_t x = Offsets[i];
			uint64_t y = x; while (Raw[y] != '>') --y;
			fwrite(Raw+y+1,1,x-y-2,db);fputc(0,db);
		}
		printf("Third write: %f s.\n",omp_get_wtime()-wtime);
		
		// Handle the H1/H2 mappings
		if (!ixPath) { // 14 (If no ixMap was provided, finish up.)
			uint32_t zeroRef = 0;
			fwrite(&zeroRef,sizeof(zeroRef),1,db);
			exit(0); 
		}
		
		wtime = omp_get_wtime();
		
		// Read the mapping file in. Up to 3 columns are used (first is index)
		// proof of concept -- read whole file, sort all 3 (perhaps parallel) by ptrs
		// Then perform classical deduplication and data structure creation 
		free(Nibs); free(KGrid); free(KIx); 
		
		// New observation: can determine a priori number and placement of bins 
		// for all possible interpolations, given sorted map. Also search/store!
		
		//gzFile map = gzopen(ixPath,"rb");
		//int mapFsz = gzread(map,wholeMap,(uint64_t)1 << 38);
		// Resume writing by calculating and writing at #14
		FILE *map = fopen(ixPath,"rb");
		if (!map) {fprintf(stderr,"Can't open map: %s\n",ixPath); exit(2);}
		
		uint64_t sz = 0; 
		fseeko(map,0,SEEK_END); sz = ftello(map); rewind(map);
		if (sz < 2) {fprintf(stderr,"ERR: map malformatted\n"); exit(2);}
		char *wholeMap = huge_calloc(sz+16); 
		size_t mapFsz = fread(wholeMap,1,sz,map);
		if (mapFsz != sz) {puts("BAD MAP!"); exit(10101);}
		
		uint64_t nL = 0;
		if (wholeMap[sz-1]!='\n') ++nL;
		
		#pragma omp parallel for reduction(+:nL)
		for (uint64_t i = 0; i < sz; ++i) 
			nL += wholeMap[i] == '\n';
		
		printf("Map contained %lu lines.\n",nL);
		//typedef struct {char *str, *h1, *h2;} map_str_t;
		//map_str_t *MapStrs = calloc(nL,sizeof(*MapStrs));
		char **RefStr = calloc(nL,sizeof(*RefStr)),
			**H1Str = calloc(nL,sizeof(*H1Str)),
			**H2Str = calloc(nL,sizeof(*H2Str));
		char *map_ptr = wholeMap;
		int ncol = 2; uint64_t numTimesH2showed = 0;
		for (uint64_t i = 0; i < nL; ++i) {
			RefStr[i] = map_ptr;
			while (*map_ptr != '\t' && *map_ptr != '\n') ++map_ptr;
			if (*map_ptr == '\n') {puts("Bad map! Need >1 columns!"); exit(2);}
			*map_ptr++ = 0;
			H1Str[i] = map_ptr;
			while (*map_ptr != '\t' && *map_ptr != '\n') ++map_ptr;
			if (*map_ptr == '\n') {*map_ptr++ = 0; ncol = 1; /* printf("1col line %lu\n",i); */ continue;}
			*map_ptr++ = 0;
			H2Str[i] = map_ptr;
			++numTimesH2showed;
			while (*map_ptr != '\n') ++map_ptr;
			*map_ptr++ = 0;
		}
		wholeMap[sz-1] = 0; // in case it's a newline
		printf("Detected %d columns (a second showed up %lu times). Parsed. [%f]\n",ncol,numTimesH2showed,omp_get_wtime()-wtime);
		
		wtime = omp_get_wtime();
		uint64_t nuniq_ref = 0, nuniq_h1 = 0, nuniq_h2 = 0;
		#pragma omp parallel sections
		{
			#pragma omp section
			{
				double wtimeL = omp_get_wtime();
				qsort(RefStr,nL,sizeof(*RefStr),strcmp_ptr);
				for (uint64_t i = 1; i < nL; ++i) // dedupe!
					if (strcmp(RefStr[i],RefStr[i-1])) 
						RefStr[nuniq_ref++] = RefStr[i-1];
				RefStr[nuniq_ref++] = RefStr[nL-1];
				printf("Name sort complete [%f]\n",omp_get_wtime()-wtimeL);
			}
			#pragma omp section
			{
				double wtimeL = omp_get_wtime();
				qsort(H1Str,nL,sizeof(*H1Str),strcmp_ptr);
				for (uint64_t i = 1; i < nL; ++i) // dedupe!
					if (strcmp(H1Str[i],H1Str[i-1])) 
						H1Str[nuniq_h1++] = H1Str[i-1];
				H1Str[nuniq_h1++] = H1Str[nL-1];
				printf("H1 sort complete [%f]\n",omp_get_wtime()-wtimeL);
			}
			#pragma omp section
			{
				double wtimeL = omp_get_wtime();
				if (ncol > 1) {
					qsort(H2Str,nL,sizeof(*H2Str),strcmp_ptr);
					for (uint64_t i = 1; i < nL; ++i) // dedupe!
						if (strcmp(H2Str[i],H2Str[i-1])) 
							H2Str[nuniq_h2++] = H2Str[i-1];
					H2Str[nuniq_h2++] = H2Str[nL-1];
				}
				printf("H2 sort complete [%f]\n",omp_get_wtime()-wtimeL);
			}
		}
		
		printf("All sorting complete. [%f]\n",omp_get_wtime()-wtime);
		printf("Unique: %lu refs, %lu H1's, %lu H2's.\n",
			nuniq_ref, nuniq_h1, nuniq_h2);
		wtime = omp_get_wtime();
		
		// Go thru each ref and match it to the list. 
		// Then lookup the h1 and h2 strings by proxy (past null).
		// Get the ids of the h1 and h2's from their respective lists.
		// Then add the pair [h1,h2] to the internal-ref-length struct.
		//hpair_t *HPairs = malloc(sizeof(*HPairs)*ns);
		uint32_t *HPairs[2] = { huge_calloc((uint64_t)ns*sizeof(*HPairs[0])),
			huge_calloc((uint64_t)ns*sizeof(*HPairs[1])) };
		#pragma omp parallel for schedule(dynamic,16)
		for (uint32_t i = 0; i < ns; ++i) {
			char *ref = Raw+Offsets[i];
			while (*ref != '>') --ref;
			uint64_t refmatch = binsearch_str_d(RefStr, ++ref, nuniq_ref);
			if (refmatch == (uint64_t)-1) {
				fprintf(stderr,"ERR: Map missing '%.*s'\n",
					(int)(Raw+Offsets[i] - ref - 1), ref); exit(2);
			}
			char *h = RefStr[refmatch];
			uint64_t h1match = 0, h2match = 0;
			while (*h) ++h;
			h1match = binsearch_str(H1Str,++h,nuniq_h1);
			if (h1match == (uint64_t)-1) {puts("INTERNAL ERROR H1"); exit(9);}
			if (ncol > 1) {
				while (*h) ++h;
				h2match = binsearch_str(H2Str,++h,nuniq_h2);
				if (h2match == (uint64_t)-1) {puts("INTERNAL ERROR H2"); exit(9);}
			}
			//HPairs[i] = (hpair_t){h1match,h2match};
			HPairs[0][i] = h1match, HPairs[1][i] = h2match;
		}
		
		// Should we make the taxonomy sublevel-map now or later? 
		// (In theory it can be made during read-in but...)
		// Update: now made during read-in!
		
		// Need to dump:
		// 1. String data (for both headers only; tracking length)
		// 2. The hpair_t dump -- update: now the two separate H arrays
		
		/*
		   12. Number of Ix1 in map [4] [0 means skip rest of file]
		   13. String size for Ix1 [8]
		   14. Ix1 strings dump [#12]
		   15. Number of Ix2 in map [4] [can be 0/skipped if no h2 map]
		   16. String size of Ix2 [8]
		   17. Ix2 strings dump [#15] 
		   //18. hpair_t dump [num ref by 8] // old
		   18. HPairs[0] dump [num ref by 4]
		   19. HPairs[1] dump [num ref by 4]
		*/
		
		stringSz = 0;
		#pragma omp parallel for reduction(+:stringSz)
		for (uint64_t i = 0; i < nuniq_h1; ++i) 
			stringSz += strlen(H1Str[i]) + 1;
		
		fileSz += 4 + 8 + stringSz;
		fwrite(&nuniq_h1,4,1,db); // # 11
		fwrite(&stringSz,sizeof(stringSz),1,db); // # 12
		for (uint64_t i = 0; i < nuniq_h1; ++i) 
			fwrite(H1Str[i],1,strlen(H1Str[i])+1,db); // #13
		
		fileSz += 4; // H2
		fwrite(&nuniq_h2,4,1,db); // #14
		stringSz = 0;
		if (ncol > 1) {
			#pragma omp parallel for reduction(+:stringSz)
			for (uint64_t i = 0; i < nuniq_h2; ++i) 
				stringSz += strlen(H2Str[i]) + 1;
			
			fwrite(&stringSz,sizeof(stringSz),1,db); // #15
			for (uint64_t i = 0; i < nuniq_h2; ++i) 
				fwrite(H2Str[i],1,strlen(H2Str[i])+1,db); // #16
		} else fwrite(&nuniq_h2,8,1,db); // 15 & 16 w/no col2
		fileSz += 8 + stringSz;
		
		fileSz += (uint64_t)ns*sizeof(*HPairs[0]);
		fwrite(HPairs[0],sizeof(*HPairs[0]),ns,db); // #17
		if (nuniq_h2) fileSz += (uint64_t)ns*sizeof(*HPairs[1]),
			fwrite(HPairs[1],sizeof(*HPairs[1]),ns,db); // #18
		
		printf("Final filesize: %lu\n",fileSz);
		
		exit(0);
	}
	
	/// Now handle the parsing and searching. 
	// TODO: make parser size-aware, not "num sequences" fixed
	
	if (ixPath) puts("WARNING: map file only applicable during DB BUILD");
	//  Read the database in.
	FILE *db = fopen(dbPath,"rb");
	if (!db) {puts("ERROR: bad input"); exit(2);}
	struct stat sb; int fno = fileno(db); fstat(fno,&sb);
	uint64_t fsz = sb.st_size, place = 0;
	char *Raw = mmap(0, fsz+16, PROT_READ, MAP_SHARED | 
		MAP_POPULATE, fno, 0);
	madvise(Raw,fsz,MADV_WILLNEED);
	double wtime = omp_get_wtime();
	if (doCopyMem) {
		puts("Copying database into local memory...");
		char *copied = huge_malloc(fsz+16);
		if (!copied) {puts("Can't do copymem on this system. Exiting..."); exit(3);}
		memcpy(copied,Raw,fsz+1);
		munmap(Raw,fsz+16);
		Raw = copied;
		printf("Copied into local memory [%f]\n",omp_get_wtime()-wtime);
		wtime = omp_get_wtime();
	}
	uint32_t ver = Raw[0] >> 4, rixSz = sizeof(rix_t), compdb = (uint8_t)Raw[0] & 15,
		PL = Raw[1], SL = Raw[2], ktSz = Raw[3];
	uint32_t numRef = *(uint32_t *)(Raw+4);
	uint64_t numK = *(uint64_t *)(Raw+8);
	// Raw[0] was stored as: fputc((VNO << 4) | comp,db); // 1 byte
	if (ver != 5) {printf("ERROR: bad DB version %d (CrossTree < v2.0)\n",ver); exit(2);}
	if (comp < -1 || comp > 2) {printf("WARNING: bad compression level %d (range is 0-2)\n",comp); comp = -1;}
	if (comp != -1 && comp != compdb) 
		printf("WARNING: forcing query compression to %d instead of DB's native %d\n",comp,compdb);
	comp = comp == -1 ? compdb : comp; // -1 = auto, 0 = no, 1+ = yes
	if (comp == 2) FltF2 = FLTF2_2, FltR2 = FLTR2_2; // use stronger filters for compressive filtering

	printf("DBv: %d, DBcomp = %d, rixSz = %d, PL = %d, SL = %d, ktSz = %d\n", 
		VNO, compdb, rixSz, PL, SL, ktSz);
	printf("Number of refs = %u, kmers = %lu\n",numRef, numK);
	if (sizeof(kmer_t)!=ktSz || sizeof(rix_t) != rixSz) 
		{puts("ERROR: wrong K or R size(s) for this DB/xtree build"); exit(2);}
	place = 16; //#1-6
	
	// Initialize stats
	uint64_t K = PL+SL;
	uint32_t kpre_shf = PL*2-2;
	kmer_t kpst_shf = SL*2-2;
	uint32_t pre_bshf = 32-(PL*2), pre_bshf_2 = pre_bshf+2;
	kmer_t pst_bshf = sizeof(kmer_t)*8-(SL*2), pst_bshf_2 = pst_bshf+2;

	// Read the prefix array. 
	uint64_t *Nibs = (uint64_t *)(Raw+place);
	uint64_t nbins = (uint64_t)1 << (2*PL);
	place += (nbins+1) * sizeof(*Nibs); //#7
	//printf("DEBUG: nbins: %lu, Nibs[nbins] = %lu\n",nbins,Nibs[nbins]);
	KPod *KGrid = (KPod *)(Raw + place);
	printf("Size of kpod = %ld\n",sizeof(*KGrid));
	place += numK*sizeof(*KGrid); //#8

	// Read the adamantine prefix/suffix arrays (if present)
	uint64_t *Adamant_Prefix = 0; // Only stored one
	kmer_t *Adamant_Suffix = 0;
	uint64_t n_adamant = *(uint64_t *)(Raw + place);
	place += sizeof(n_adamant); //#9
	//printf("Adamant DB size = %lu (current position = %lu)\n",n_adamant,place); fflush(stdout);
	if (n_adamant) { // #10 Read adamantium DBs
		// Read the prefix and suffix arrays
		Adamant_Suffix = (kmer_t *)(Raw + place);
		place += n_adamant*sizeof(*Adamant_Suffix); //#10 (part 1)
		Adamant_Prefix = (uint64_t *)(Raw + place);
		place += (nbins+1)*sizeof(*Adamant_Prefix); //#10 (part 2)
		printf("Adamantly read adamant DB of %u adamantine elements of adamantium\n",n_adamant); 
		if (n_adamant > 0) printf("And it was super hard to do so...\n");
		else if (adamantium) puts("WARNING: Empty Adamantium DB; disabling Adamantium Engine"), adamantium = 0;
	}

	// Read the ref sizes
	uint32_t *RefSizes = (uint32_t *)(Raw + place);
	place += numRef*sizeof(*RefSizes); //#11
	
	// Read the ref names
	uint64_t stringSz = *(uint64_t *)(Raw + place);
	//printf("String size = %lu\n",stringSz); fflush(stdout);
	place += sizeof(stringSz); //#9
	char *RefRaw = Raw + place;
	place += stringSz; //#10
	char **RefNames = malloc(sizeof(*RefNames)*numRef); // defer
	//printf("String size = %lu\n",stringSz);
	//printf("String 1 = %s\n",RefRaw);
	
	// Read the h1 and h2 lists
	uint32_t nuniq_h1 = 0, nuniq_h2 = 0;
	char *H1Raw = 0, *H2Raw = 0, **HStr[2] = {0,0}; // **H1Str = 0, **H2Str = 0;
	//hpair_t *HPairs = 0;
	uint32_t *HPairs[2] = { 0,0 };
	nuniq_h1 = *(uint32_t *)(Raw+place);
	//printf("nuniq_h1 = %u\n",nuniq_h1);
	place += sizeof(nuniq_h1); //#11
	if (nuniq_h1) {
		HStr[0] = malloc((uint64_t)nuniq_h1*sizeof(*HStr[0]));
		stringSz = *(uint64_t *)(Raw + place);
		place += sizeof(stringSz); //#12
		H1Raw = Raw + place;
		place += stringSz; //#13
		
		nuniq_h2 = *(uint32_t *)(Raw+place);
		place += sizeof(nuniq_h2); //#14
		if (nuniq_h2) HStr[1] = malloc((uint64_t)nuniq_h2*sizeof(*HStr[1]));
		stringSz = *(uint64_t *)(Raw + place);
		place += sizeof(stringSz); //#15
		
		H2Raw = Raw + place;
		place += stringSz; //#16
		HPairs[0] = (uint32_t *)(Raw + place);
		place += (uint64_t)numRef * sizeof(*HPairs[0]); // #17
		if (nuniq_h2) HPairs[1] = (uint32_t *)(Raw + place);
		place += !nuniq_h2? 0 : (uint64_t)numRef * sizeof(*HPairs[1]); // #18
	}
	uint32_t NUniqH[2] = {nuniq_h1, nuniq_h2};
	printf("Read file of size = %lu (h1: %u, h2: %u)\n",place,nuniq_h1,nuniq_h2);
	
	// Make data structures to contain bin mappings
	// (first think how you want to count and create the interpolated
	// bins to be accessed -- ideally you'd specify an interpolation 
	// and get the unique (interpolated) index back.
	uint32_t **LBins[2] = {calloc(4096,sizeof(*LBins[0])), calloc(4096,sizeof(*LBins[1]))};
	//uint32_t **Tbins_h1 = calloc(4096,sizeof(*Tbins_h1));
	//uint32_t **Tbins_h2 = calloc(4096,sizeof(*Tbins_h2));
	
	// Parse the strings in parallel
	if (!nuniq_h1) 
		printf("WARNING: No taxonomy was included during DB formation\n");
	//FILE *debug = fopen("debug.txt","wb"); // DEBUG ONLY
	#pragma omp parallel sections
	{
		#pragma omp section
		{
		for (uint32_t i = 0; i < numRef; ++i) // Ref names
			RefNames[i] = RefRaw, RefRaw = strchr(RefRaw,0)+1;
		

		/* FILE *log = fopen("Dbg.txt","wb");
		for (uint32_t i = 0; i < numRef; ++i)
			fprintf(log,"%u\t%s\n",i,RefNames[i]);
		fflush(log); fclose(log); */
		//exit(10101);
		}

		#pragma omp section
		{
			//char **H1Str = HStr[0]; //uint32_t **Tbins_h1 = LBins[0];
			for (uint32_t i = 0; i < nuniq_h1; ++i) // H1 names
				HStr[0][i] = H1Raw, H1Raw = strchr(H1Raw,0)+1;
			for (uint32_t i = 0; i < nuniq_h1; ++i) {
				char *ref = HStr[0][i], *ptr = ref-1;
				int lv = 0;
				while (ptr = strchr(ptr+1,';')) {
					//fprintf(debug,"Searching: %s until %lu...\n --> %.*s\n",ref,ptr-ref,(int)(ptr-ref),ref);
					int64_t find = binsearch_str_L(HStr[0],ref,nuniq_h1,ptr-ref);
					//fprintf(debug,"%.*s\t%ld\n",(int)(ptr-ref),ref,find);
					if (!LBins[0][lv]) {
						LBins[0][lv] = calloc(nuniq_h1,sizeof(*LBins[0][lv]));
						for (uint32_t j = 0; j < nuniq_h1; ++j) LBins[0][lv][j]=-1;
					}
					LBins[0][lv++][i] = find;
				}
			}
		}
		#pragma omp section
		{
			for (uint32_t i = 0; i < nuniq_h2; ++i) // H2 names
				HStr[1][i] = H2Raw, H2Raw = strchr(H2Raw,0)+1;
			for (uint32_t i = 0; i < nuniq_h2; ++i) {
				char *ref = HStr[1][i], *ptr = ref-1;
				int lv = 0;
				while (ptr = strchr(ptr+1,';')) {
					int64_t find = binsearch_str_L(HStr[1],ref,nuniq_h2,ptr-ref);
					if (!LBins[1][lv]) {
						LBins[1][lv] = calloc(nuniq_h2,sizeof(*LBins[1][lv]));
						for (uint32_t j = 0; j < nuniq_h2; ++j) LBins[1][lv][j]=-1;
					}
					LBins[1][lv++][i] = find;
				}
			}
		}
	}

	printf("Database read in %f seconds\n",omp_get_wtime()-wtime);
	
	// debug: print everything out.
	/* FILE *debug = fopen("debug.txt","wb");
	for (uint32_t i = 0; i < numRef; ++i) 
		fprintf(debug,"%s\t%s\t%s\n",RefNames[i],H1Str[HPairs[i].h1],H2Str[HPairs[i].h2]);
	fprintf(debug,"AND NOW THE H1 ACTION BOYS:\n");
	for (uint32_t i = 0; i < nuniq_h1; ++i) 
		fprintf(debug,"%s\t%u\t%u\n",H1Str[i],Tbins_h1[i+1],Tbins_h1[i+1]-Tbins_h1[i]);
	fprintf(debug,"AND NOW THE H2 ACTION BOYS:\n");
	if (nuniq_h2) for (uint32_t i = 0; i < nuniq_h2; ++i) 
		fprintf(debug,"%s\t%u\t%u\n",H2Str[i],Tbins_h2[i+1],Tbins_h2[i+1]-Tbins_h2[i]);
	
	// Find all the combos in each grid!!
	fprintf(debug,"AND NOW THE SEARCH IS ON!\n");
	for (uint32_t i = 0; i < nuniq_h1; ++i) {
		char *ref = H1Str[i], *ptr = ref;
		while (ptr = strchr(ptr,';')) {
			//printf("Searching: %s until %lu...\n",ref,ptr-ref-1);
			int64_t find = binsearch_str_L(H1Str,ref,nuniq_h1,ptr-ref-1);
			int nsemi = 0; for (char *p = ptr; *p; ++p) nsemi += *p==';';
			fprintf(debug,"%.*s\t%ld\t%ld\n",(int)(ptr-ref),ref,find,Tbins_h1[find+1]-nsemi);
			//printf("%.*s\t%ld\n",(int)(ptr-ref-1),ref,find);
			++ptr;
		}
		//int64_t find = binsearch_str(H1Str,ref,nuniq_h1);
		int64_t find = binsearch_str_L(H1Str,ref,nuniq_h1,9999);
		int64_t find2 = binsearch_str(H1Str,ref,nuniq_h1);
		if (find2 != find) {puts("ERROR FIND"); exit(6);}
		fprintf(debug,"%s\t%ld\t%ld\n",ref,find,Tbins_h1[find+1]);
	}
	
	exit(1); */
	
	uint32_t *QueryAligns = calloc(numK,sizeof(*QueryAligns)),
		*FullQueryAligns = calloc((uint64_t)numRef,sizeof(*FullQueryAligns));
	
	
	// Open those queries up and align 'em
	KPod *EndPod = KGrid + Nibs[nbins] + 1;
	uint64_t n_raw = 0, n_filt = 0, n_matched = 0;
	
	wtime = omp_get_wtime();
	printf("\n* Beginning alignment...\n");
	
	gzFile in;
	if (!strcmp(seqPath,"-")) in = gzdopen(fileno(stdin),"rb");
	else in = gzopen(seqPath,"rb");
	if (!in) {puts("ERROR: bad input fast[a/q][.gz]"); exit(2);}
	uint8_t *lineO = calloc(16 + (1 << 28),1), *line = lineO+16, 
		*head = malloc(1 << 20);
	uint32_t qChunk = 1<<16;
	uint8_t **QBucket = malloc(sizeof(*QBucket)*qChunk);
	char **HBucket = malloc(sizeof(*HBucket)*qChunk);
	*QBucket = calloc((uint64_t)131072*qChunk,sizeof(**QBucket));
	*HBucket = calloc((uint64_t)131072*qChunk,sizeof(**HBucket));
	
	// Make a bucket for temporary LCA & ref votes, per query
	uint64_t maxQsz = 1 << 27; //27; // arbitrary? Set a limit?
	uint64_t masterBnSz = (uint64_t)1 << 31;
	uint64_t masterLstSz = (uint64_t)1 << 31;

	typedef struct {uint32_t p; uint64_t s;} SBin_t;
	
	SBin_t **SBins = malloc(sizeof(*SBins)*threads);   // For ixing matches
	int32_t **RBins = malloc(sizeof(*RBins)*threads),  // For tally stores
		**TBins = malloc(sizeof(*TBins)*threads);
	#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		//printf("I'm thread %d / %d\n",tid, threads);
		int counter = 0;
		counter += !!(SBins[tid] = huge_malloc(maxQsz*sizeof(*SBins)));
		counter += !!(TBins[tid] = huge_malloc(maxQsz*sizeof(*TBins)));
		counter += !!(RBins[tid] = huge_calloc((uint64_t)numRef*sizeof(*RBins)));
		if (counter != 3) printf("Error allocating bins on thread %d (counter = %d)\n",tid,counter);
	}
	// For the capitalist bins (variable-length bins depending on hits). Thread-local.
	uint64_t **C_ixs = huge_malloc(sizeof(*C_ixs)*threads); //[3];   // for rix, h1, h2 cap arrays
	uint64_t **C_szs = huge_malloc(sizeof(*C_szs)*threads); //[3];   // the current size of the bins
	uint64_t ***C_bins = huge_malloc(sizeof(*C_bins)*threads); //[3]; // the bin arrays themselves
	//uint64_t C_MAX = (uint64_t)1 << 21;
	uint64_t C_INIT = (uint64_t)1 << 12;
	#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		C_ixs[tid] = calloc(3,sizeof(**C_ixs));
		C_szs[tid] = calloc(3,sizeof(**C_szs)); // enables thread resizing
		C_bins[tid] = calloc(3,sizeof(**C_bins));
		for (int i = 0; i < 3; ++i) 
			C_bins[tid][i] = malloc(C_INIT*sizeof(***C_bins)),
			C_szs[tid][i] = C_INIT-1;
	}
	// For "master" bins -- number of bins equal to number of queries. Single global storage.
	//typedef union {uint32_t a[4]; __uint128_t n; uint8_t b[16];} MasterBin_t;
	MasterBin_t *MasterBin = calloc(masterBnSz,sizeof(*MasterBin)); // bin IX // 32
	if (!MasterBin) {printf("ERROR: failed to allocate master bin!\n"); exit(3);}
	// a[0] = rix, a[1] = h1, a[2] = h2, a[3] = ...?
	uint8_t  **MasterIx = calloc(3,sizeof(*MasterIx)); // which thread's bin
	uint64_t **MasterList = calloc(3,sizeof(*MasterList)); 
	for (int i = 0; i < 3; ++i) {

		int counter = 0;
		counter += !!(MasterList[i] = calloc(masterLstSz,sizeof(*MasterList[i]))); //31
		counter += !!(MasterIx[i] = calloc(masterLstSz,sizeof(*MasterIx[i]))); //31
		if (counter != 2) printf("Error allocating master list (counter = %d)\n",counter);
	}

	FILE *outq = 0;
	if (perqPath) {
		outq = fopen(perqPath,"wb");
		if (!outq) {puts("ERROR: can't open per-q output file!"); exit(2);}
	}
	
	uint32_t shifty = 32 - 2*PL; // shift to zero out extra letters in the prefix left over due to the type being larger than the prefix
	kmer_t shifty2 = sizeof(kmer_t)*8 - 2*SL; // shift to zero out extra letters in the suffix left over due to the type being larger than the suffix
	uint32_t rc_rshift_pfx = 2*(PL-1); // shift to put letter in first position of prefix
	kmer_t rc_rshift_sfx = 2*(SL-1); // shift to put letter in first position of suffix
	wtime = omp_get_wtime();
	uint64_t nq, NQ = 0, nAligns = 0;
	while (nq = get_queries(in,QBucket,HBucket,head,line,qChunk,INT32_MAX)) {
		printf("Processed %lu queries\r",NQ); fflush(stdout);
		if (nq + NQ >= INT32_MAX) {puts("Exceeded 2B queries; stopping"); nq = 0; break;}
		#pragma omp parallel
		{
			int tid = omp_get_thread_num();
			SBin_t *SBin = SBins[tid];
			int32_t *RBin = RBins[tid]; uint32_t *TBin = TBins[tid];
			uint64_t *C_ix = C_ixs[tid], *C_sz = C_szs[tid],
				**C_bin = C_bins[tid];
			#pragma omp for schedule(dynamic,1) reduction(+:n_raw,n_filt,n_matched)
			for (uint64_t q = 0; q < nq; ++q) {
				
				uint8_t *Inf = QBucket[q];
				char *Qhed = HBucket[q];
				uint64_t x = 0, y = strlen(Inf);
				if (y >= maxQsz) 
					y = maxQsz-1,
					printf("\nQuery %lu exceed maximum length of %lu; truncating\n",q+NQ,maxQsz-1);
				
				n_raw += y-K+1; // Increment the number of raw kmers by the number of kmers in the read
				
				char *seq = Inf; // Inf is the read
				int64_t lastAmbig = -1;
				uint32_t pfx = 0; kmer_t sfx = 0; // stores binarized prefix and suffix of the k-mer
				uint32_t pfxRC = 0; kmer_t sfxRC = 0; // stores binarized prefix and suffix of the reverse complement of the k-mer
				uint32_t tix = 0; // the index where the current k-mer is stored (from this read)
				for (int64_t j = 0; j < y; ++j) { // Start at the prefix length, and go to the end of the sequence
					uint64_t slen = 0, seed = -1;
					uint32_t c = 0;  // c is the character we're looking at
					if (j >= SL) c = CONV[31 & seq[j-SL]]; // Get binary conversion of the letter, lag by suffix len (we don't want to overlap it)
						if (c==AMBIG) lastAmbig = j-SL, c = 0; // Mark last ambiguous position
					pfx = pfx << 2 | c; // Shift the prefix left by 2 bits and tack the new character onto the end
					pfx = pfx << shifty >> shifty; // remove the leading bits that are not part of the prefix

					// Now get the first character of the suffix, which is the current character
					kmer_t k = CONV[31 & seq[j]]; // the read head is always aligned with the end of the k-mer (also the end of the suffix)
					if (k==AMBIG) lastAmbig = lastAmbig > j ? lastAmbig : j, k = 0; // Mark last ambiguous position
					sfx = sfx << 2 | k; // Shift the suffix left by 2 bits and tack the new character onto the end
					sfx = sfx << shifty2 >> shifty2; // remove the leading bits that are not part of the suffix

					// Now get the next character for the reverse complement's prefix and suffix
					c = RCONV[31 & seq[j]]; // the read head is always aligned with the end of the k-mer (also the end of the suffix)
					c = c==AMBIG ? 0 : c; // Ambiguity already marked, so just cast to 0 ('A') to allow it to keep going
					pfxRC = pfxRC >> 2 | (uint32_t)c << rc_rshift_pfx; // Shift the prefix right by 2 bits and tack the new character onto the beginning
					if (j >= PL) k = RCONV[31 & seq[j-PL]]; // the read head is always aligned with the end of the k-mer (also the end of the suffix)
					k = k==AMBIG ? 0 : k; // Ambiguity already marked, so just cast to 0 ('A') to allow it to keep going
					sfxRC = sfxRC >> 2 | (kmer_t)k << rc_rshift_sfx; // Shift the suffix right by 2 bits and tack the new character onto the beginning
					// If the read head is past the last ambiguous position, our current k-mer is complete, so look it up
					if (j >= lastAmbig + (int64_t)K) { 
						// process the kmer
						++n_filt;
						uint32_t final_pfx = pfx; kmer_t final_sfx = sfx;
						if (pfxRC < pfx || pfxRC == pfx && sfxRC < sfx) final_pfx = pfxRC, final_sfx = sfxRC;

						// Use the comp logic as in the previous code: CA_ or _TG
						if (!comp || j > K && FltF1[CONV[31 & seq[j-K]]] && FltF2[CONV[31 & seq[j-K-1]]] || 
								j < y - 1 && FltR1[CONV[31 & seq[j+1]]] && FltR2[CONV[31 & seq[j+2]]] ) { 
							// We now only work with final_pfx and final_sfx (the canonical k-mer)
							slen = Nibs[final_pfx+1]-Nibs[final_pfx];
							if (!slen) continue; // If the k-mer is not in the prefix index, skip it
							
							seed = LBS_k(KGrid+Nibs[final_pfx],slen,final_sfx);
							if (seed==(uint64_t)-1) continue; // If the k-mer is not in the suffix index, skip it
	
							SBin[tix++] = (SBin_t){final_pfx,seed+Nibs[final_pfx]};
							++n_matched;
						}
					}
				}
				
				// Resize the query tally bins if needed (rix, h1, h2)
				if (doRedist) for (int j = 0; j < 3; ++j) if (tix+3 >= C_sz[j]-C_ix[j]) 
					//{printf("OOM: [%lu] C_ix[%d]\n",NQ+q,j); exit(3);}
					C_bin[j] = C_bins[tid][j] = realloc(C_bin[j],(C_sz[j]=2*(C_sz[j]+tix+3))*sizeof(*C_bin[j]));
						//madvise(C_bin[j], C_sz[j], MADV_HUGEPAGE);
				
				if (doRedist) for (int j = 0; j < 3; ++j) 
					MasterList[j][NQ+q] = C_ix[j],
					MasterIx[j][NQ+q] = tid;
				MasterBin[q+NQ].n = (__uint128_t)-1; // initialize MasterBin (ref, h1, h2, ?)
				if (!tix) {
					if (perqPath) fprintf(outq,"%s\tNo matches found\n",Qhed);
					if (doRedist) for (int j = 0; j < 3; ++j) 
						C_bin[j][C_ix[j]++] = -1;
					continue; // no alignments, no worries
				}
				
				/// Per-query processing (reporting, taxonomy tally, etc)
				// Vote for which reference to store
				int tempix = 0;
				// Go through and tally refs into storage bin RBin (init'ed to 0)
				for (uint32_t i = 0; i < tix; ++i) {
					uint32_t pfx = SBin[i].p;
					uint64_t ix = SBin[i].s;
					
					uint64_t hardstop = Nibs[pfx+1];
					rix_t prev_rix = -1; // set to non-this for dupe detection
					kmer_t prev_sfx = KGrid[ix].sfx;
					for (uint64_t j = ix; j < hardstop; ++j) {
						if (KGrid[j].sfx != prev_sfx) break;
						rix_t rix = KGrid[j].rix;
						//printf("Q %s: [%u/%u, piece %lu], %u / %s\n",Qhed,i,tix,j-ix,rix,RefNames[rix]);
						if (rix == prev_rix) continue; // don't double-count refs
						if (!RBin[rix]) TBin[tempix++] = rix;
						++RBin[rix];
						prev_rix = rix;
					}
				} // All unique references for this query are now in TBin by name, and RBin by count.

				// Resize the query tally bins if needed (rix, h1, h2)
				if (doRedist) for (int j = 0; j < 3; ++j) if (tempix+3 >= C_sz[j]-C_ix[j]) 
					//{printf("OOM: [%lu] C_ix[%d]\n",NQ+q,j); exit(3);}
					C_bin[j] = C_bins[tid][j] = realloc(C_bin[j],(C_sz[j]=2*(C_sz[j]+tempix+3))*sizeof(*C_bin[j]));
						//madvise(C_bin[j], C_sz[j], MADV_HUGEPAGE);
				
				// Consider all references that matched, uniquely
				uint32_t max = 0, max2 = 0, sum = 0;
				rix_t maxRix = -1, maxRix2 = -1;
				for (int i = 0; i < tempix; ++i) { // tally to determine max and 2nd max ref
					uint32_t rix = TBin[i];
					if (RBin[rix] > max || RBin[rix]==max && rix < maxRix) 
						max2 = max, maxRix2 = maxRix,
						max = RBin[rix], maxRix = rix;
					else if (RBin[rix] > max2)
						max2 = RBin[rix], maxRix2 = rix;
				}

				// Early terminate if no ref, or one ref got N or more assignments (controled by cmdline):
				if (!tempix || maxRix == -1 || max < nUniqMatches) {
					if (perqPath) fprintf(outq,"%s\tNo matches found\n",Qhed);
					if (doRedist) for (int j = 0; j < 3; ++j) 
						C_bin[j][C_ix[j]++] = -1;
					for (int i = 0; i < tempix; ++i) RBin[TBin[i]] = 0;
					continue; // no alignments, no worries
				}
				
				if (covPath) for (uint32_t i = 0; i < tix; ++i) { // coverage of k-mers
					uint64_t ix = SBin[i].s;
					uint64_t hardstop = Nibs[SBin[i].p+1];
					kmer_t prev_sfx = KGrid[ix].sfx;
					for (uint64_t j = ix; j < hardstop; ++j) {
						if (KGrid[j].sfx != prev_sfx) break;
						rix_t rix = KGrid[j].rix;
						if (doForage || RBin[rix] == max || (doHalfForage && RBin[rix] >= 2*max/3 /* && RBin[rix] >= max2 */ )) { // if forage (meaning "grab all refs" for coverage table), or if this ref is tied for max
							#pragma omp atomic
							++QueryAligns[doFullLCA? ix : j]; // TODO: change "ix" to "j" ? ix is the "first unique match" found by the left binary search. 
							// we just made it toggleable, so we can have it both ways. To enable the old 'ix' behavior, specify --shallow-lca
						}
					}
				}
				for (int i = 0; i < tempix; ++i) { // ref-level set and clear
					uint32_t rix = TBin[i];
					if (RBin[rix] == max) {
						if (covPath) {
							#pragma omp atomic
							++FullQueryAligns[rix]; // # all tied refs
						}
						if (doRedist) C_bin[0][C_ix[0]++] = rix;
					} else if (covPath && (doForage || (doHalfForage && RBin[rix] >= 2*max/3 /* && RBin[rix] >= max2 */))) {
						#pragma omp atomic
						++FullQueryAligns[rix]; // # all refs
					}
					RBin[rix] = 0; // clear the bin so that reference counting doesn't double-count
				}
				if (doRedist) C_bin[0][C_ix[0]++] = -1;
				
				uint32_t finalRix = maxRix;
				MasterBin[q+NQ].a[0] = finalRix;
				char *finalT1 = "", *finalT2 = "";
				//uint32_t f1_l = UINT16_MAX, f2_l = UINT16_MAX; // full length
				char *finalT[2] = {finalT1,finalT2};
				int32_t finalL[2] = {UINT16_MAX,UINT16_MAX};
				//if (!max2 || (!doFullLCA && max > max2 && (double)max/tix >= conf)) { // TODO: revisit this earlyterm (specifically conf)
				if (!max2 || (max > max2 && (double)max/tix >= conf)) { // TODO: revisit this earlyterm (specifically conf)
					// early terminate -- call this taxon and don't do LCA
					if (HStr[0]) finalT[0] = HStr[0][HPairs[0][maxRix]];
					if (taxPath) MasterBin[q+NQ].a[1] = HPairs[0][maxRix];
					if (HStr[1]) {
						finalT[1] = HStr[1][HPairs[1][maxRix]];
						if (taxPath) MasterBin[q+NQ].a[2] = HPairs[1][maxRix];
					}
					if (doRedist) {
						C_bin[1][C_ix[1]++] = HStr[0] ? HPairs[0][maxRix] : -1;
						C_bin[1][C_ix[1]++] = -1;
						if (HStr[1]) C_bin[2][C_ix[2]++] = HPairs[1][maxRix],
							C_bin[2][C_ix[2]++] = -1;
					}
				} else for (int H = 0; H < 2 && HStr[H]; ++H) { // Interpolate the taxonomies. 
					// Consider max level first!
					tempix = 0;
					// Go through and tally refs into storage bin RBin (init'ed to 0)
					for (uint32_t i = 0; i < tix; ++i) {
						uint32_t pfx = SBin[i].p;
						uint64_t ix = SBin[i].s;
						
						uint64_t hardstop = Nibs[pfx+1];
						kmer_t prev_sfx = KGrid[ix].sfx;
						for (uint64_t j = ix; j < hardstop; ++j) {
							if (KGrid[j].sfx != prev_sfx) break;
							rix_t rix = KGrid[j].rix;
							rix_t h = HPairs[H][rix];
							if (!RBin[h]) TBin[tempix++] = h;
							if (RBin[h] >= 0) RBin[h] = -RBin[h] - 1;
						}
						for (int j = 0; j < tempix; ++j)
							RBin[TBin[j]] = RBin[TBin[j]] < 0 ? -RBin[TBin[j]] : RBin[TBin[j]];
					}
					int32_t h_max1 = 0, h_maxIx1 = -1,
						h_max2 = 0, h_maxIx2 = -1;
					// Consider all references that matched, uniquely
					for (int i = 0; i < tempix; ++i) { // tally to determine max and 2nd max ref
						uint32_t rix = TBin[i];
						if (RBin[rix] > h_max1 || RBin[rix]==h_max1 && rix < h_maxIx1) 
							h_max2 = h_max1, h_maxIx2 = h_maxIx1,
							h_max1 = RBin[rix], h_maxIx1 = rix;
						else if (RBin[rix] > h_max2)
							h_max2 = RBin[rix], h_maxIx2 = rix;
						//RBin[rix] = 0; // clear the bin
					}
					for (int i = 0; i < tempix; ++i) { // add best ref to redistribute pile; clear Ref count hash (RBin)
						uint32_t rix = TBin[i];
						if (doRedist && RBin[rix] == h_max1) // TODO: revisit adding a specificity control in here (i.e. only if !h_max2 or h_max1/h_max2 > 1.5...)
							C_bin[H+1][C_ix[H+1]++] = rix;
						RBin[rix] = 0; // clear the bin
					}
					if (doRedist) C_bin[H+1][C_ix[H+1]++] = -1;
					// Report if good....
					//if (perqPath) fprintf(outq,"-->%s\t[%u]\t%s\t%u\t%s\t%u\n",Qhed,tix,H1Str[h1_maxIx1],h1_max1,H1Str[h1_maxIx2],h1_max2);
						
					if (!h_max2 || (!doFullLCA && h_max1 > h_max2 && (double)h_max1/tix >= conf))  // TODO: revisit this early terminate (particularly conf)
					//if (!h_max2 || (h_max1 > h_max2 && (double)h_max1/tix >= conf))  // TODO: revisit this early terminate (particularly conf)
						finalT[H] = HStr[H][h_maxIx1];
					else { // We are in for pain. Do a full aufbau
						uint32_t agreed = tix; // Start out all agreeing at 0 semicolons!
						uint32_t ag_thres = conf*tix; 
						uint32_t winner = -1, winLv = -1;
						// TODO: implement blank-aware aufbau (remove blank "" taxa from consideration and total)
						//printf("%u matches in query %s\n",tix,Qhed);
						for (int semi = 1; agreed >= ag_thres; ++semi) {
							if (!LBins[H][semi-1]) break;
							agreed = 0; 
							int tempix = 0;
							for (uint32_t i = 0; i < tix; ++i) {
								uint64_t hardstop = Nibs[SBin[i].p+1];
								kmer_t prev_sfx = KGrid[SBin[i].s].sfx;
								for (uint64_t j = SBin[i].s; j < hardstop; ++j) {
									if (KGrid[j].sfx != prev_sfx) break;
									uint32_t h = HPairs[H][KGrid[j].rix];
									char *ref = HStr[H][h];
									if (LBins[H][semi-1][h]==-1) continue;
									uint32_t find = LBins[H][semi-1][h];
									//printf("[%u:%lu] Searching [lv %d]: %.*s\n --> Found: %ld: %s\n",
									//	i,j,semi,(int)(ptr-ref),ref,find,H1Str[find]);
									if (!RBin[find]) TBin[tempix++] = find;
									if (RBin[find] >= 0) RBin[find] = -RBin[find] - 1;
								}
								for (int j = 0; j < tempix; ++j)
									RBin[TBin[j]] = RBin[TBin[j]] < 0 ? -RBin[TBin[j]] : RBin[TBin[j]];
							}
							// The winner will be the count that is >= ag_thres.
							uint32_t local_max = 0, local_max2 = 0;
							uint32_t local_winner = 0;
							for (int i = 0; i < tempix; ++i) { // go thru unique options 
								uint32_t thisTax = TBin[i];
								if (RBin[thisTax] >= ag_thres) {
									if (RBin[thisTax] > local_max)
										local_max2 = local_max,
										local_max = RBin[thisTax],
										local_winner = thisTax;
									else if (RBin[thisTax] > local_max2) 
										local_max2 = RBin[thisTax];
								}
								RBin[thisTax] = 0; // also reset RBin
							}
							if (local_max > local_max2 && local_max >= ag_thres) 
								agreed = local_max, winner=local_winner, winLv = semi;
							//printf("Purity at lv %d: %u/%u [2nd: %u] (%s)\n --> Winner %d @ lv %d [%s]\n",
							//	semi,local_max,tix,local_max2,local_winner? H1Str[local_winner] : "[NONE]",winner,winLv,H1Str[winner]);
						}
						
						// All levels have been explored. Now let's expand the winner's string. 
						if (winner != -1) {
							
							if (perqPath) {
								int L = 0; char *p = HStr[H][winner]-1;
								for (int s = 0; s < winLv; ++s) 
									p = strchr(p+1,';');
								L = p - HStr[H][winner];
								finalL[H] = L;
								finalT[H] = HStr[H][winner];
							}
							if (taxPath) MasterBin[q+NQ].a[H+1] = winner + winLv * NUniqH[H];
						}
					}
				}
				
				DONE_TAXA:NULL;
				#pragma omp atomic
				++nAligns;
				if (perqPath) fprintf(outq,"%s\t%s\t[%u,%u]\t%.*s\t%.*s\t%u\n",Qhed,
					finalRix!=-1? RefNames[finalRix] : "",max,max2,finalL[0],finalT[0],
					finalL[1],finalT[1],tix);
			}
		}
		NQ += nq;
	}
	gzclose(in);
	printf(" - Total k-mers: %lu, non-error: %lu, matched: %lu\n",
		n_raw, n_filt, n_matched);
	printf("--> Successfully aligned %lu/%lu queries [%f].\n",nAligns,NQ,omp_get_wtime()-wtime);
	
	
	if (doRedist) {
		wtime = omp_get_wtime();
		printf("\n* Performing capitalist redistribution...\n");
		// The goal is to go multi-pass through each cap bin
		uint64_t Sizes[3] = {(refPath? numRef: 0),NUniqH[0],NUniqH[1]};
		int Pass[3] = {0};
		for (int i = 0; i < 3; ++i) {
			if (!Sizes[i]) continue;
			uint64_t *List = MasterList[i];
			uint8_t *Ix = MasterIx[i];
			uint64_t *Tally = huge_calloc(Sizes[i]*sizeof(*Tally)),
				*NextTally = huge_calloc(Sizes[i]*sizeof(*NextTally));
			if (!Tally || !NextTally) {puts("ERROR: OOM TALLY"); exit(3);}
			uint64_t debugSum_q = 0;
			for (uint64_t q = 0; q < NQ; ++q) {
				//if (MasterBin[q].a[0]==-1) continue;
				uint64_t *Bin = C_bins[Ix[q]][i]+List[q];
				debugSum_q += *Bin != -1;
				while (*Bin != -1) ++Tally[*Bin++];
			}
			//for (uint64_t j = 0; j < Sizes[i]; ++j) NextTally[j] = Tally[j];

			// Skip over blanks (disable with "--redistribute-blanks")
			// First, figure out which (if any) taxa are blanks. If multiple, stop at first.
			//printf("Beginning assignment, header %d\n",i); fflush(stdout);
			uint64_t firstIx = -1, next = 0;
			
			if (i > 0) {
				#pragma omp parallel for reduction (min:firstIx)
				for (uint64_t j = 0; j < Sizes[i]; ++j) 
					if (!HStr[i-1][j][0]) firstIx = firstIx < j ? firstIx : j;
				next = firstIx + 1; if (next >= Sizes[i]) next = 0;
				//printf("--> Found blank taxonomy string at H %d entry ix %lu (next: %s)\n",
				//	i,firstIx,HStr[i-1][next]); fflush(stdout);
			}
			
			uint64_t debugSum = 0;
			for (uint64_t j = 0; j < Sizes[i]; ++j) debugSum += Tally[j];
			//printf("Sum %lu [qsum %lu]; Doing passes now...\n",debugSum,debugSum_q); fflush(stdin); 
			uint64_t changes = -1, conv = NQ/100000, maxPass = 100;
			if (doFastRedist) maxPass = 1;
			for (int pass = 0; pass < maxPass && changes > conv; ++pass) {
				changes = 0;
				for (uint64_t q = 0; q < NQ; ++q) {
					//if (MasterBin[q].a[0]==-1) continue;
					uint64_t *Bin = C_bins[Ix[q]][i]+List[q];
					uint64_t maxTally = 0, whichMax = -1;
					while (*Bin != -1) {
						//if (Tally[*Bin] > maxTally) 
						if (Tally[*Bin] > maxTally && (*Bin != firstIx || whichMax == -1)) 
							maxTally = Tally[*Bin], whichMax = *Bin;
						++Bin;
					}
					
					if (whichMax != -1) ++NextTally[whichMax];
				}
				for (uint64_t j = 0; j < Sizes[i]; ++j) 
					changes += abs((int64_t)Tally[j] - (int64_t)NextTally[j]);
				//printf(" - pass: %d, changes: %lu\n",p, changes);
				
				uint64_t *t = Tally; Tally = NextTally; NextTally = t;
				memset(NextTally,0,sizeof(*NextTally)*Sizes[i]);
				++Pass[i];
			}
			
			//printf("Assigning distros...\n"); fflush(stdout);
			
			// Now assign final distributions
			for (uint64_t q = 0; q < NQ; ++q) {
				uint64_t *Bin = C_bins[Ix[q]][i]+List[q];
				uint64_t maxTally = 0, whichMax = -1;
				while (*Bin != -1) {
					if (Tally[*Bin] > maxTally && (*Bin != firstIx || whichMax == -1)) 
						maxTally = Tally[*Bin], whichMax = *Bin;
					++Bin;
				}
				MasterBin[q].a[i] = whichMax;
			}
			free(Tally); free(NextTally);
		}
		printf("--> Capitalist redistribution complete (%d,%d,%d) [%f]\n",
			Pass[0],Pass[1],Pass[2],omp_get_wtime()-wtime);
	}
	
	if (refPath) {
		wtime = omp_get_wtime();
		printf("\n* Printing reference tally table...\n");
		FILE *taxout = fopen(refPath,"wb");
		if (!taxout) {printf("ERROR: couldn't write to %s!\n",refPath); exit(2);}
		
		uint64_t *TaxTally = huge_calloc((uint64_t)numRef*sizeof(*TaxTally));
		for (uint64_t i = 0; i < NQ; ++i) // multithread?
			if (MasterBin[i].a[0] != -1) ++TaxTally[MasterBin[i].a[0]];
		for (uint32_t i = 0; i < numRef; ++i) if (TaxTally[i]) 
			fprintf(taxout,"%s\t%lu\n",RefNames[i],TaxTally[i]);
		free(TaxTally);
		fclose(taxout);
		printf("--> Reference table written [%f]\n",omp_get_wtime()-wtime);
	}
	
	if (taxPath) {
		wtime = omp_get_wtime();
		printf("\n* Printing taxonomy tally table...\n");
		FILE *taxout = fopen(taxPath,"wb");
		if (!taxout) {printf("ERROR: couldn't write to %s!\n",taxPath); exit(2);}
		// One option is to identify the max value for H1 and H2's values in MasterBin,
		// and allocate accordingly
		uint32_t maxH[2] = {0};
		for (uint64_t q = 0; q < NQ; ++q) {
			if (MasterBin[q].a[1] != -1 && MasterBin[q].a[1] > maxH[0])
				maxH[0] = MasterBin[q].a[1];
			if (MasterBin[q].a[2] != -1 && MasterBin[q].a[2] > maxH[1])
				maxH[1] = MasterBin[q].a[2];
		}
		//printf("Maximum tid for H1: %u, H2: %u\n",maxH[0],maxH[1]);
		
		// Indices >= NUniqH signify interpolation at lv = i/NUniqH
		
		// Create a bucket to tally the taxa
		for (int H = 0; H < 2; ++H) {
			if (!maxH[H]) continue;
			//printf("Number of unique headers to consider: %u, maxH = %u\n",NUniqH[H],maxH[H]);
			uint64_t *TaxTally = huge_calloc((uint64_t)(maxH[H]+1)*sizeof(*TaxTally));
			for (uint64_t i = 0; i < NQ; ++i) // multithread?
				if (MasterBin[i].a[H+1] != -1) ++TaxTally[MasterBin[i].a[H+1]];
			uint32_t lv = 0, nextLv = NUniqH[H];
			for (uint32_t i = 0; i <= maxH[H]; ++i) if (TaxTally[i]) { // 
				while (i >= nextLv) // Control depth of taxonomy reporting
					nextLv += NUniqH[H], ++lv;
				if (!lv) fprintf(taxout,"%s\t%lu\n",HStr[H][i],TaxTally[i]);
				else { // the interpolation case
					char *s = HStr[H][i-(nextLv-NUniqH[H])], *orig = s;
					for (uint32_t semi = 0; semi < lv; ++s)
						semi += *s==';';
					fprintf(taxout,"%.*s\t%lu\n",(uint32_t)(s-orig-1),orig,TaxTally[i]);
				}
			}
			free(TaxTally);
		}
		fclose(taxout);
		printf("--> Taxonomy written [%f]\n",omp_get_wtime()-wtime);
	} // end taxonomy tally output
	
	if (orthogPath) {
		wtime = omp_get_wtime();
		printf("\n* Printing orthogonal tally table...\n");
		if (!NUniqH[0] || !NUniqH[1]) 
			{printf("--> ERROR: Orthogonalizing requires 2 taxonomies.\n"); exit(2);}
		FILE *taxout = fopen(orthogPath,"wb");
		if (!taxout) {printf("ERROR: couldn't write to %s!\n",orthogPath); exit(2);}
		
		uint32_t maxH[2] = {0};
		for (uint64_t q = 0; q < NQ; ++q) {
			if (MasterBin[q].a[1] != -1 && MasterBin[q].a[1] > maxH[0])
				maxH[0] = MasterBin[q].a[1];
			if (MasterBin[q].a[2] != -1 && MasterBin[q].a[2] > maxH[1])
				maxH[1] = MasterBin[q].a[2];
		}
		
		//#define PRIME 1299827
		#define PRIME 4969
		uint64_t HashTable[PRIME+2] = {0};
		
		for (uint64_t q = 0; q < NQ; ++q) {
			if (MasterBin[q].a[0]==-1) continue;
			uint64_t val = *(uint64_t *)(MasterBin[q].b+4);
			++HashTable[(val % PRIME)+2];
		}
		for (uint64_t i = 1; i <= PRIME+1; ++i) HashTable[i] += HashTable[i-1];
		printf(" - Number added: %lu\n",HashTable[PRIME+1]);
		MasterBin_t **SortedIX = huge_calloc(HashTable[PRIME+1]*sizeof(*SortedIX));
		for (uint64_t q = 0; q < NQ; ++q) {
			if (MasterBin[q].a[0]==-1) continue;
			uint64_t val = *(uint64_t *)(MasterBin[q].b+4);
			SortedIX[HashTable[(val % PRIME)+1]++] = MasterBin + q;
		}
			
		for (uint64_t h = 0; h < PRIME; ++h) {
			uint64_t range = HashTable[h+1]-HashTable[h];
			MasterBin_t **ThisBinP = SortedIX+HashTable[h];
			if (range) qsort(ThisBinP,range,sizeof(*ThisBinP),binCmp);
		}
		// print result
		for (uint64_t h = 0; h < PRIME; ++h) {
			uint64_t range = HashTable[h+1]-HashTable[h];
			MasterBin_t **ThisBinP = SortedIX+HashTable[h];
			if (range) {
				uint64_t last = *(uint64_t *)(ThisBinP[0]->b+4);
				uint64_t tally = 0;
				for (uint64_t i = 0; i < range; ++i) {
					uint64_t val = *(uint64_t *)(ThisBinP[i]->b+4);
					if (val != last || i == range-1) { // commit
						uint32_t h1 = ThisBinP[i-1]->a[1], 
							h2 = ThisBinP[i-1]->a[2],
							lv1=h1/NUniqH[0], lv2=h2/NUniqH[1],
							L1 = 0, L2 = 0;
						for (int s = 0; s < lv1; ++L1) 
							s += HStr[0][h1][L1]==';';
						for (int s = 0; s < lv2; ++L2) 
							s += HStr[1][h2][L2]==';';
						L1 = L1 ?: UINT16_MAX;
						L2 = L2 ?: UINT16_MAX;
						fprintf(taxout,"%.*s\t%.*s\t%lu\n",L1-1,HStr[0][h1],L2-1,HStr[1][h2],tally);
						tally = 0;
					}
					++tally; 
					last = val; 
				}
			}
		}
		
		fclose(taxout);
		printf("--> Taxonomy written [%f]\n",omp_get_wtime()-wtime);
		
	}
	
	if (covPath) {
		wtime = omp_get_wtime();
		printf("\n* Printing coverage table...\n");
		uint64_t *TotK_m = huge_calloc((uint64_t)numRef*sizeof(*TotK_m));
		uint64_t *TotUniq_m = huge_calloc((uint64_t)numRef*sizeof(*TotUniq_m));
		uint64_t *FoundK_m = huge_calloc((uint64_t)numRef*sizeof(*FoundK_m));
		uint64_t *FoundUniq_m = huge_calloc((uint64_t)numRef*sizeof(*FoundUniq_m));
		uint32_t *PropK_m = huge_calloc((uint64_t)numRef*sizeof(*PropK_m));
		uint32_t *PropUniq_m = huge_calloc((uint64_t)numRef*sizeof(*PropUniq_m));

		//uint32_t *TotUniq_dedupe_m = huge_calloc((uint64_t)numRef*sizeof(*TotUniq_dedupe_m)); // new
		uint32_t *FoundUniq_dedupe_m = huge_calloc((uint64_t)numRef*sizeof(*FoundUniq_dedupe_m)); // new
		//uint32_t *PropUniq_dedupe_m = huge_calloc((uint64_t)numRef*sizeof(*PropUniq_dedupe_m)); // new
		// define a TotUniq_adamantium_m, FoundUniq_adamantium_m, PropUniq_adamantium_m
		uint32_t *TotUniq_adamantium_m = 0, *FoundUniq_adamantium_m = 0, *PropUniq_adamantium_m = 0;
		if (adamantium) { // whether to do adamantium hardened k-space scouring
			TotUniq_adamantium_m = huge_calloc((uint64_t)numRef*sizeof(*TotUniq_adamantium_m)); // new
			FoundUniq_adamantium_m = huge_calloc((uint64_t)numRef*sizeof(*FoundUniq_adamantium_m)); // new
			PropUniq_adamantium_m = huge_calloc((uint64_t)numRef*sizeof(*PropUniq_adamantium_m)); // new
		}
		#pragma omp parallel
		{
			int tid = omp_get_thread_num();
			uint64_t *TotK = TotK_m, *TotUniq = TotUniq_m;
			
			uint64_t *FoundK = FoundK_m, *FoundUniq = FoundUniq_m;
			uint32_t *PropK=PropK_m, *PropUniq=PropUniq_m;
			uint32_t //*TotUniq_dedupe=TotUniq_dedupe_m, 
				*FoundUniq_dedupe = FoundUniq_dedupe_m; // new
			//uint32_t *PropUniq_dedupe = PropUniq_dedupe_m; // new
			uint32_t *TotUniq_adamantium = TotUniq_adamantium_m, *FoundUniq_adamantium = FoundUniq_adamantium_m, *PropUniq_adamantium = PropUniq_adamantium_m; // new

			//char kmer_str[64] = {0}, rc_str[64] = {0};
			
			if (tid) 
				TotK = huge_calloc((uint64_t)numRef*sizeof(*TotK)),
				TotUniq = huge_calloc((uint64_t)numRef*sizeof(*TotUniq)),
				FoundK = huge_calloc((uint64_t)numRef*sizeof(*FoundK)),
				FoundUniq = huge_calloc((uint64_t)numRef*sizeof(*FoundUniq)),
				PropK = huge_calloc((uint64_t)numRef*sizeof(*PropK)),
				PropUniq = huge_calloc((uint64_t)numRef*sizeof(*PropUniq)),
				FoundUniq_dedupe = huge_calloc((uint64_t)numRef*sizeof(*FoundUniq_dedupe));
			if (tid && adamantium)
				TotUniq_adamantium = huge_calloc((uint64_t)numRef*sizeof(*TotUniq_adamantium)),
				FoundUniq_adamantium = huge_calloc((uint64_t)numRef*sizeof(*FoundUniq_adamantium)),
				PropUniq_adamantium = huge_calloc((uint64_t)numRef*sizeof(*PropUniq_adamantium));
				//TotUniq_dedupe = huge_calloc((uint64_t)numRef*sizeof(*TotUniq_dedupe)),
				//PropUniq_dedupe = huge_calloc((uint64_t)numRef*sizeof(*PropUniq_dedupe)); // new
			
			#pragma omp for schedule(dynamic,1024)
			for (uint64_t i = 0; i < nbins; ++i) {
				if (Nibs[i+1] <= Nibs[i]) continue; // empty
				uint32_t ambig = 0;
				uint64_t end = Nibs[i+1], nd = 1, range = end-Nibs[i];
				uint64_t mv = 0; // max count of this kmer
				kmer_t thisK = KGrid[Nibs[i]].sfx+1;
				
				for (uint64_t j = Nibs[i]; j < end; j += nd) { // for each kmer suffix within this prefix
					rix_t rix = KGrid[j].rix;
					if (rix >= numRef) {
						printf("ERROR: rix at bin %lu, nib %lu = %lu, >= numRef [%lu]. nbins = %lu\n",
							i, j, rix, numRef,nbins);
						continue;
					}
					
					// If new k-mer, check if ambig, store max value
					if (KGrid[j].sfx != thisK) {
						thisK = KGrid[j].sfx; 
						ambig = 0; mv = QueryAligns[j];
						for (uint64_t k = j+1; k < end && KGrid[k].sfx == thisK; ++k)
							mv = mv > QueryAligns[k] ? mv : QueryAligns[k],
							ambig |= KGrid[k].rix ^ rix;
					}

					// Find number of in-ref copies
					nd = 1;
					for (uint64_t k = j+1; k < end && 
					KGrid[k].sfx == thisK && KGrid[k].rix == rix; ++k) ++nd;
					
					// Increment the appropriate variables
					if (!ambig) {
						TotUniq[rix] += nd;
						FoundUniq[rix] += mv;
						PropUniq[rix] += mv < nd ? mv : nd;

						/* if (mv > 0) {
							char prefix_str[64] = {0}, suffix_str[64] = {0};
							num_to_dna(i, PL, prefix_str);
							num_to_dna(thisK, SL, suffix_str);
							printf("%s %s nd: %lu, mv: %lu\n", prefix_str, suffix_str, nd, mv);
						} */
						

						FoundUniq_dedupe[rix] += mv/nd; // new
						//TotUniq_dedupe[rix]++; // new 
						//PropUniq_dedupe[rix] += mv > 0;
						if (adamantium) {
							// First look up this kmer in the adamantium table
							uint64_t ix = skLBS(Adamant_Suffix + Adamant_Prefix[i],Adamant_Prefix[i+1]-Adamant_Prefix[i], thisK);
							//printf("Looking up %lu in adamantium table for %lu (range from %lu to %lu): %lu\n",thisK,i,Adamant_Prefix[i],Adamant_Prefix[i+1],ix);	

							if (ix == -1) {
								TotUniq_adamantium[rix] += nd; // new
								FoundUniq_adamantium[rix] += (uint32_t)( (double)mv/(double)nd + 0.5 ); // new
								PropUniq_adamantium[rix] += mv < nd ? mv : nd; // new
							}
						}
					}
					TotK[rix] += nd;
					FoundK[rix] += mv;
					PropK[rix] += mv < nd ? mv : nd;
				}
			}
			#pragma omp critical
			if (tid) {
				for (uint32_t i = 0; i < numRef; ++i) {
					TotK_m[i] += TotK[i], TotUniq_m[i] += TotUniq[i],
					FoundK_m[i] += FoundK[i], FoundUniq_m[i] += FoundUniq[i],
					PropK_m[i] += PropK[i], PropUniq_m[i] += PropUniq[i];
					FoundUniq_dedupe_m[i] += FoundUniq_dedupe[i];
					if (adamantium)
						TotUniq_adamantium_m[i] += TotUniq_adamantium[i], // new
					FoundUniq_adamantium_m[i] += FoundUniq_adamantium[i], // new
					PropUniq_adamantium_m[i] += PropUniq_adamantium[i]; // new
					//TotUniq_dedupe_m[i] += TotUniq_dedupe[i]; // new
					//PropUniq_dedupe_m[i] += PropUniq_dedupe[i]; // new
				}
			}
		}
		
		printf(" - Done with coverage loop. [%f]\n",omp_get_wtime() - wtime);
		printf(" - Beginning accumulation...\n");
		
		uint64_t totFoundKUniq = TotUniq_m[numRef-1];
		for (uint64_t i = 1; i < numRef; ++i) totFoundKUniq+=TotUniq_m[i-1];
		printf(" - Total k-mers found: %lu (%lu unique)\n",numK,totFoundKUniq);
		
		//wtime = omp_get_wtime();
		
		double qkAvg = (double)n_raw / NQ;

		FILE *out = fopen(covPath,"wb");
		if (!out) {puts("ERROR: can't open output file!"); exit(2);}
		setvbuf(out, 0, _IOFBF, 1<<22);
		fprintf(out, "Reference\tRef_len\tTotal_kmers\tTotal_unique\tKmers_found\t");
		fprintf(out, "Unique_kmers_found\tKmers_covered\tUnique_kmers_covered\t");

		fprintf(out, "Proportion_covered\tUnique_proportion_covered\tReads_covered\t");
		fprintf(out, "Exp\tSuspicion\tCoverage_est\tNum_reads_est"); //\tStrain\n");
		if (adamantium) fprintf(out, "\tAdamantine_tot\tAdamantium_found\tAdamantium_covered"),
			fprintf(out, "\tAdamantium_prop\tXCov_adamantium\tReads_covered_adamantium\n");
		else fprintf(out,"\n");
		//FILE *log = fopen("DbgPost.txt","wb");
		for (uint32_t i = 0; i < numRef; ++i) {
			//fprintf(log,"%u\t%s\n",i,RefNames[i]);
			if (!FoundK_m[i]) continue;
			
			double xcov = (double)FoundUniq_dedupe_m[i]/TotUniq_m[i], ucov = (double)PropUniq_m[i]/TotUniq_m[i];
			double //ucov_dd = (double)PropUniq_dedupe_m[i]/TotUniq_dedupe_m[i], 
				pcov = (double)PropK_m[i]/TotK_m[i];
			fprintf(out,"%s\t%u\t%lu\t%lu\t%u\t%u\t%u\t%u", RefNames[i],RefSizes[i],TotK_m[i], TotUniq_m[i], 
				//TotUniq_dedupe_m[i], 
				FoundK_m[i], FoundUniq_m[i],PropK_m[i],PropUniq_m[i]); //, ucov_dd);
			
			fprintf(out,"\t%f\t%f\t%u",
				pcov, ucov, FullQueryAligns[i]);

			double expc = 1 - exp(-xcov), susp = ucov ? log10(expc / ucov) : 0;
			fprintf(out,"\t%f\t%f\t%f\t%u", expc, susp > 0 ? susp : 0, xcov, (uint32_t)(xcov*(RefSizes[i])/qkAvg + 0.5));
				//,pcov/ucov > 0.9 && susp < .055 && ucov > .10 && (ucov > .5 || ucov*PropUniq_m[i] >= 5));
			if (adamantium) {
				double xcov_adamantium = (double)FoundUniq_adamantium_m[i]/TotUniq_adamantium_m[i], 
					ucov_adamantium = (double)PropUniq_adamantium_m[i]/TotUniq_adamantium_m[i];
				fprintf(out,"\t%u\t%u\t%u\t%f\t%f\t%u\n", TotUniq_adamantium_m[i], FoundUniq_adamantium_m[i], PropUniq_adamantium_m[i], 
					ucov_adamantium,xcov_adamantium, (uint32_t)(xcov_adamantium*(RefSizes[i])/qkAvg + 0.5));
			} else fprintf(out,"\n");
		}
		printf("--> Coverage table written [%f]\n",omp_get_wtime()-wtime);
		fclose(out);
	} // end coverage output
	
	printf("\nAll outputs written.\n");
	return 0;
}
