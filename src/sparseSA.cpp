#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <limits.h>

#include "sparseSA.hpp"

// LS suffix sorter (integer alphabet). 
extern "C" { void suffixsort(int *x, int *p, int n, int k, int l); };

pthread_mutex_t cout_mutex = PTHREAD_MUTEX_INITIALIZER;

sparseSA::sparseSA(string &S_, vector<string> &descr_, vector<long> &startpos_, bool __4column, long K_) : 
  descr(descr_), startpos(startpos_), S(S_) {
  _4column = __4column;

  // Get maximum query sequence description length.
  maxdescrlen = 0;
  for(long i = 0; i < (long)descr.size(); i++) {
    if(maxdescrlen < (long)descr[i].length())  maxdescrlen = descr[i].length();
  }
  K = K_;

  // Increase string length so divisible by K. 
  // Don't forget to count $ termination character. 
  if(S.length() % K != 0) {
    long appendK = K - S.length() % K ;
    for(long i = 0; i < appendK; i++) S += '$';
  }
  // Make sure last K-sampled characeter is this special character as well!!
  for(long i = 0; i < K; i++) S += '$'; // Append "special" end character. Note: It must be lexicographically less.
  N = S.length();

  if(K > 1) {
    long bucketNr = 1;
    int *intSA = new int[N/K+1];  for(int i = 0; i < N/K; i++) intSA[i] = i; // Init SA.
    int* t_new = new int[N/K+1];
    long* BucketBegin = new long[256]; // array to save current bucket beginnings
    radixStep(t_new, intSA, bucketNr, BucketBegin, 0, N/K-1, 0); // start radix sort
    t_new[N/K] = 0; // Terminate new integer string.
    delete[] BucketBegin;

    // Suffix sort integer text.
    cerr << "# suffixsort()" << endl;
    suffixsort(t_new, intSA, N/K, bucketNr, 0);
    cerr << "# DONE suffixsort()" << endl;

    delete[] t_new;

    // Translate suffix array. 
    SA.resize(N/K);
    for (long i=0; i<N/K; i++) SA[i] = (unsigned int)intSA[i+1] * K;
    delete[] intSA;

    // Build ISA using sparse SA. 
    ISA.resize(N/K);             
    for(long i = 0; i < N/K; i++) { ISA[SA[i]/K] = i; }
  }
  else {
    SA.resize(N);
    ISA.resize(N);
    int char2int[UCHAR_MAX+1]; // Map from char to integer alphabet.

    // Zero char2int mapping.
    for (int i=0; i<=UCHAR_MAX; i++) char2int[i]=0;

    // Determine which characters are used in the string S.
    for (long i = 0; i < N; i++) char2int[(int)S[i]]=1;

    // Count the size of the alphabet. 
    int alphasz = 0; 
    for(int i=0; i <= UCHAR_MAX; i++) {
      if (char2int[i]) char2int[i]=alphasz++;
      else char2int[i] = -1;
    }

    // Remap the alphabet. 
    for(long i = 0; i < N; i++) ISA[i] = (int)S[i]; 
    for (long i = 0; i < N; i++) ISA[i]=char2int[ISA[i]] + 1; 
    // First "character" equals 1 because of above plus one, l=1 in suffixsort(). 
    int alphalast = alphasz + 1;
    
    // Use LS algorithm to construct the suffix array.
    int *SAint = (int*)(&SA[0]);
    suffixsort(&ISA[0], SAint , N-1, alphalast, 1);
  }

  cerr << "N=" << N << endl;

  // Adjust to "sampled" size. 
  logN = (long)ceil(log(N/K) / log(2.0));

  LCP.resize(N/K);
  cerr << "N/K=" << N/K << endl;
  // Use algorithm by Kasai et al to construct LCP array.
  computeLCP();  // SA + ISA -> LCP
  LCP.init();

  NKm1 = N/K-1;
}

// Uses the algorithm of Kasai et al 2001 which was described in
// Manzini 2004 to compute the LCP array. Modified to handle sparse
// suffix arrays and inverse sparse suffix arrays.
void sparseSA::computeLCP() {
  long h=0;
  for(long i = 0; i < N; i+=K) { 
    long m = ISA[i/K]; 
    if(m==0) LCP.set(m, 0); // LCP[m]=0;
    else {
      long j = SA[m-1];
      while(i+h < N && j+h < N && S[i+h] == S[j+h])  h++;
      LCP.set(m, h); //LCP[m] = h;
    }
    h = max(0L, h - K);
  }
}


// Implements a variant of American flag sort (McIlroy radix sort).
// Recurse until big-K size prefixes are sorted. Adapted from the C++
// source code for the wordSA implementation from the following paper:
// Ferragina and Fischer. Suffix Arrays on Words. CPM 2007.
void sparseSA::radixStep(int *t_new, int *SA, long &bucketNr, long *BucketBegin, long l, long r, long h) {
  if(h >= K) return;
  // first pass: count
  vector<long> Sigma(256, 0); // Sigma counts occurring characters in bucket
  for (long i = l; i <= r; i++) Sigma[ S[ SA[i]*K + h ] ]++; // count characters
  BucketBegin[0] = l; for (long i = 1; i < 256; i++) { BucketBegin[i] = Sigma[i-1] + BucketBegin[i-1]; } // accumulate

  // second pass: move (this variant does *not* need an additional array!)
  unsigned char currentKey = 0;    // character of current bucket
  long end = l-1+Sigma[currentKey]; // end of current bucket
  long pos = l;                     // 'pos' is current position in bucket
  while (1) {
    if (pos > end) { // Reached the end of the bucket.
      if (currentKey == 255) break; // Last character?
      currentKey++; // Advance to next characer.
      pos = BucketBegin[currentKey]; // Next bucket start.
      end += Sigma[currentKey]; // Next bucket end.
    }
    else {
      // American flag sort of McIlroy et al. 1993. BucketBegin keeps
      // track of current position where to add to bucket set.
      int tmp = SA[ BucketBegin[ S[ SA[pos]*K + h ] ] ]; 
      SA[ BucketBegin[ S[ SA[pos]*K + h] ]++ ] = SA[pos];  // Move bucket beginning to the right, and replace 
      SA[ pos ] = tmp; // Save value at bucket beginning.
      if (S[ SA[pos]*K + h ] == currentKey) pos++; // Advance to next position if the right character.
    }
  }
  // recursively refine buckets and calculate new text:
  long beg = l; end = l-1;
  for (long i = 1; i < 256; i++) { // step through Sigma to find bucket borders
    end += Sigma[i];
    if (beg <= end) {
      if(h == K-1) {
	for (long j = beg; j <= end; j++) {
	  t_new[ SA[j] ] = bucketNr; // set new text
	}
	bucketNr++;
      }
      else {
	radixStep(t_new, SA, bucketNr, BucketBegin, beg, end, h+1); // recursive refinement
      }
      beg = end + 1; // advance to next bucket
    }
  }
}

// Binary search for left boundry of interval.
long sparseSA::bsearch_left(char c, long i, long s, long e) {
  if(c == S[SA[s]+i]) return s;
  long l = s, r = e;
  while (r - l > 1) {
    long m = (l+r) / 2;
    if (c <= S[SA[m] + i]) r = m;
    else l = m;
  }
  return r;
}

// Binary search for right boundry of interval.
long sparseSA::bsearch_right(char c, long i, long s, long e) {
  if(c == S[SA[e]+i]) return e;
  long l = s, r = e;
  while (r - l > 1) {
    long m = (l+r) / 2;
    if (c < S[SA[m] + i]) r = m;
    else l = m;
  }
  return l;
}


// Simple top down traversal of a suffix array.
bool sparseSA::top_down(char c, long i, long &start, long &end) {
  if(c < S[SA[start]+i]) return false;
  if(c > S[SA[end]+i]) return false;
  long l = bsearch_left(c, i, start, end);
  long l2 = bsearch_right(c, i, start, end);
  start = l; end = l2;
  return l <= l2;
}

// Top down traversal of the suffix array to match a pattern.  NOTE:
// NO childtab as in the enhanced suffix array (ESA).
bool sparseSA::search(string &P, long &start, long &end) {
  start = 0; end = N - 1;
  long i = 0;
  while(i < (long)P.length()) {
    if(top_down(P[i], i, start, end) == false) {
      return false;
    }
    i++;
  }
  return true;
}


// Traverse pattern P starting from a given prefix and interval
// until mismatch or min_len characters reached.
void sparseSA::traverse(string &P, long prefix, interval_t &cur, int min_len) {
  if(cur.depth >= min_len) return;

  while(prefix+cur.depth < (long)P.length()) {
    long start = cur.start; long end = cur.end;
    // If we reach a mismatch, stop.
    if(top_down_faster(P[prefix+cur.depth], cur.depth, start, end) == false) return;

    // Advance to next interval.
    cur.depth += 1; cur.start = start; cur.end = end;
    
    // If we reach min_len, stop.
    if(cur.depth == min_len) return;
  }
}


// Given SA interval apply binary search to match character c at
// position i in the search string. Adapted from the C++ source code
// for the wordSA implementation from the following paper: Ferragina
// and Fischer. Suffix Arrays on Words. CPM 2007.
bool sparseSA::top_down_faster(char c, long i, long &start, long &end) {
  long l, r, m, r2=end, l2=start, vgl;
  bool found = false;
  long cmp_with_first = (long)c - (long)S[SA[start]+i];
  long cmp_with_last = (long)c - (long)S[SA[end]+i];
  if(cmp_with_first < 0) { 
    l = start+1; l2 = start; // pattern doesn't occur!  
  }
  else if(cmp_with_last > 0) { 
    l = end+1; l2 = end; 
    // pattern doesn't occur!  
  }
  else {
    // search for left border:
    l = start; r = end;
    if (cmp_with_first == 0) {
      found = true; r2 = r;
    }
    else {
      while (r - l > 1) {
	m = (l+r) / 2;
	vgl = (long)c - (long)S[SA[m] + i];
	if (vgl <= 0) {
	  if (!found && vgl == 0) {
	    found = true;
	    l2 = m; r2 = r; // search interval for right border
	  }
	  r = m;
	}
	else l = m;
      }
      l = r;
    }
    // search for right border (in the range [l2:r2])
    if (!found) { 
      l2 = l - 1; // pattern not found => right border to the left of 'l' 
    } 
    if (cmp_with_last == 0) { 
      l2 = end; // right border is the end of the array 
    }
    else {
      while (r2 - l2 > 1) {
	m = (l2 + r2) / 2;
	vgl = (long)c - (long)S[SA[m] + i];
	if (vgl < 0) r2 = m;
	else l2 = m;
      }
    }
  }
  start = l;
  end = l2;
  return l <= l2;
}


// Suffix link simulation using ISA/LCP heuristic.
bool sparseSA::suffixlink(interval_t &m) {
  m.depth -= K;
  if( m.depth <= 0) return false;
  m.start = ISA[SA[m.start] / K + 1];  
  m.end = ISA[SA[m.end] / K + 1]; 
  return expand_link(m);
}

// For a given offset in the prefix k, find all MEMs.
void sparseSA::findMEM(long k, string &P, vector<match_t> &matches, int min_len, bool print) {
  if(k < 0 || k >= K) { cerr << "Invalid k." << endl; return; }
  // Offset all intervals at different start points.
  long prefix = k;
  interval_t mli(0,N/K-1,0); // min length interval
  interval_t xmi(0,N/K-1,0); // max match interval
  
  // Right-most match used to terminate search.
  int min_lenK = min_len - (K-1);

  while( prefix <= (long)P.length() - (K-k)) {
    traverse(P, prefix, mli, min_lenK);    // Traverse until minimum length matched.
    if(mli.depth > xmi.depth) xmi = mli;
    if(mli.depth <= 1) { mli.reset(N/K-1); xmi.reset(N/K-1); prefix+=K; continue; }

    if(mli.depth >= min_lenK) {   
      traverse(P, prefix, xmi, P.length()); // Traverse until mismatch.
      collectMEMs(P, prefix, mli, xmi, matches, min_len, print); // Using LCP info to find MEM length.
      // When using ISA/LCP trick, depth = depth - K. prefix += K. 
      prefix+=K;	
      if( suffixlink(mli) == false ) { mli.reset(N/K-1); xmi.reset(N/K-1); continue; }
      suffixlink(xmi);
    }
    else {
      // When using ISA/LCP trick, depth = depth - K. prefix += K. 
      prefix+=K;	
      if( suffixlink(mli) == false ) { mli.reset(N/K-1); xmi.reset(N/K-1); continue; }
      xmi = mli;
    }
  }
  if(print) print_match(match_t(), matches);   // Clear buffered matches.
}


// Use LCP information to locate right maximal matches. Test each for
// left maximality.
void sparseSA::collectMEMs(string &P, long prefix, interval_t mli, interval_t xmi, vector<match_t> &matches, int min_len, bool print) {
  // All of the suffixes in xmi's interval are right maximal.
  for(long i = xmi.start; i <= xmi.end; i++) find_Lmaximal(P, prefix, SA[i], xmi.depth, matches, min_len, print);

  if(mli.start == xmi.start && mli.end == xmi.end) return;

  while(xmi.depth >= mli.depth) {
    // Attempt to "unmatch" xmi using LCP information.
    if(xmi.end+1 < N/K) xmi.depth = max(LCP[xmi.start], LCP[xmi.end+1]);
    else xmi.depth = LCP[xmi.start];

    // If unmatched XMI is > matched depth from mli, then examine rmems.
    if(xmi.depth >= mli.depth) {
      // Scan RMEMs to the left, check their left maximality..
      while(LCP[xmi.start] >= xmi.depth) { 
	xmi.start--; 
	find_Lmaximal(P, prefix, SA[xmi.start], xmi.depth, matches, min_len, print);
      }
      // Find RMEMs to the right, check their left maximality.
      while(xmi.end+1 < N/K && LCP[xmi.end+1] >= xmi.depth) { 
	xmi.end++;  
	find_Lmaximal(P, prefix, SA[xmi.end], xmi.depth, matches, min_len, print);
      }
    }
  }
}


// Finds left maximal matches given a right maximal match at position i.
void sparseSA::find_Lmaximal(string &P, long prefix, long i, long len, vector<match_t> &matches, int min_len, bool print) {
  // Advance to the left up to K steps.
  for(long k = 0; k < K; k++) {
    // If we reach the end and the match is long enough, print.
    if(prefix == 0 || i == 0) {
      if(len >= min_len) {
	if(print) print_match(match_t(i, prefix, len), matches);
	else matches.push_back(match_t(i, prefix, len));
      }
      return; // Reached mismatch, done.
    }
    else if(P[prefix-1] != S[i-1]){
      // If we reached a mismatch, print the match if it is long enough.
      if(len >= min_len) {
	if(print) print_match(match_t(i, prefix, len), matches);
	else matches.push_back(match_t(i, prefix, len));
      }
      return; // Reached mismatch, done.
    }
    prefix--; i--; len++; // Continue matching.
  }
}


// Print results in format used by MUMmer v3.  Prints results
// 1-indexed, instead of 0-indexed.
void sparseSA::print_match(match_t m) {
  if(_4column == false) {
    printf("%8ld  %8ld  %8ld\n", m.ref + 1, m.query + 1, m.len);
  }
  else {
    long refseq=0, refpos=0;
    from_set(m.ref, refseq, refpos); // from_set is slow!!!
    // printf works faster than count... why? I don't know!!
    printf("  %s", descr[refseq].c_str());
    for(long j = 0; j < maxdescrlen - (long)descr[refseq].length() + 1; j++) putchar(' '); 
    printf(" %8ld  %8ld  %8ld\n", refpos + 1L, m.query + 1L, m.len);
  }
}

// This version of print match places m_new in a buffer. The buffer is
// flushed if m_new.len <= 0 or it reaches 1000 entries.  Buffering
// limits the number of locks on cout.
void sparseSA::print_match(match_t m_new, vector<match_t> &buf) {
  if(m_new.len > 0)  buf.push_back(m_new);
  if(buf.size() > 1000 || m_new.len <= 0) {
    pthread_mutex_lock(&cout_mutex);
    for(long i = 0; i < (long)buf.size(); i++) print_match(buf[i]);
    pthread_mutex_unlock(&cout_mutex);
    buf.clear();
  }
}

void sparseSA::print_match(string meta, vector<match_t> &buf, bool rc) {
  pthread_mutex_lock(&cout_mutex);
  if(!rc) printf("> %s\n", meta.c_str());
  else printf("> %s Reverse\n", meta.c_str());
  for(long i = 0; i < (long)buf.size(); i++) print_match(buf[i]);
  pthread_mutex_unlock(&cout_mutex);
  buf.clear();
}

// Finds maximal almost-unique matches (MAMs) These can repeat in the
// given query pattern P, but occur uniquely in the indexed reference S.
void sparseSA::findMAM(string &P, vector<match_t> &matches, int min_len, bool print) {
  interval_t cur(0, N-1, 0);
  long prefix = 0; 
  while(prefix < (long)P.length()) {
    // Traverse SA top down until mismatch or full string is matched.
    traverse(P, prefix, cur, P.length());
    if(cur.depth <= 1) { cur.depth = 0; cur.start = 0; cur.end = N-1; prefix++; continue; }
    if(cur.size() == 1 && cur.depth >= min_len) { 
      if(is_leftmaximal(P, prefix, SA[cur.start])) {
	// Yes, it's a MAM.
	match_t m; m.ref = SA[cur.start]; m.query = prefix; m.len = cur.depth;
	if(print) print_match(m);
	else  matches.push_back(m); 
      }
    }
    do {
      cur.depth = cur.depth-1;
      cur.start = ISA[SA[cur.start] + 1];  
      cur.end = ISA[SA[cur.end] + 1]; 
      prefix++;
      if( cur.depth == 0 || expand_link(cur) == false ) { cur.depth = 0; cur.start = 0; cur.end = N-1; break; }
    } while(cur.depth > 0 && cur.size() == 1);
  }
}

// Returns true if the position p1 in the query pattern and p2 in the
// reference is left maximal.
bool sparseSA::is_leftmaximal(string &P, long p1, long p2) {
  if(p1 == 0 || p2 == 0) return true;
  else return P[p1-1] != S[p2-1];
}


struct by_ref { bool operator() (const match_t &a, const match_t &b) const { if(a.ref == b.ref) return a.len > b.len; else return a.ref < b.ref; }  };

// Maximal Unique Match (MUM) 
void sparseSA::MUM(string &P, vector<match_t> &unique, int min_len, bool print) {
  // Find unique MEMs.
  vector<match_t> matches;
  MAM(P, matches, min_len, false);

  // Adapted from Stephan Kurtz's code in cleanMUMcand.c in MUMMer v3.20. 
  long currentright, dbright = 0;
  bool ignorecurrent, ignoreprevious = false;
  sort(matches.begin(), matches.end(), by_ref());
  for(long i = 0; i < (long)matches.size(); i++) {
    ignorecurrent = false;
    currentright = matches[i].ref + matches[i].len - 1;
    if(dbright > currentright) 
      ignorecurrent = true;
    else {
      if(dbright == currentright) {
	ignorecurrent = true;
	if(!ignoreprevious && matches[i-1].ref == matches[i].ref) 
	  ignoreprevious = true;
      }
      else {
	dbright = currentright;
      }
    }
    if(i > 0 && !ignoreprevious) {
      if(print)	print_match(matches[i-1]);
      else unique.push_back(matches[i-1]);
    }
    ignoreprevious = ignorecurrent;
  }
  if(!ignoreprevious) {
    if(matches.size() > 0) {
      if(print) print_match(matches[matches.size()-1]);
      else unique.push_back(matches[matches.size()-1]);
    }
  }

}

struct thread_data {
  vector<long> Kvalues; // Values of K this thread should process.
  sparseSA *sa; // Suffix array + aux informaton
  int min_len; // Minimum length of match.
  string *P; // Query string.
};

void *MEMthread(void *arg) {
  thread_data *data = (thread_data*)arg;
  vector<long> &K = data->Kvalues;
  sparseSA *sa = data->sa;

  // Find MEMs for all assigned offsets to this thread.

  vector<match_t> matches; // place-holder
  matches.reserve(2000);   // TODO: Use this as a buffer again!!!!!! 

  for(long k = 0; k < (long)K.size(); k++) {  sa->findMEM(K[k], *(data->P), matches, data->min_len, true); }

  pthread_exit(NULL);
}

void sparseSA::MEM(string &P, vector<match_t> &matches, int min_len, bool print, int num_threads) {
  if(min_len < K) return;
  if(num_threads == 1) {
    for(int k = 0; k < K; k++) { findMEM(k, P, matches, min_len, print); }
  }
  else if(num_threads > 1) {
    vector<pthread_t> thread_ids(num_threads);  
    vector<thread_data> data(num_threads);

    // Make sure all num_threads are joinable.
    pthread_attr_t attr;  pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // Distribute K-values evenly between num_threads.
    int t = 0;
    for(int k = 0; k < K; k++) {
      data[t].Kvalues.push_back(k);
      t++;
      if(t == num_threads) t = 0;
    }
    // Initialize additional thread data.
    for(int i = 0; i < num_threads; i++) { data[i].sa = this; data[i].min_len = min_len;  data[i].P = &P; }
    // Create joinable threads to find MEMs.
    for(int i = 0; i < num_threads; i++) pthread_create(&thread_ids[i], &attr, MEMthread, (void *)&data[i]);
    // Wait for all threads to terminate.
    for(int i = 0; i < num_threads; i++) pthread_join(thread_ids[i], NULL);    
  }
}

