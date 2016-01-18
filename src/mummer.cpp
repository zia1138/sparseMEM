#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include "sparseSA.hpp"
#include "fasta.hpp"

#include <getopt.h>
#include <time.h>
#include <sys/time.h>

#include <cctype> // std::tolower(), uppercase/lowercase conversion

// NOTE use of special characters ~, `, and $ !!!!!!!!

using namespace std;

void usage(string prog);

enum mum_t { MUM, MAM, MEM };

int min_len = 20;
mum_t type = MAM;
bool rev_comp = false, _4column = false, nucleotides_only = false;
int K = 1, num_threads = 1, query_threads = 1;
sparseSA *sa;
string query_fasta;

struct query_arg {
  int skip0;
  int skip;
};

void *query_thread(void *arg_) {
  query_arg *arg = (query_arg *)arg_;

  string meta, line;
  ifstream data(query_fasta.c_str());

  vector<match_t> matches;

  bool print = arg->skip == 1;

  long seq_cnt = 0;

  if(!data.is_open()) { cerr << "unable to open " << query_fasta << endl; exit(1); }

  // Collect meta data.
  while(!data.eof()) {
    getline(data, line); // Load one line at a time.
    if(line.length() == 0) continue;
    if(line[0] == '>') {
      long start = 1, end = line.length() - 1;
      trim(line, start, end);
      for(long i = start; i <= end; i++) {
	if( line[i] == ' ') break; // Behave like MUMmer 3 cut off meta after first space.
	meta += line[i];
      }
      cerr << "# " << meta << endl;
      break;
    }
  }
  string *P = new string;
  while(!data.eof()) {
    getline(data, line); // Load one line at a time.
    if(line.length() == 0) continue;
    long start = 0, end = line.length() - 1;
    // Meta tag line and start of a new sequence.
    // Collect meta data.
    if(line[0] == '>') {
      if(meta != "") {
	if(seq_cnt % arg->skip == arg->skip0) {
	  // Process P.
	  cerr << "# P.length()=" << P->length() << endl;
	  if(print) printf("> %s\n", meta.c_str());
	  if(type == MAM) sa->MAM(*P, matches, min_len, print);
	  else if(type == MUM) sa->MUM(*P, matches, min_len, print);
	  else if(type == MEM) sa->MEM(*P, matches, min_len, print, num_threads);
	  if(!print) sa->print_match(meta, matches, false); 
	  if(rev_comp) {
	    reverse_complement(*P, nucleotides_only);
	    if(print) printf("> %s Reverse\n", meta.c_str());
	    if(type == MAM) sa->MAM(*P, matches, min_len, print);
	    else if(type == MUM) sa->MUM(*P, matches, min_len, print);
	    else if(type == MEM) sa->MEM(*P, matches, min_len, print, num_threads);
	    if(!print) sa->print_match(meta, matches, true); 
	  }
	}
	seq_cnt++;
        delete P; P = new string; meta = ""; 
      }
      start = 1;
      trim(line, start, end);
      for(long i = start; i <= end; i++) {
	if(line[i] == ' ') break; // Behave like MUMmer 3 cut of meta after first space.
	meta += line[i];
      }
      cerr << "# " << meta << endl;
    }
    else { // Collect sequence data.
      trim(line, start,end);
      for(long i = start; i <= end; i++) {
	char c = std::tolower(line[i]);
	if(nucleotides_only) {
	  switch(c) {
	  case 'a': case 't': case 'g': case 'c': break;
	  default:
	    c = '~';
	  }
	}
	*P += c;
      }
    }
  }
  // Handle very last sequence.
  if(meta != "") {
    if(seq_cnt % arg->skip == arg->skip0) {
      cerr << "# P.length()=" << P->length() << endl;
      if(print) printf("> %s\n", meta.c_str());

      if(type == MAM) sa->MAM(*P, matches, min_len, print);
      else if(type == MUM) sa->MUM(*P, matches, min_len, print);
      else if(type == MEM) sa->MEM(*P, matches, min_len, print, num_threads);
      if(!print) sa->print_match(meta, matches, false); 
      if(rev_comp) {
	reverse_complement(*P, nucleotides_only);
	if(print) printf("> %s Reverse\n", meta.c_str());
	if(type == MAM) sa->MAM(*P, matches, min_len, print);
	else if(type == MUM) sa->MUM(*P, matches, min_len, print);
	else if(type == MEM) sa->MEM(*P, matches, min_len, print, num_threads);
	if(!print) sa->print_match(meta, matches, true); 
      }
    }
  }
  delete P;
  
  pthread_exit(NULL);
}


int main(int argc, char* argv[]) {
  // Collect parameters from the command line.
  while (1) {
    static struct option long_options[] = { 
      {"l", 1, 0, 0}, // 0
      {"mumreference", 0, 0, 0}, // 1
      {"b", 0, 0, 0}, // 2
      {"maxmatch", 0, 0, 0}, // 3
      {"mum", 0, 0, 0}, // 4
      {"mumcand", 0, 0, 0},  // 5
      {"F", 0, 0, 0}, // 6
      {"k", 1, 0, 0}, // 7
      {"threads", 1, 0, 0}, // 8
      {"n", 0, 0, 0}, // 9
      {"qthreads", 1, 0, 0}, // 10
      {0, 0, 0, 0} 
    };
    int longindex = -1;
    int c = getopt_long_only(argc, argv, "", long_options, &longindex);
    if(c == -1) break; // Done parsing flags.
    else if(c == '?') { // If the user entered junk, let him know. 
      cerr << "Invalid parameters." << endl;
      usage(argv[0]);
    }
    else {
      // Branch on long options.
      switch(longindex) { 
      case 0: min_len = atol(optarg); break;
      case 1: type = MAM; break;
      case 2: rev_comp = true;	break;
      case 3: type = MEM; break;
      case 4: type = MUM; break;
      case 5: type = MAM; break;
      case 6: _4column = true; break;
      case 7: K = atoi(optarg); break;
      case 8: num_threads = atoi(optarg); break;
      case 9: nucleotides_only = true; break;
      case 10: query_threads = atoi(optarg) ; break;
      default: break; 
      }
    }
  }
  if (argc - optind != 2) usage(argv[0]);

  if(K != 1 && type != MEM) { cerr << "-k option valid only for -maxmatch" << endl; exit(1); }
  if(num_threads <= 0) { cerr << "invalid number of threads specified" << endl; exit(1); }

  string ref_fasta = argv[optind]; 
  query_fasta = argv[optind+1];

  string ref;
  
  vector<string> refdescr; 
  vector<long> startpos;

  load_fasta(ref_fasta, ref, refdescr, startpos);

  // Automatically use 4 column format if there are multiple reference sequences.
  if(startpos.size() > 1) _4column = true;

  sa = new sparseSA(ref, refdescr, startpos, _4column, K);

  pthread_attr_t attr;  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  vector<query_arg> args(query_threads);
  vector<pthread_t> thread_ids(query_threads);  

  // Initialize additional thread data.
  for(int i = 0; i < query_threads; i++) { 
    args[i].skip = query_threads;
    args[i].skip0 = i;
  }

  // Create joinable threads to find MEMs.
  for(int i = 0; i < query_threads; i++) 
    pthread_create(&thread_ids[i], &attr, query_thread, (void *)&args[i]);

  // Wait for all threads to terminate.
  for(int i = 0; i < query_threads; i++) 
    pthread_join(thread_ids[i], NULL);    
  
  delete sa;
}


void usage(string prog) {
  cerr << "Usage: " << prog << " [options] <reference-file> <query-file>" << endl;
  cerr << "Implemented MUMmer v3 options:" << endl;
  cerr << "-mum           compute maximal matches that are unique in both sequences" << endl;
  cerr << "-mumreference  compute maximal matches that are unique in the reference-" << endl;
  cerr << "               sequence but not necessarily in the query-sequence (default)" << endl;
  cerr << "-mumcand       same as -mumreference" << endl;
  cerr << "-maxmatch      compute all maximal matches regardless of their uniqueness" << endl;
  cerr << "-l             set the minimum length of a match" << endl;
  cerr << "               if not set, the default value is 20" << endl;
  cerr << "-b             compute forward and reverse complement matches" << endl;
  cerr << "-F             force 4 column output format regardless of the number of" << endl;
  cerr << "               reference sequence inputs"  << endl;
  cerr << "-n             match only the characters a, c, g, or t" << endl;
  cerr << endl;
  cerr << "Additional options:" << endl;
  cerr << "-k             sampled suffix positions (one by default)" << endl;
  cerr << "-threads       number of threads to use for -maxmatch, only valid k > 1 " << endl;
  cerr << "-qthreads      number of threads to use for queries " << endl;
  cerr << endl;
  cerr << "Example usage:" << endl;
  cerr << endl;
  cerr << "./mummer -maxmatch -l 20 -b -n -k 3 -threads 3 query.fa ref.fa" << endl;
  cerr << "Find all maximal matches on forward and reverse strands" << endl;
  cerr << "of length 20 or greater, matching only a, c, t, or g." << endl;
  cerr << "Index every 3rd position in the ref.fa and use 3 threads to find MEMs." << endl;
  cerr << "Fastest method for one long query sequence." << endl;
  cerr << endl;
  cerr << "./mummer -maxmatch -l 20 -b -n -k 3 -qthreads 3 query.fa ref.fa" << endl;
  cerr << "Same as above, but now use a single thread for every query sequence in" << endl;
  cerr << "query.fa. Fastest for many small query sequences." << endl;
  
  exit(1);
}
