#include <cstdio>

#include <fstream>
#include <unordered_map>
#include <iostream> 
#include <string>
#include <sstream>

using namespace std;

struct bed {
  int key;

  string chr, id;
  int start, end;
  int reads;
  char strand;

  bed(): key(0), chr(""), id(""), start(0), end(0), reads(0), strand('.') {}

  bed(int key, string chr, string id,
      int start, int end, int reads, char strand) :
      key(key), chr(chr), id(id), start(start), end(end), 
      reads(reads), strand(strand) {}  
};

unordered_map<int, bed> rnas, dnas;

int main() {
  ifstream in ("GSM2396700_mESC_merged.ghits.pkbin.net.txt");
  ofstream rna ("GSM2396700_mESCRNAs.bed"); 
  ofstream dna ("GSM2396700_mESCDNAs.bed"); 
  ofstream gph ("GSM2396700_mESCgraph.bed"); 
 
  string line;
  
  string rna_chr, gene_id;
  int rna_start, rna_end, rna_reads;
  char rna_strand;

  string dna_chr;
  int dna_bin;
  
  double val;
  
  int rna_key = 1, dna_key = 1;

  if (in.is_open()) {
    while (getline(in, line)) {
      stringstream ss(line);
      ss >> rna_chr >> rna_start >> rna_end >> gene_id >> rna_reads >> rna_strand
         >> dna_chr >> dna_bin >> val;

      if (rnas.find(rna_start) == rnas.end()) {
        bed new_rna(rna_key++, rna_chr, gene_id, rna_start, rna_end, rna_reads, rna_strand);

        rnas[rna_start] = new_rna;
        rna << rna_chr << '\t' << rna_start << '\t' << rna_end 
          << '\t' << gene_id << '\t' << rna_reads << '\t' << rna_strand << '\n';
      }

      if (dnas.find(dna_bin) == dnas.end()) {
        bed new_dna(dna_key++, dna_chr, " ", dna_bin, dna_bin + 999, 0, '+');

        dnas[dna_bin] = new_dna;

        dna << dna_chr << '\t' << dna_bin << '\t' << dna_bin + 999 << '\n';
      }   

      gph << rnas[rna_start].key << ' ' << dnas[dna_bin].key << ' ' << val << endl;


    } 
  }

  return 0;
}
