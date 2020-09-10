#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
#include <set>
#include <cstdlib>
#include <math.h>
#include <mutex>
#include <future>
#include <chrono>
#include <atomic>
#include <algorithm>
//#include "bamtools/include/api/BamReader.h"
#include "bamtools/src/api/BamReader.h"
#include <boost/math/special_functions/binomial.hpp>
#include <boost/dynamic_bitset.hpp>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

//#include "vcflib/include/Variant.h"

using namespace std;
using namespace BamTools;
//using namespace vcflib;


mutex mtx;

const string STDERR_INFO_SEPARATOR = "\n***************************\n";
static const unsigned int ALPHABET_SIZE = 4;
static constexpr char ALPHABET[ALPHABET_SIZE] = {'A', 'C', 'G', 'T'};
string RNA_SEQ_BAM_FILE, VCF_FILE, GTF_FILE, GENE_NAMES_TABLE_FILE, RNAMAP_REP_FILE, VCF_FILTER, EXPR_FILE, GENE_EXPR_FILE;
unsigned int THREAD_NUM = 1;
bool VERBOSE = false, USE_ONLY_FIRST_MATE = false, USE_ONLY_HI_QUALITY_ALIGNMENTS = false, FIX_VCF = false, PHASE_ANALYSIS = false;
double HQ_THRESHOLD = 0.05, MATCH_FRAC = 0.9;

void prefetch_variants(string&, map<string, vector<bool> >&);

class exon
{
	private:
	unsigned int start, end;
	public:
	exon(const unsigned int, const unsigned int);
	const string report() const;
	const unsigned int get_start() const;
	const unsigned int get_end() const;
};

class transcript
{
	private:
	vector<exon> EXONS;	
	string tr_id, gene_id;
	struct pos_pointer{unsigned int chr_ref; unsigned int pos; unsigned int T;};
	vector<pos_pointer> POS_VAR_V;
	double expr, expr_w;
	mutable vector<double> V_prob[2];
	mutable unsigned int count_var;
	mutable bool aspec_flag;
	public:
	enum {ALLELE_SPECIFIC, NOT_ALLELE_SPECIFIC};
	static constexpr double NO_EXPRESSION = 0;
	transcript(const string*, const string*);
	void add_exon(const unsigned int, const unsigned int);
	const string report() const;
	void add_pos(unsigned int, unsigned int, unsigned int);
	const vector<exon>* get_exons_ptr() const;
	const string * get_id() const;
	const unsigned int get_pos_ptrv_size() const;
	const unsigned int get_pos_ptrv_chr_ref_at(unsigned int) const;
	const unsigned int get_pos_ptrv_pos_at(unsigned int) const;
	const double get_expr() const;
//	const double get_expr_weight() const;
	const string get_gene_id() const;
	void set_expr(double);
	void set_expr_weight(double);
	void add_var_count(bool) const;
	const bool get_aspec_flag() const;
	const unsigned int get_count_var() const;
	void add_pos_prob(double*) const;
	const double get_prod_prob(unsigned int) const;
	const string all_prob_report(unsigned int) const;
	const string snp_count_report() const;
};

class gene
{
	private:
	struct snp_counts {unsigned int big_c, small_c;};
	map<string, transcript> TRANSCRIPTS;
	vector<snp_counts> SNP_C;
	const RefVector *refv;
	string gene_id, chr;
	char strand;
	int chr_ref;
	void set_chr_ref();
	double tot_expr;
	double expression;
	double min_pv, max_pv, tot_pv;
	mutable unsigned int omo_count, etero_count;
	public:
	static constexpr int INVALID_CHR = -1;
	gene(const string *, const string*, const char, const RefVector*);	
	const map<string, transcript> * get_transcripts_ptr() const;
	map<string, transcript> * get_transcripts_ptr() ;
	map<string, transcript>::iterator add_transcript_or_exon(const string*, const unsigned int, const unsigned int);
	const string report() const;
	const int get_chr_ref() const;
	const bool has_transcript_with_var() const;
	const string id() const;
	void set_expression_weight();
	void add_etero_var_to_count() const;
	void add_omo_var_to_count() const;
	const string report_var_counts() const;
	void add_snp_c(unsigned int, unsigned int);
	void add_snp_c(unsigned int, unsigned int, short int);
	const unsigned int get_big_c() const;
	const unsigned int get_small_c() const;
	const double gtest_1() const;
	void gtest_2();
	const double get_expr() const;
	void set_expr(double);
	const double get_tot_pv() const;
	const double get_min_pv() const;
	const double get_max_pv() const;
};

class gtf_file
{
	private:
	map<string, map<string, gene> > GENES;	
	map<string, string> TR_TO_GENE;
	map<string, map<string, transcript>::iterator> ID_TO_TR_MAP;
	map<string, double> EXPRESSION, GENE_EXPRESSION;
	const RefVector *refv;
	void clean_ids(string *);
	void read_tr_to_gene_file(const string*);
	void read_expression(const string &, map<string, double> &);
	void set_expression();
	void set_expression_gene();
	void set_expression_weight();
	public:
	gtf_file(const string*, const string *, const RefVector*);
	const string* get_gene_name(const string*) const;
	void report() const;
	const map<string, map<string, gene> >* get_genes_ptr() const;
	map<string, map<string, gene> >* get_genes_ptr();
};

class variant
{
	private:
	enum {A,C,G,T,N};
	bool PHASED;
	short VAR_CHR;
//	static const unsigned int ALPHABET_SIZE = 4;
	unsigned int type, zygosity;
	map<string, string> GT;
	string chr, id, ref, alt, filter, info, gtfields, gt;
	const RefVector *refv;
	unsigned int pos, internal_id;
	int chr_ref;
	double qual;
	void set_type();
	void set_chr_ref();
	void set_zygosity();
	void set_tr_hyp_prob_table();
	const unsigned int get_char_int(const string *) const;
	public:
	enum {SNP, INS, DEL, OTHER, FILTERED};
	enum {REF_OMO_Z, ALT_OMO_Z, ETERO_Z, DOUBLE_ETERO_Z, UNKNOWN_Z, UNPHASED};
	static constexpr int INVALID_CHR = -1;
	static constexpr unsigned int type_num = 5, zyg_types = 6;
	variant(const string *, const RefVector *, const unsigned int);
	unsigned int get_type() const;
	unsigned int get_pos() const;
	unsigned int get_zygosity() const;
	const int get_chr_ref() const;
	const string get_chr() const;
	const string get_ref() const;
	const string get_alt() const;
	const unsigned int get_ref_int() const;
	const unsigned int get_alt_int() const;
	const unsigned int get_internal_id() const;
	const short int get_phase() const;
	const string* zygosity_str() const;
	void gt_report();
	void report();
};

class position
{
	private:
	enum {A,C,G,T};
	enum Pos_type {ETERO_ZYG_EXPR, ALLELE_SPECIFIC_EXPR, RNA_EDIT, OMOZYG_EXPR, DUBIOUS_EXPR, UNKNOWN_EXPR};
	enum {AS_1, AS_0, ET};
//	struct tcomb{double b[2];};
//	map<unsigned int, vector<vector<tcomb> > >  *VCOMB;
	map<unsigned int, vector<boost::dynamic_bitset<> > > *BCOMB;
//	vector<double> PV;
	Pos_type p_type;
	bool MULTIPLE_GENE, MULTIPLE_GENE_SET; 
	mutable bool PERF_OMO, PERF_OMO_SET;
	atomic<unsigned int> N_ATOM[ALPHABET_SIZE];
	vector<const transcript *> VPTR;
	const variant *var_ptr;
	void var_diagnostic(unsigned int);
	void build_b_comb(map<unsigned int, vector<boost::dynamic_bitset<> > >::iterator);
	void set_multiple_gene();
	void tr_expr_weights(vector<double> &) const;
	void get_counts(unsigned int, unsigned int, double &, double &, bool &) const;
//	bool var_flag;
	public:
	enum {OMO_Z, ETERO_Z, UNKNOWN_Z};
	position(char);
//	~position();
	position(const position &obj);
	void assign_nuc_atomic(char);
	const string report_s_atomic() const;
	const string report_s_transcripts() const;
	void add_transcript(const transcript *);
	void set_var(const variant *);
	const variant* get_var() const;
	const bool tr_marker() const;
	const string* zygosity_str() const;
	const int get_zygosity() const;
	const string* var_diagnostic_s(unsigned int); 
	const Pos_type var_diagnostic_p(unsigned int);
	const double marker_p(unsigned int, unsigned int, double*) const;
	void compute_probability_binom(vector<double> &, vector<double> &, double *, unsigned int, unsigned int);
	double compute_probability_gtest(vector<double> &, vector<double> &, vector<double> &, unsigned int &, unsigned int, unsigned int) const;
//	void build_v_comb();
	const string report_b_comb();
	const string report_b_comb_at(unsigned int) const;
	const unsigned int get_count_at(const unsigned int) const;
	const bool multiple_gene();
	const unsigned int trv_size() const;
	const bool perf_omo() const;
	void add_var_to_tr(unsigned int) const;
	void pass_prob_values_to_tr(vector<double> &) const;
	const unsigned int get_ref_count() const;
	const unsigned int get_alt_count() const;
	const short int get_phase() const;
};

class gmap
{
	private:
	vector<map<unsigned int, position > > GENOME;
	public:
	gmap(unsigned int);
	gmap();
	void init(unsigned int);
	void assign(unsigned int, unsigned int, char);
	void report(const RefVector*);
	void report_detailed(const RefVector*);
	const position* get_position_ptr(unsigned int, unsigned int) const;
	position* get_position_ptr(unsigned int, unsigned int);
};

class rnaseq
{
	private: 
	static const unsigned int BUF_ALIGNMENTS = 75000;
	map<string, vector<bool> >& var_map;
	RefVector refv; 
	string rna_seq_bam_file;
	gmap GMAP;
	BamReader bam;
	vector<bool> busy_vector;
	vector<future<bool> > FUTV;
	int busy_check();
	void open_bam();
	void process_alignment(BamAlignment *, gmap *);
	void build_gmap_from_bam();
	void build_gmap_from_bam_multi();
        bool build_gmap_from_bam_multi_sub(vector<BamAlignment>*, gmap*);
	bool alignment_quality_check(BamAlignment *);
	public:
	rnaseq(const string*, map<string, vector<bool> >&);
	const RefVector* refv_ptr ();
	const position* get_position_ptr(unsigned int, unsigned int) const;
	position* get_position_ptr(unsigned int, unsigned int);
	const string get_file_name();
};


class vcf_file
{
	private:
	vector<vector<variant> > VV_MAP;
	vector<string> header;
	vector<unsigned int> type_count;
	ifstream in;
	public:
	vcf_file(const string *, const RefVector*);
	void type_report();
	unsigned int VV_MAP_size();
	unsigned int V_MAP_size(unsigned int);
	const variant* get_variant_at(unsigned int, unsigned int) const;
//	const vector<vector<variant> >* get_VV_MAP_ptr();
};

class analysis
{
	private:
	rnaseq RNA;
	vcf_file VCF;
	gtf_file GTF;
	void assign_transcript_to_pos();
	void add_snp_counts_to_genes();
	void flag_var_positions();
	map<gene*, set<position*> > GENE_POS;
	public:	
	analysis(const string*, const string*, const string*, const string*, map<string, vector<bool> >& var_map);
	void simple_snp_report(const string );
	void gene_level_analysis(const string) const;
	void simple_transcript_report(const string );
	void gene_report(const string);
	void marker_var_transcript_report(const string );
	void tr_final_prob_report(const string) const;
	void tr_and_gene_var_report(const string, const string) const;
//	const string* var_diagnostic(unsigned int, unsigned int) const; 
	const string position_full_report(const variant *, position*);
	const string position_small_report(const variant*, position*);
	const string asonly_transcripts() const;
};


void command_line_parser(int, char **);
void display_help();

