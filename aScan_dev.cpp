#include "aScan_dev.h"

int main(int argc, char **argv)
{
	if(argc == 1)
		display_help();

	command_line_parser(argc, argv);
	map<string, vector<bool> > var_map;

//	rnaseq R(&RNA_SEQ_BAM_FILE);
	prefetch_variants(VCF_FILE, var_map);

	analysis A(&RNA_SEQ_BAM_FILE, &VCF_FILE, &GTF_FILE, &GENE_NAMES_TABLE_FILE, var_map);

	return EXIT_SUCCESS;
}

void prefetch_variants(string& VCF_FILE, map<string, vector<bool> >& var_map)
{
	ifstream in(VCF_FILE.c_str());	

	if(!in)
        {
                cerr << "\nCan't find VCF file: " << VCF_FILE << endl;
                exit(EXIT_FAILURE);
        }

	cerr << "\nPrefetching VCF file...";

	string line;
	
	while(getline(in,line))
	{
		if(line.empty())
			continue;

		if(line[0] == '#')
			continue;	

		istringstream str(line);
		string chr;
		unsigned int pos;

		str >> chr >> pos;

		pos--;

		if(FIX_VCF)
			chr.insert(0, "chr");

		map<string, vector<bool> >::iterator mi = var_map.find(chr);

		if(mi == var_map.end())
			mi = var_map.emplace(make_pair(chr, vector<bool>())).first;

		if(mi->second.size() < pos + 1)
			mi->second.resize(pos + 1, false);

		mi->second.at(pos) = true;
	}

	in.close();
	
	return;
}

void display_help()
{
	cerr << "\nSyntax: aScan --rna rnaseq_bamfile --vcf vcf_file --gtf gtf_file --only_first_mate [-nt transcript/gene_correspondence_file] [-p thread_num] [--filter filter_key_word]\n" << endl;
	
	exit(EXIT_SUCCESS);
}

void command_line_parser(int argc, char **argv)
{
	for(int i = 1; i < argc; i++)
        {
                string buf = argv[i];

                if(buf == "--rna" || buf == "-r")
                {
                        if(i < argc - 1)
                                RNA_SEQ_BAM_FILE = argv[++i];

                        continue;
                }
		else if(buf == "--vcf" || buf == "-V")
                {
                        if(i < argc - 1)
                                VCF_FILE = argv[++i];

                        continue;
                }	
		else if(buf == "-h" || buf == "--help")
		{
			display_help();
		}
		else if(buf == "--gtf" || buf == "-g")
                {
                        if(i < argc - 1)
                                GTF_FILE = argv[++i];

                        continue;
                }
		else if(buf == "--texpr" || buf == "-te")
                {
                        if(i < argc - 1)
                                EXPR_FILE = argv[++i];

                        continue;
                }
		else if(buf == "--gexpr" || buf == "-ge")
                {
                        if(i < argc - 1)
                                GENE_EXPR_FILE = argv[++i];

                        continue;
                }
		else if(buf == "--nametable" || buf == "-nt")
                {
                        if(i < argc - 1)
                                GENE_NAMES_TABLE_FILE = argv[++i];

                        continue;
                }
		else if(buf == "-p")
		{
			if( i < argc - 1)
				THREAD_NUM = atoi(argv[++i]);

			continue;
				
		}
		else if(buf == "--filter" || buf == "-f")
		{
			if( i < argc - 1)
                                VCF_FILTER = argv[++i];

			continue;
		}
		else if(buf == "-rep")
                {
                        if( i < argc - 1)
                                RNAMAP_REP_FILE = argv[++i];

                        continue;

                }
		else if(buf == "-v")
		{
			VERBOSE = true;
			continue;
		}	
		else if(buf == "--phase" || buf == "-phase")
                {
                        PHASE_ANALYSIS = true;
                        continue;
                }
		else if(buf == "-hq")
		{
			USE_ONLY_HI_QUALITY_ALIGNMENTS = true;
			continue;
		}
		else if(buf == "--only_first_mate")
		{
			USE_ONLY_FIRST_MATE = true;
			continue;
		}
		else if(buf == "--fix_vcf")
                {
                        FIX_VCF = true;
                        continue;
                }
	
         /*       else if(buf == "-b")
                {
                        if(i < argc - 1)
                                bedfile = argv[++i];

                        continue;
                }
                else if(buf == "-shift")
                {
                        if(i < argc - 1)
                                SHIFT = atoi(argv[++i]);

                        if(SHIFT > WARN_MAX_SHIFT)
                                cerr << "\nWarning: the set shift seems way too big! shift = " << SHIFT << endl;

                        continue;
                }
                else if(buf == "-zoom")
                {
                        if(i < argc - 1)
                                zoomfile = argv[++i];

                        unlink(zoomfile.c_str());

                        continue;
                }
*/
		else
                {
                        cerr << "\n\nUnknown option: " << buf << endl << endl;
                        exit(EXIT_FAILURE);
                }
	}

	return;
}

// POSITION CLASS

position::position(char C):MULTIPLE_GENE(false), MULTIPLE_GENE_SET(false), PERF_OMO(false), PERF_OMO_SET(false)
{
//	static atomic<unsigned int> NATOM[4];
	var_ptr = NULL;
	p_type = UNKNOWN_EXPR;

	static map<unsigned int, vector<boost::dynamic_bitset<> > > bcomb;
//	static map<unsigned int, vector<vector<tcomb> > > vcomb;

	BCOMB = &bcomb;
//	VCOMB = &vcomb;
//	COLL_COMB = &coll_comb;
	
	for(unsigned int i = 0; i < ALPHABET_SIZE; i++)
	{
//		N[i] = 0;
		N_ATOM[i].store(0);
	}
	
//	assign_nuc(C);	
	assign_nuc_atomic(C);

	return;
}

const short int position::get_phase() const
{
	return var_ptr->get_phase();
}

const unsigned int position::get_ref_count() const
{
	return N_ATOM[var_ptr->get_ref_int()];
}

const unsigned int position::get_alt_count() const
{
        return N_ATOM[var_ptr->get_alt_int()];
}

const unsigned int position::trv_size() const
{
	return VPTR.size();
}
/*
position::~position()
{
	return;
}
*/
position::position(const position &obj)
{
	for(unsigned int i = 0; i < ALPHABET_SIZE; i++)
	{
//		N[i] = obj.N[i];
		N_ATOM[i].store(obj.N_ATOM[i].load());
	}
	
	VPTR = obj.VPTR;
	var_ptr = obj.var_ptr;
	p_type = obj.p_type;
//	VCOMB = obj.VCOMB;
	BCOMB = obj.BCOMB;
	MULTIPLE_GENE = obj.MULTIPLE_GENE;
	MULTIPLE_GENE_SET = obj.MULTIPLE_GENE_SET;
	PERF_OMO = obj.PERF_OMO;
	PERF_OMO_SET = obj.PERF_OMO_SET;
//	PV = obj.PV;

	return;
}
/*
position::position(const unsigned int *arr)
{
	for(unsigned int i = 0; i < ALPHABET_SIZE; i++)
		N[i] = arr[i];

	return;
}

void position::add(const unsigned int *arr)
{
	for(unsigned int i = 0; i < ALPHABET_SIZE; i++)
                N[i] += arr[i];

        return;
}
*/
/*const unsigned int * position::n_arr()
{
	return N;
}*/
/*
inline void position::assign_nuc(char c)
{
	switch (c)
	{
		case 'A':
			N[A]++;
			break;
		case 'C':	
			N[C]++;
			break;
		case 'G':
			N[G]++;
			break;
		case 'T':
			N[T]++;
			break;
		default:
			break;
	}

	return;
}
*/
inline void position::assign_nuc_atomic(char c)
{
        switch (c)
        {
                case 'A':
                        N_ATOM[A]++;
                        break;
                case 'C':
                        N_ATOM[C]++;
                        break;
                case 'G':
                        N_ATOM[G]++;
                        break;
                case 'T':
                        N_ATOM[T]++;
                        break;
                default:
                        break;
        }

        return;
}

/*string position::report_s()
{
	ostringstream out;

	for(unsigned int i = 0; i < ALPHABET_SIZE; i++)
		out << N[i] << '\t';

	return out.str();
		
}*/

const string position::report_s_atomic() const
{
        ostringstream out;

        for(unsigned int i = 0; i < ALPHABET_SIZE; i++)
                out << '\t' << N_ATOM[i].load();

        return out.str();

}

const string position::report_b_comb() 
{
//	build_v_comb();
	ostringstream out;

/*	unsigned int nxcount = 0;

	for(unsigned int i = 0; i < VPTR.size(); i++)
		if(VPTR.at(i)->get_expr() <= transcript::NO_EXPRESSION)
			nxcount++;
*/
//	cerr << endl << VPTR.size() << '\t' << nxcount;

	auto bi = BCOMB->find(VPTR.size());
//	auto vi = VCOMB->find(VPTR.size());
	
	out << '\t' << bi->second.size() << '\t' << VPTR.size() << '\t' << pow(2, VPTR.size());

	if(bi->second.size() != pow(2, VPTR.size()))
		out << "ERROR!";

//	out << "\t---\t" << vi->second.size() << '\t' << VPTR.size() << '\t' << pow(2, VPTR.size());

	for(unsigned int i = 0; i < bi->second.size(); i++)
	{
		out << '\t';

		for(unsigned int j = 0; j < bi->second.at(i).size(); j++)
			out << bi->second[i][j] << ',';
	}

	return out.str();
}

void position::add_var_to_tr(unsigned int p) const
{
	auto bi = BCOMB->find(VPTR.size());

	for(unsigned int i = 0; i < bi->second.at(p).size(); i++)	
		VPTR.at(i)->add_var_count(bi->second[p][i]);

	return;
}

const string position::report_b_comb_at(unsigned int p) const
{
	ostringstream out;

	auto bi = BCOMB->find(VPTR.size());

	if(bi == BCOMB->end())
		return out.str();

	for(unsigned int i = 0; i < bi->second.at(p).size(); i++)
		out << bi->second[p][i] << ',';

	return out.str();
}

void position::add_transcript(const transcript *TR)
{
	VPTR.push_back(TR);	
}

const string position::report_s_transcripts() const
{
	ostringstream out;

	for(unsigned int i = 0; i < VPTR.size(); i++)
		out << *VPTR.at(i)->get_id() << ',';

	return out.str();	
}

void position::set_var(const variant *V)
{
	var_ptr = V;
	return;
}

const variant * position::get_var() const
{
	return var_ptr;
}

const bool position::tr_marker() const
{
	if(VPTR.size() == 1)
		return true;
	else
		return false;
}

const string* position::zygosity_str() const 
{
	static const vector<string> ZYG = {"RNA_OMOZYG", "RNA_ETEROZYG", "RNA_???"};
	int Z = -1;
	vector<double> frac(ALPHABET_SIZE, 0); 
	double sum = 0;

	for(unsigned int i = 0; i < ALPHABET_SIZE; i++)
	{
		frac[i] = N_ATOM[i].load();	
		sum += N_ATOM[i].load();
	}

	for(unsigned int i = 0; i < ALPHABET_SIZE; i++)
	{
		frac[i] /= sum;

		if(frac[i] != 0)
			Z++;
	}

	switch (Z)
	{
		case 0:
		return &ZYG[Z];
		break;
		case 1:
		return &ZYG[Z];
		break;
		default:
		return &ZYG[2];
		break;
	}

}

const int position::get_zygosity() const
{
        int Z = -1;
        vector<double> frac(ALPHABET_SIZE, 0);
        double sum = 0;

        for(unsigned int i = 0; i < ALPHABET_SIZE; i++)
        {
                frac[i] = N_ATOM[i].load();
                sum += N_ATOM[i].load();
        }

        for(unsigned int i = 0; i < ALPHABET_SIZE; i++)
        {
                frac[i] /= sum;

                if(frac[i] != 0)
                        Z++;
        }

        switch (Z)
        {
                case 0:
                return OMO_Z;
                break;
                case 1:
                return ETERO_Z;
                break;
                default:
                return UNKNOWN_Z;
                break;
        }

}

void position::var_diagnostic(unsigned int vcf_z)
{
	unsigned int rna_z = get_zygosity();
	
	if(vcf_z == variant::ETERO_Z && rna_z == position::ETERO_Z)
		p_type = ETERO_ZYG_EXPR;
        else if(vcf_z == variant::ETERO_Z && rna_z == position::OMO_Z)
		p_type = ALLELE_SPECIFIC_EXPR;
        else if(vcf_z == variant::ALT_OMO_Z && rna_z == position::ETERO_Z)
		p_type = RNA_EDIT;
        else if(vcf_z == variant::ALT_OMO_Z && rna_z == position::OMO_Z)
		p_type = OMOZYG_EXPR;
        else
                p_type = DUBIOUS_EXPR;	

	return;
}

const string* position::var_diagnostic_s(unsigned int vcf_z)  
{
	static const vector<string> D = {"ETERO_ZYG_EXPR", "ALLELE_SPECIFIC_EXPR", "RNA_EDIT", "OMOZYG_EXPR", "DUBIOUS_EXPR", "UNKNOWN_EXPR"};

	if(p_type == UNKNOWN_EXPR)
		var_diagnostic(vcf_z);

	return &D[p_type];
}

const position::Pos_type position::var_diagnostic_p(unsigned int vcf_z)
{
	if(p_type == UNKNOWN_EXPR)
                var_diagnostic(vcf_z);

	return p_type;
}

void position::get_counts(unsigned int p1, unsigned int p2, double &x, double &n, bool &v1p) const
{
	if(N_ATOM[p1] >= N_ATOM[p2])
        {
                x = N_ATOM[p1].load();
                n = x + N_ATOM[p2].load();
                v1p = false;
        }
        else
        {
                x = N_ATOM[p2].load();
                n = x + N_ATOM[p1].load();
                v1p = true;
        }

	if(n == x)
		PERF_OMO = true;

	PERF_OMO_SET = true;

	return;
}

double position::compute_probability_gtest(vector<double> &V_V1, vector<double> &V, vector<double> &chi_raw, unsigned int &bgv, unsigned int p1, unsigned int p2) const
{
	double x, n;
        bool v1p;

	get_counts(p1, p2, x, n, v1p);
	auto bi = BCOMB->find(VPTR.size());  //BCOMB of VPTR.size() exists since compute_probability_binom is called before compute_probability_gtest

//	if(x == n)
	if(perf_omo())
	{
		double p = 1.0 - gsl_ran_binomial_pdf(n, 0.5, n);  //IF ONE POSITION IS 0 THAN SWITCH TO BINOMIAL TEST WITH P = 0.5
		bgv = 0;
	
		V.push_back(p);   //ONLY THE FIRST BOOL VECTOR (0,0,0...0) HAS PROB. P. 
		chi_raw.push_back(0); // AND CHI 0 

	}
	
	{	
		unsigned int i;

		if(perf_omo())
			i = 1;
		else
			i = 0;

		for(; i < bi->second.size(); i++)	
		{
			double a1, a2;

			a1 = V_V1.at(i) * n;	
			a2 = (1.0 - V_V1.at(i)) * n;

			double tmp;
		
			if(perf_omo())
				tmp = (x * log(x / a1));
			else

				tmp = (x * log(x / a1)) + ((n -x) * log((n-x) / a2 ));

//			V.push_back((x * log(x / a1)) + ((n -x) * log((n-x) / a2 )));
	
			tmp *= 2;

			chi_raw.push_back(tmp);
			V.push_back(gsl_cdf_chisq_Q(tmp, 1));
		}
	}

	double max = V.at(0);
	bgv = 0;

	for(unsigned int i = 1; i < V.size(); i++)
	{
		if(V.at(i) > max)
		{
			max = V.at(i);
			bgv = i;
		}
	}

	return V.at(bgv);
//	return gsl_cdf_chisq_Q(2 * min, 1);
}

void position::compute_probability_binom(vector<double> &PV, vector<double> &V_V1, double *Att, unsigned int p1, unsigned int p2) 
{
	static constexpr double PAS[2] = {0.99, 0.5};
	vector<double> weights;

	if(p1 >= ALPHABET_SIZE || p2 >= ALPHABET_SIZE)
        {
                cerr << "\n\nFATAL ERROR!!";
                exit(EXIT_FAILURE);
        }

	if(VPTR.empty())
                return;

	auto bi = BCOMB->find(VPTR.size());

	if(bi == BCOMB->end())
	{	
		bi = BCOMB->insert(make_pair(VPTR.size(), vector<boost::dynamic_bitset<> >())).first;
		build_b_comb(bi);
	}

	double x, n;
	bool v1p;

	get_counts(p1, p2, x, n, v1p);

	tr_expr_weights(weights);

	for(unsigned int i = 0; i < bi->second.size(); i++)
	{
		double V1 = 0;		

		for(unsigned int j = 0; j < bi->second.at(i).size(); j++)
		{
			double frac = weights.at(j);
			V1 += frac * PAS[bi->second.at(i)[j]];	
		}

		PV.push_back(gsl_ran_binomial_pdf(x, V1, n));
		V_V1.push_back(V1);
	}

	double max = 0; 
	unsigned int vpos = 0;	

	for(unsigned int i = 0; i < PV.size(); i++)
	{
		if(max < PV.at(i))
		{
			max = PV.at(i);
			vpos = i;	
		}
	}

	Att[v1p] = V_V1.at(vpos) * n;
	Att[!v1p] = (1 - V_V1.at(vpos)) * n;

	return;
}

void position::tr_expr_weights(vector<double> &W) const
{
	double sum = 0; 

	for(auto vi = VPTR.begin(); vi != VPTR.end(); ++vi)
		sum += (*vi)->get_expr();

	for(auto vi = VPTR.begin(); vi != VPTR.end(); ++vi)
		W.push_back((*vi)->get_expr() / sum);

	return;	
}

const unsigned int position::get_count_at(const unsigned int p) const
{
	if(p >= ALPHABET_SIZE)
	{
		cerr << "\n\nFATAL ERROR!!";
                exit(EXIT_FAILURE);
	}

	return N_ATOM[p].load();
}

const double position::marker_p(unsigned int p1, unsigned int p2, double *rec_p) const
{
	if(p1 >= ALPHABET_SIZE || p2 >= ALPHABET_SIZE)
	{
		cerr << "\n\nFATAL ERROR!!";
		exit(EXIT_FAILURE);
	}

	long double boost_bin_coeff, prob, x, n, rec;	

	if(N_ATOM[p1] >= N_ATOM[p2])
	{
		x = N_ATOM[p1].load();
		n = x + N_ATOM[p2].load();
	}	
	else
	{
		x = N_ATOM[p2].load();
		n = x + N_ATOM[p1].load();
	}


//	bin_coeff = tgamma(n + 1) / (tgamma(x + 1) * tgamma(n - x + 1));

//	cerr << endl << "bin_coeff " << n << ' ' << x << " = " << bin_coeff << '\t' ;

	typedef boost::math::policies::policy<boost::math::policies::overflow_error<boost::math::policies::ignore_error> > my_policy;

	boost_bin_coeff = boost::math::binomial_coefficient<long double> (n, x, my_policy());

//	cerr << boost_bin_coeff;
	
	prob = boost_bin_coeff * pow(0.5, x) * pow(0.5, n - x);

	rec = boost_bin_coeff * pow(0.99, x) * pow(0.01, n - x);

	*rec_p = rec;

//	cerr << endl << rec << '\t' << *rec_p << endl;

	return prob;
}
/*
void position::build_v_comb()   //DA ELIMINARE
{
	if(VPTR.empty())
		return;

	auto bi = BCOMB->find(VPTR.size());

	if(bi != BCOMB->end())
		return;

//	map<unsigned int, vector<vector<tcomb> > >::iterator vi;
	bi = BCOMB->insert(make_pair(VPTR.size(), vector<boost::dynamic_bitset<> >())).first;
//	vi = VCOMB->insert(make_pair(VPTR.size(), vector<vector<tcomb> >())).first; 
			
	build_b_comb(bi);

	vi->second.resize(bi->second.size());
	
	for(unsigned int i = 0; i < bi->second.size(); i++)
	{
		for(unsigned int t = 0; t < VPTR.size(); t++)
		{
			tcomb COMB;
			
			if(bi->second.at(i)[t])
			{
				COMB.b[0] = PAS[AS_0];	
				COMB.b[1] = PAS[AS_1];
			}
			else
			{
				for(unsigned int k = 0; k < 2; k++)
                                        COMB.b[k] = PAS[ET];
			}

			vi->second.at(t).push_back(COMB);
		}		
	}

	auto ci = COLL_COMB->insert(make_pair(VPTR.size(), vector<tcomb>())).first;

	ci->second.resize(bi->second.size());

	for(unsigned int i = 0; i < bi->second.size(); i++)
	{
		ci->second.at(i).b[0] = ci->second.at(i).b[1] = 1;

		for(unsigned int j = 0; j < vi->second.at(i).size(); j++)
		{
			ci->second.at(i).b[0] *= vi->second.at(i).at(j).b[0];
			ci->second.at(i).b[1] *= vi->second.at(i).at(j).b[1];
		}
	}

	return;
}
*/
//void position::build_b_comb(vector<bool> &VB, unsigned int l)
void position::build_b_comb(map<unsigned int, vector<boost::dynamic_bitset<> > >::iterator bi)
{
	unsigned long int p = pow(2,bi->first);

	for(unsigned long int i = 0; i < p; i++) 	
		bi->second.push_back(boost::dynamic_bitset<> (bi->first,i));

	return;
}

const bool position::multiple_gene()
{
	if(!MULTIPLE_GENE_SET)
		set_multiple_gene();
	
	return MULTIPLE_GENE;
}

const bool position::perf_omo() const
{
	if(!PERF_OMO_SET)
		cerr << "\nWarning, perfect omozygosity requested but never set\n";

	return PERF_OMO;
}

void position::set_multiple_gene()
{
	if(VPTR.empty())
		return;

	const string gid = VPTR.at(0)->get_gene_id();	

	for(unsigned int i = 1; i < VPTR.size(); i++)
	{
		if(gid != VPTR.at(i)->get_gene_id())
		{
			MULTIPLE_GENE = true;
			break;
		}
	}

	MULTIPLE_GENE_SET = true;

	return;
}

void position::pass_prob_values_to_tr(vector<double> &CHI_RAW) const
{
/*
	if(!perf_omo())
	{
		for(auto vi = G_V.cbegin(); vi != G_V.cend(); ++vi)
			P_V.push_back(gsl_cdf_chisq_Q(*vi * 2, 1));
	}
	else
	{
		for(auto vi = G_V.cbegin(); vi != G_V.cend(); ++vi)
			P_V.push_back(*vi);
	}
*/
	auto bi = BCOMB->find(VPTR.size()); //BCOMB OF APPROPRIATE SIZE ALREADY EXISTENT SINCE THIS RUNS AFTER COMPUTE_PROBABILITY	

	for(unsigned int i = 0; i < VPTR.size(); i++)
	{
		double P[2] = {0,0};

		for(unsigned int bv = 0; bv < bi->second.size(); bv++)
				P[bi->second.at(bv)[i]] += CHI_RAW.at(bv);			

/*		for(unsigned int j = 0; j < 2; j++)
			if(P[j] > 1)
				P[j] = 1;   //AVOID PROBABLITIES GREATER THAN 1
*/
		for(unsigned int j = 0; j < 2; j++)
			P[j] = gsl_cdf_chisq_Q(P[j], bi->second.size() / 2);


		VPTR.at(i)->add_pos_prob(P);
	}

	return;
}

// END OF POSITION CLASS

// RNASEQ CLASS

rnaseq::rnaseq(const string *File, map<string, vector<bool> >& Var_map):var_map(Var_map)
{
	rna_seq_bam_file = *File;

	open_bam();

	GMAP.init(refv.size());

	FUTV.resize(THREAD_NUM);
	busy_vector.resize(THREAD_NUM, false);
	
	if(THREAD_NUM == 1)
		build_gmap_from_bam();
	else
		build_gmap_from_bam_multi();


	if(VERBOSE)
		GMAP.report(&refv);

	GMAP.report_detailed(&refv);

	bam.Close();

	return;
}

const RefVector* rnaseq::refv_ptr()
{
	return &refv;
}

const string rnaseq::get_file_name()
{
	return rna_seq_bam_file;
}

const position* rnaseq::get_position_ptr(unsigned int ref_sequence, unsigned int pos) const
{
	return GMAP.get_position_ptr(ref_sequence, pos);
}

position* rnaseq::get_position_ptr(unsigned int ref_sequence, unsigned int pos) 
{
        return GMAP.get_position_ptr(ref_sequence, pos);
}

void rnaseq::open_bam()
{
	bam.Open(rna_seq_bam_file.c_str());

	if(!bam.IsOpen())
        {
                cerr << "\nCan't find file: " << rna_seq_bam_file << "..." << endl;
                exit(1);
        }

	refv = bam.GetReferenceData();

	return;
}

void rnaseq::build_gmap_from_bam()
{
/*	BamReader bam;
	
	bam.Open(rna_seq_bam_file.c_str());

	if(!bam.IsOpen())
        {
        	cerr << "\nCan't find file: " << rna_seq_bam_file << "..." << endl;
                exit(1);
        }

	refv = bam.GetReferenceData();
*/
	cerr << "\nReading " << rna_seq_bam_file << " and building genomic map...\n";

//	cerr << endl << "var_map.size() = " << var_map.size() << endl;


//	gmap_init();

	BamAlignment ba;

	while(bam.GetNextAlignment(ba))
		process_alignment(&ba, &GMAP);
	
//	bam.Close();

	return;
}
/*
void rnaseq::process_alignment(BamAlignment *ba)
{
	if(ba->IsPaired())
	{
		if(ba->IsMapped() && ba->IsMateMapped() && ba->IsProperPair())
		{

//			cerr << endl << ba->QueryBases << '\t' << ba->AlignedBases << '\t' << refv[ba->RefID].RefName << '\t' << ba->Position << '\t' << ba->GetEndPosition(false, true) << '\t' << ba->GetEndPosition(false, true) - ba->Position << '\t' << ba->AlignedBases.size() << endl;

			for(unsigned int i = 0; i < ba->AlignedBases.size(); i++)
			{
				if(toupper(ba->AlignedBases.at(i)) == 'N')
					continue;

				map<unsigned int, position>::iterator pi = GMAP[ba->RefID].find(ba->Position + i);

				if(pi == GMAP[ba->RefID].end())
				{
					mtx.lock();
					pi = GMAP[ba->RefID].insert(make_pair(ba->Position + i, position(toupper(ba->AlignedBases.at(i))))).first;
					mtx.unlock();
				}
				else
					pi->second.assign_nuc(toupper(ba->AlignedBases.at(i)));
			}
		}
	}

	return;
}*/

bool rnaseq::alignment_quality_check(BamAlignment *ba)
{
//	if(ba->QueryBases.size() != ba->AlignedBases.size())
//		return false;

	//if(ba->QueryBases.find("N") != string::npos || ba->AlignedBases.find("N") != string::npos)
	//	return false;

//	if(ba->QueryBases.find("n") != string::npos || ba->AlignedBases.find("n") != string::npos)
  //              return

	double mm_count = 0, match_count = 0;

//	cerr << "\nQuery:" << endl << ba->QueryBases << endl << "\nAligned:" << endl << ba->AlignedBases << endl;

//	cerr << endl << "CIGAR" << endl;

	unsigned int a_i = 0, q_i = 0;

	for(unsigned int i = 0; i < ba->CigarData.size(); i++)
	{
//		cerr << ba->CigarData.at(i).Type << '-' << ba->CigarData.at(i).Length << '\t';

//		if(ba->CigarData.at(i).Type == 'M')
		switch (ba->CigarData.at(i).Type)
		{
			case 'M':
			{
				for(unsigned int j = 0; j < ba->CigarData.at(i).Length; j++)
				{
					if(toupper(ba->QueryBases.at(q_i)) != toupper(ba->AlignedBases.at(a_i)))
						mm_count++;
					else
						match_count++;

					a_i++;
					q_i++;
				}

				break;
			}	
//		else if(ba->CigarData.at(i).Type == 'S')
			case 'S':
			{
				q_i += ba->CigarData.at(i).Length;
				break;
			}
//		else if(ba->CigarData.at(i).Type == 'N')
			case 'N':
			{
                        	a_i += ba->CigarData.at(i).Length;
				break;
			}
//		else if(ba->CigarData.at(i).Type == 'D')
			case 'D':
			{
				a_i += ba->CigarData.at(i).Length;
				break;
			}
//		else 
			default:
			{
				return false;
				break;
			}
		}
	}

//	cerr << "mcount = " << mm_count << '\t' << "match = " << match_count << endl;


/*	for(unsigned int i = 0; i < ba->QueryBases.size(); i++)
	{
		if(toupper(ba->QueryBases.at(i)) == 'N' || toupper(ba->AlignedBases.at(i)) == 'N')
			continue;

		if(toupper(ba->QueryBases.at(i)) != toupper(ba->AlignedBases.at(i)))
			mm_count++;
	}
*/
//	if((double)(mm_count / (double)ba->QueryBases.size()) > HQ_THRESHOLD)
	if(mm_count / (mm_count + match_count) > HQ_THRESHOLD || (double)((mm_count + match_count) / (double)ba->QueryBases.size()) < MATCH_FRAC)
		return false;

	return true;
}

void rnaseq::process_alignment(BamAlignment *ba, gmap *G)
{
//	cerr << endl << ba->QueryBases << '\n' << ba->AlignedBases << endl;

	if(ba->RefID == -1)
		return;

	if(USE_ONLY_HI_QUALITY_ALIGNMENTS)
		if(!alignment_quality_check(ba))
			return;

	auto vi = var_map.find(refv.at(ba->RefID).RefName);  //USE PREFETCH VCF TO SPEED UP

	if(vi == var_map.end())
		return;

//	cerr << endl << "vi->second.size() = " << vi->second.size();

//	cerr << endl << "refv.size() = " << refv.size() << '\t' << "RefID = " << ba->RefD;

	if(ba->IsPaired())
	{
		if(!USE_ONLY_FIRST_MATE || ba->IsFirstMate())
		{
			if(ba->IsMapped() && ba->IsMateMapped() && ba->IsProperPair())
			{
				for(unsigned int i = 0; i < ba->AlignedBases.size(); i++)
				{
					if(toupper(ba->AlignedBases.at(i)) == 'N')
						continue;


		/*		map<unsigned int, position>::iterator pi = GMAP[ba->RefID].find(ba->Position + i);

				if(pi == GMAP[ba->RefID].end())
				{
					mtx.lock();
					pi = GMAP[ba->RefID].insert(make_pair(ba->Position + i, position(toupper(ba->AlignedBases.at(i))))).first;
					mtx.unlock();
				}
				else
					pi->second.assign_nuc(toupper(ba->AlignedBases.at(i)));
		*/
					unsigned int pos = ba->Position + i;

					if(pos >= vi->second.size())
						break;

					if(!vi->second.at(pos))
						continue;           //////

					G->assign(ba->RefID, pos, ba->AlignedBases.at(i));
				}
			}
		}
	}
	else
	{
		if(ba->IsMapped())
		{
			for(unsigned int i = 0; i < ba->AlignedBases.size(); i++)
			{
				if(toupper(ba->AlignedBases.at(i)) == 'N')
					continue;

				unsigned int pos = ba->Position + i;

				if(pos >= vi->second.size())
					break;

                              	if(!vi->second.at(pos))
                                	continue;                //////

				G->assign(ba->RefID, pos, ba->AlignedBases.at(i));
			}
		}
	}

	return;
}
/*
void rnaseq::gmap_init()
{
	GMAP.resize(refv.size(), map<unsigned int, position>());

	return;
}
*/
void rnaseq::build_gmap_from_bam_multi()
{
	chrono::milliseconds span (1);
	vector<vector<BamAlignment> > AV (THREAD_NUM, vector<BamAlignment> (BUF_ALIGNMENTS));
	vector<gmap> LOC_GMAP(THREAD_NUM, gmap(refv.size()));

/*
	BamReader bam;

	bam.Open(rna_seq_bam_file.c_str());

        if(!bam.IsOpen())
        {
                cerr << "\nCan't find file: " << rna_seq_bam_file << "..." << endl;
                exit(1);
        }

	refv = bam.GetReferenceData();
*/
	cerr << "\nReading " << rna_seq_bam_file << " and building genomic map using " << THREAD_NUM << " threads...";

//	gmap_init();

	bool aflag = true;

	while(aflag)
	{
		BamAlignment ba;
		int tr = busy_check();

		while(tr < 0)
		{
			for(unsigned int f = 0; f < FUTV.size(); f++)
			{
				future_status boh = FUTV.at(f).wait_for(span); 

//				if(FUTV.at(f).wait_for(span) == future_status::ready)	
				if(boh == future_status::ready)
				{
//					cerr << endl << "Thread N " << f+1 << " has finished and is now ready for a new job " << endl;
					busy_vector.at(f) = FUTV.at(f).get();	
					AV[f].clear();
					AV[f].resize(BUF_ALIGNMENTS);
					break;
				}
					
			}

			tr = busy_check();
		}

		for(unsigned int buf_count = 0; buf_count < BUF_ALIGNMENTS && aflag; )
		{
			aflag = bam.GetNextAlignmentCore(ba);	

                	if(ba.RefID >= 0 && aflag)
			{
				AV.at(tr).at(buf_count) = ba;
				buf_count++;
			}
		}
	
		busy_vector.at(tr) = true;
//		cerr << endl << "Starting thread N " << tr+1 << endl;
//		FUTV.at(tr) = async (launch::async, &rnaseq::build_gmap_from_bam_multi_sub, this, &AV[tr], &LOC_GMAP[tr]);
		FUTV.at(tr) = async (launch::async, &rnaseq::build_gmap_from_bam_multi_sub, this, &AV[tr], &GMAP);
	}

//	bam.Close();
	
	for(unsigned int flast = 0; flast < FUTV.size(); flast++)
	{
		if(busy_vector.at(flast))
			busy_vector.at(flast) = FUTV.at(flast).get();

//		cerr << endl << "Thread N " << flast+1 << " has finished " << endl;
	}

	AV.clear();

/*	cerr << endl << "MERGING..." << endl;
	
	for(unsigned int i = 0; i < LOC_GMAP.size(); i++)
		GMAP.merge(LOC_GMAP.at(i));
*/
	return;
}


bool rnaseq::build_gmap_from_bam_multi_sub(vector<BamAlignment> *V, gmap *LOC_GMAP)
{
	for(unsigned int i = 0; i < V->size(); i++)
        {
                V->at(i).BuildCharData();
                process_alignment(&V->at(i), LOC_GMAP);
        }

	return false;
}

int rnaseq::busy_check()
{
	for(int i = 0; i < (int)busy_vector.size(); i++)
		if(!busy_vector.at(i))
			return i;

	return -1;
}
/*
void rnaseq::build_gmap_from_bam_multi3()
{
	vector<vector<BamAlignment> > AV;
	BamReader bam;

	bam.Open(rna_seq_bam_file.c_str());

        if(!bam.IsOpen())
        {
                cerr << "\nCan't find file: " << rna_seq_bam_file << "..." << endl;
                exit(1);
        }

        cerr << "\nReading " << rna_seq_bam_file << "..." << endl;

	refv = bam.GetReferenceData();

	gmap_init();

	AV.resize(refv.size(), vector<BamAlignment>());

	BamAlignment ba;

	while(bam.GetNextAlignmentCore(ba))
		if(ba.RefID >= 0)
			AV.at(ba.RefID).push_back(ba);

	cerr << "\nBuilding genomic map using " << THREAD_NUM << " threads..." << endl;

	vector<thread> VT;

	for(unsigned int i = 0; i < AV.size(); )
	{
		VT.clear();
		
		for(unsigned int j = 0; j < THREAD_NUM && i < AV.size(); j++)
		{
			VT.push_back(thread (&rnaseq::build_gmap_from_bam_multi3_sub, this, &AV[i]));
			i++;
		}

		for(unsigned int tn = 0; tn < VT.size(); tn++)
			VT[tn].join();
	}

	bam.Close();

	return;
}

void rnaseq::build_gmap_from_bam_multi3_sub(vector<BamAlignment> *V)
{
	for(unsigned int i = 0; i < V->size(); i++)
	{
		V->at(i).BuildCharData();
		process_alignment(&V->at(i));	
	}

	V->clear();
	V->resize(0);
	
	return;
}

void rnaseq::build_gmap_from_bam_multi2()
{
	BamReader bam;

	bam.Open(rna_seq_bam_file.c_str());

	if(!bam.IsOpen())
        {
                cerr << "\nCan't find file: " << rna_seq_bam_file << "..." << endl;
                exit(1);
        }

	cerr << "\nReading " << rna_seq_bam_file << " and building genomic map...\n" << endl;

	refv = bam.GetReferenceData();

	bam.Close();
	
	vector<thread> TV;

	for(unsigned int i = 0; i < refv.size(); i++)
	{
		TV.clear();

		for(unsigned int t = 0; t < THREAD_NUM; t++)
		{
			if(i < refv.size())
			{
				TV.push_back(thread (&rnaseq::build_gmap_from_bam_multi2_sub, this, i));
				i++;
			}

			for(unsigned int j = 0; j < TV.size(); j++)
                                TV.at(j).join();
		}
	}

	return;
}

void rnaseq::build_gmap_from_bam_multi2_sub(unsigned int refid)
{
	BamReader bam;

//	mtx.lock();
        bam.Open(rna_seq_bam_file.c_str());	
//	mtx.unlock();

	if(!bam.IsOpen())
        {       
//		mtx.lock();
                cerr << "\nCan't find file: " << rna_seq_bam_file << "..." << endl;
//		mtx.unlock();
                exit(1);
        }

	BamAlignment ba;
	bool aflag = true;
        
 //     while(bam.GetNextAlignment(ba))
	while(aflag)
	{
//		mtx.lock();
		aflag = bam.GetNextAlignment(ba);		

//		if((unsigned int)ba.RefID == refid && aflag)
  //              	process_alignment(&ba);
//		mtx.unlock();
	}
        
//	mtx.lock();
        bam.Close();
//	mtx.unlock();

	return;
}
*/
/*
void rnaseq::build_gmap_from_bam_multi()
{
	BamReader bam;

	bam.Open(rna_seq_bam_file.c_str());

	if(!bam.IsOpen())
        {
        	cerr << "\nCan't find file: " << rna_seq_bam_file << "..." << endl;
                exit(1);
        }

	cerr << "\nReading " << rna_seq_bam_file << " and building genomic map...\n" << endl;

	refv = bam.GetReferenceData();
	vector<BamAlignment> bav (THREAD_NUM, BamAlignment());
	vector<thread> TV;	
	bool aflag = true;

	while(aflag)
	{
		for(unsigned int i = 0; i < THREAD_NUM; i++)
		{
			TV.clear();
				
			aflag = bam.GetNextAlignment(bav[i]);

			if(aflag)
				TV.push_back(std::thread (&rnaseq::process_alignment, this, &bav.at(i)));

			for(unsigned int t = 0; t < TV.size(); t++)
                        	TV.at(t).join();
		}
	}
	

	bam.Close();

	return;
}
*/

// END OF RNASEQ CLASS

// START OF GMAP CLASS

gmap::gmap(unsigned int size)
{
	init(size);

	return;
}

const position* gmap::get_position_ptr(unsigned int ref_sequence, unsigned int pos) const
{
	map<unsigned int, position>::const_iterator mi = GENOME.at(ref_sequence).find(pos);

	if(mi != GENOME.at(ref_sequence).end())
		return &mi->second;
	else
		return NULL;
}

position* gmap::get_position_ptr(unsigned int ref_sequence, unsigned int pos) 
{
        map<unsigned int, position>::iterator mi = GENOME.at(ref_sequence).find(pos);

        if(mi != GENOME.at(ref_sequence).end())
                return &mi->second;
        else
                return NULL;
}

gmap::gmap()
{
	return;
}

void gmap::init(unsigned int size)
{
	GENOME.resize(size, map<unsigned int, position >());

	return;
}

void gmap::assign(unsigned int ref, unsigned int pos, char base)
{
	map<unsigned int, position>::iterator mi = GENOME.at(ref).find(pos);

	if(mi == GENOME.at(ref).end())
	{
		pair<unsigned int, position> PAIR = make_pair(pos, position(base));
		mtx.lock();
		mi = GENOME.at(ref).find(pos);     //DOUBLE CHECK, AVOID DATA RACE <--- HOLY GRAIL
		if(mi == GENOME.at(ref).end())	
		{
//			pair<unsigned int, position> PAIR = make_pair(pos, position(base));
			GENOME.at(ref).emplace(PAIR);
			mtx.unlock();
		}
		else
		{
			mtx.unlock();
			mi->second.assign_nuc_atomic(base);
		}
	}
	else
		mi->second.assign_nuc_atomic(base);
//		mi->second.assign_nuc(base);

/*
	pair<unsigned int position> PP = make_pair(pos, position(base));

	mtx.lock();

//	pair<map<unsigned int,position>::iterator, bool> PAIR = GENOME.at(ref).insert(make_pair(pos, position(base)));
	pair<map<unsigned int,position>::iterator, bool> PAIR = GENOME.at(ref).insert(PP);

	mtx.unlock();

	if(!PAIR.second)
		PAIR.first->second.assign_nuc_atomic(base);
*/

	return;
}

void gmap::report(const RefVector *refv)
{
	cerr << STDERR_INFO_SEPARATOR << "CHR\tNUMBER_OF_POSITIONS\n\n";

	for(unsigned int i = 0; i < GENOME.size(); i++)
		cout << refv->at(i).RefName << '\t' << GENOME[i].size() << endl;

	 cerr << STDERR_INFO_SEPARATOR;

	return;
}

void gmap::report_detailed(const RefVector *refv)
{
	ofstream out;

	if(RNAMAP_REP_FILE.empty())
		return;
	else
		out.open(RNAMAP_REP_FILE.c_str(), ios::out);
	
	if(!out)
	{
		cerr << "\nCan't open " << RNAMAP_REP_FILE << " for writing." << endl;
		return;
	}

	cerr << endl << "Writing " << RNAMAP_REP_FILE << " ..." << endl;

	for(unsigned int i = 0; i < GENOME.size(); i++)
	{
		map<unsigned int, position>::iterator mi = GENOME.at(i).begin();

		while(mi != GENOME.at(i).end())
		{
			out << refv->at(i).RefName << '\t' << mi->first << mi->second.report_s_atomic() << endl;
			mi++;
		}
	}	

	out.close();

	return;
}
/*

void gmap::merge(gmap &M)
{
	if(GENOME.size() != M.GENOME.size())
	{
		cerr << "\nFatal error: trying to merge GMAPS with different number of chromosomes... " << endl;
		exit(EXIT_FAILURE);
	}

	pair<map<unsigned int, position>::iterator, bool> PAIR;

	for(unsigned int i = 0; i < GENOME.size(); i++)	
	{
		map<unsigned int, position>::iterator mi = M.GENOME.at(i).begin();	

		while(mi != M.GENOME.at(i).end())
		{
			PAIR = GENOME.at(i).insert(make_pair(mi->first, position(mi->second.n_arr())));

			if(!PAIR.second)
				PAIR.first->second.add(mi->second.n_arr());	


			mi++;
		}
	}

	return;
}*/
// END OF GMAP CLASS

// BEGIN OF VCF_FILE CLASS

vcf_file::vcf_file(const string *FILE_VCF, const RefVector *refv)
{
	type_count.resize(variant::type_num, 0);

	in.open(FILE_VCF->c_str());
	string line;

	if(!in)
	{
		cerr << "\nCan't find VCF file: " << *FILE_VCF << endl;
		exit(EXIT_FAILURE);
	}

	cerr << "\nReading " << *FILE_VCF << "...";

	VV_MAP.resize(THREAD_NUM, vector<variant>());
	unsigned int vp = 0;
	unsigned int var_counter = 0; // TO BE DEALED WITH FOR MULTITHREAD
	
	while(getline(in,line))
	{	
		if(line.empty())
			continue;
		
		if(line[0] == '#')
			continue;

		VV_MAP.at(vp).emplace_back(variant(&line, refv, var_counter));	

//		cerr << endl << "TYPE = " << VV_MAP.at(vp).back().get_type() << endl;

		type_count.at(VV_MAP.at(vp).back().get_type())++;

		vp++;
		var_counter++;

		if(vp == THREAD_NUM)
			vp = 0;
	}	

	in.close();

	if(VERBOSE)
		type_report();

	return;
}

/*const vector<vector<variant> > * vcf_file::get_VV_MAP_ptr()
{
	return &VV_MAP;
}*/

unsigned int vcf_file::VV_MAP_size()
{
	return VV_MAP.size();
}

unsigned int vcf_file::V_MAP_size(unsigned int VV)
{
	return VV_MAP.at(VV).size();
}

const variant * vcf_file::get_variant_at(unsigned int VV , unsigned int V) const
{
	return &VV_MAP.at(VV).at(V);
}

void vcf_file::type_report()
{
	cerr << STDERR_INFO_SEPARATOR << "VARIANT_TYPE\tNUM" << endl; 

	static const string type_name[variant::type_num] = {"SNP", "INS", "DEL", "OTHER", "FILTERED"};

	for(unsigned int i = 0; i < type_count.size(); i++)
	{
		cerr << endl;

		if(i < variant::type_num)
			cerr << type_name[i];
		else
			cerr << "???";

		cerr << '\t' << type_count[i];
	}

	cerr << STDERR_INFO_SEPARATOR << endl;

	return;

}

// END OF VCF_FILE CLASS

// BEGIN OF ANALYSIS CLASS

analysis::analysis(const string *File_RNA, const string *File_VCF, const string *File_GTF, const string *File_Gene_Table, map<string, vector<bool> >& var_map):RNA(File_RNA, var_map),VCF(File_VCF, RNA.refv_ptr()), GTF(File_GTF, File_Gene_Table, RNA.refv_ptr())
{
	flag_var_positions();
	assign_transcript_to_pos();
	add_snp_counts_to_genes();

	string tr_report, snp_report, marker_var_tr_report, g_report, tr_final, tr_var_report, gene_var_report, gene_level;
	tr_report = snp_report = marker_var_tr_report = g_report = tr_final = tr_var_report = gene_var_report = gene_level = RNA.get_file_name();
	snp_report += "_snp_report.txt";
	tr_report += "_tr_report.txt";
	marker_var_tr_report += "_marker_report.txt";
	g_report += "_gene_report.txt";
	tr_final += "_tr_as_report.txt";
	tr_var_report += "_tr_var_report.txt";
	gene_var_report += "_gene_var_report.txt";
	gene_level += "_gene_level.txt";

//	simple_snp_report(snp_report);	
	gene_level_analysis(gene_level);
//	simple_transcript_report(tr_report);
//	gene_report(g_report);
//	tr_final_prob_report(tr_final);
//	tr_and_gene_var_report(gene_var_report, tr_var_report);
//	marker_var_transcript_report(marker_var_tr_report);
//	cout << "List of AS transcripts:\n" << asonly_transcripts() << endl;

	return;
}

void analysis::gene_level_analysis(const string out_file) const
{
	ofstream out(out_file.c_str());
	multimap <double, string> mm_out_buf;

        if(!out)
        {
                cerr << "\nCan not open " << out_file << " for writing... ";
                exit(EXIT_FAILURE);
        }

        cerr << endl << "Writing report file: " << out_file << endl;
	out << "GENE_NAME\tHZ_POS_N\tHZ_POS\tH_HAPLO\tL_HAPLO\tH_COV\tL_COV\tPV\tFDR\n";

	for(auto gpi = GENE_POS.cbegin(); gpi != GENE_POS.cend(); ++gpi)
	{
		ostringstream out_buf;

		out_buf << gpi->first->id() << '\t' << gpi->second.size() << '\t';

		out_buf << gpi->first->get_chr() << ':';

		set<position*>::const_iterator pptr = gpi->second.cbegin();

		map<unsigned int, position*> pos_map; //required to sort

		while(pptr != gpi->second.cend())
		{
			pos_map.insert(make_pair((*pptr)->get_var()->get_pos(), *pptr));
			pptr++;
		}

		map<unsigned int, position*>::const_iterator pptr_m =  pos_map.cbegin();

//		while(pptr != gpi->second.cend())
		while(pptr_m != pos_map.cend())
		{
//			out_buf << (*pptr)->get_var()->get_pos() << '|';
			out_buf << pptr_m->first << '|';

//			pptr++;
			pptr_m++;
		}

		out_buf << '\t';

//
//	 	pptr = gpi->second.cbegin();        

		pptr_m =  pos_map.cbegin();

//		while(pptr != gpi->second.cend())
		while(pptr_m != pos_map.cend())
                {
//			if((*pptr)->get_var()->is_phased())
//				out_buf << (*pptr)->get_var()->get_phase() << '|';
			if(pptr_m->second->get_var()->is_phased())
 	                       out_buf << pptr_m->second->get_var()->get_phase() << '|';
			else
			{
	//			if((*pptr)->get_ref_count() >= (*pptr)->get_alt_count())
				if(pptr_m->second->get_ref_count() >= pptr_m->second->get_alt_count())
					out_buf << 0 << '/';
				else
					out_buf << 1 << '/';	
			}

                  //      pptr++;
                  	  pptr_m++;
                }

		out_buf << '\t'; 

//		pptr = gpi->second.cbegin();
		pptr_m =  pos_map.cbegin();

         //       while(pptr != gpi->second.cend())
         	while(pptr_m != pos_map.cend())
                {
                //        if((*pptr)->get_var()->is_phased())
                	if(pptr_m->second->get_var()->is_phased())
                                out_buf << !pptr_m->second->get_var()->get_phase() << '|';
                        else
                        {
                        //        if((*pptr)->get_ref_count() >= (*pptr)->get_alt_count())
                        	if(pptr_m->second->get_ref_count() >= pptr_m->second->get_alt_count())
                                        out_buf << 1 << '/';
                                else
                                        out_buf << 0 << '/';
                        }

                  //      pptr++;
                  	  pptr_m++;
                }

                out_buf << '\t';

//
		if(gpi->first->get_big_c() > gpi->first->get_small_c())   // FOR PHASED ANALYSIS
			out_buf << gpi->first->get_big_c() << '\t' << gpi->first->get_small_c();
		else
			out_buf << gpi->first->get_small_c() << '\t' << gpi->first->get_big_c();

		gpi->first->gtest_2();

		out_buf /*<< '\t' << gpi->first->get_min_pv() << '\t' << gpi->first->get_max_pv()*/  << '\t' << gpi->first->get_tot_pv();

		mm_out_buf.insert(make_pair(gpi->first->get_tot_pv(), out_buf.str()));
		
	}

	auto mm_ci = mm_out_buf.cbegin();
	double rank = 1;
	vector<double> prov_fdr;

	while(mm_ci !=  mm_out_buf.cend())
	{
		double mult = (double)mm_out_buf.size() / rank;	
		prov_fdr.push_back(mm_ci->first * mult); 

//		out << mm_ci->second << endl; 

		rank++;
		mm_ci++;
	}

	mm_ci = mm_out_buf.cbegin();
	unsigned int vpos = 0;

	while(mm_ci !=  mm_out_buf.cend())
	{
		double FDR = prov_fdr.at(vpos);

		if(vpos < prov_fdr.size() - 1)
		{
			if(FDR > prov_fdr.at(vpos + 1))
				FDR = prov_fdr.at(vpos + 1);
		}

		if(FDR > 1)
			FDR = 1;

		out << mm_ci->second << '\t' << FDR << endl;;

		vpos++;
		mm_ci++;
	}

	out.close();

	return;
}

void analysis::add_snp_counts_to_genes()
{
	for(auto gpi = GENE_POS.begin(); gpi != GENE_POS.end(); ++gpi)
	{
		if(PHASE_ANALYSIS)
		{
			for(auto spi = gpi->second.begin(); spi != gpi->second.end(); ++spi)
				gpi->first->add_snp_c((*spi)->get_ref_count(), (*spi)->get_alt_count(), (*spi)->get_phase());
		}
		else
		{
			for(auto spi = gpi->second.begin(); spi != gpi->second.end(); ++spi)
                                gpi->first->add_snp_c((*spi)->get_ref_count(), (*spi)->get_alt_count());
		}
	}

	return;
}

void analysis::tr_and_gene_var_report(const string tr_var_report_file, const string gene_var_report_file) const
{
	ofstream out1(tr_var_report_file.c_str()), out2(gene_var_report_file.c_str()); 	

	cerr << endl << "Writing report files: " << tr_var_report_file << " and " << gene_var_report_file << endl;

	if(!out1)
        {
                cerr << "\nCan not open " << tr_var_report_file << " for writing... ";
                exit(EXIT_FAILURE);
        }

	if(!out2)
        {
                cerr << "\nCan not open " << gene_var_report_file << " for writing... ";
                exit(EXIT_FAILURE);
       	} 

	const map<string, map<string, gene> > *CHR = GTF.get_genes_ptr();

	map<string, map<string, gene> >::const_iterator ci = CHR->cbegin();

	while(ci != CHR->cend())
        {
		map<string, gene>::const_iterator gi = ci->second.cbegin();

		while(gi != ci->second.cend())
		{
			if(!gi->second.has_transcript_with_var())
                        {
                                gi++;
                                continue;
                        }

			out1 << gi->first << '\t' << gi->second.report_var_counts() << endl;

		        const map<string, transcript> *TRANSCRIPTS = gi->second.get_transcripts_ptr();

                        map<string, transcript>::const_iterator ti = TRANSCRIPTS->cbegin();

                        while(ti != TRANSCRIPTS->cend())
			{
				if(ti->second.get_expr() > transcript::NO_EXPRESSION)
					out2 << gi->first << '\t' << ti->first << '\t' << ti->second.snp_count_report() << endl;

				ti++;
			}

			gi++;
		}	

		ci++;
	}

	out1.close();
	out2.close();

	return;
}

void analysis::tr_final_prob_report(const string out_file) const
{
	ofstream out(out_file.c_str());

        if(!out)
        {
                cerr << "\nCan not open " << out_file << " for writing... ";
                exit(EXIT_FAILURE);
        }

	cerr << endl << "Writing report file: " << out_file << endl;

        const map<string, map<string, gene> > *CHR = GTF.get_genes_ptr();

        map<string, map<string, gene> >::const_iterator ci = CHR->cbegin();

        while(ci != CHR->cend())
        {
                map<string, gene>::const_iterator gi = ci->second.cbegin();

                while(gi != ci->second.cend())
                {
			if(!gi->second.has_transcript_with_var())
			{
				gi++;
				continue;	
			}

			const map<string, transcript> *TRANSCRIPTS = gi->second.get_transcripts_ptr();

                	map<string, transcript>::const_iterator ti = TRANSCRIPTS->cbegin();

                        while(ti != TRANSCRIPTS->cend())
			{
				if(ti->second.get_count_var() == 0)
				{
					ti++;
					continue;
				}

				out << gi->first << '\t' << ti->first << '\t' << ti->second.get_count_var() << '\t' << ti->second.get_prod_prob(transcript::ALLELE_SPECIFIC) << '\t' <<  ti->second.get_prod_prob(transcript::NOT_ALLELE_SPECIFIC) << "\t---\t" << ti->second.all_prob_report(transcript::ALLELE_SPECIFIC) << '\t' << ti->second.all_prob_report(transcript::NOT_ALLELE_SPECIFIC)<< endl;

				ti++;
			}

			gi++;
		}

		ci++;
	}

	out.close();

	return;
}

void analysis::flag_var_positions()
{
	for(unsigned int i = 0; i < VCF.VV_MAP_size(); i++)
        {
                for(unsigned int j = 0; j < VCF.V_MAP_size(i); j++)
                {
                        const variant *V = VCF.get_variant_at(i,j);

			if(V->get_type() == variant::SNP && V->get_chr_ref() != variant::INVALID_CHR)
			{	
				position *POS = RNA.get_position_ptr((unsigned int)V->get_chr_ref(), V->get_pos());

				if(POS != NULL)
					POS->set_var(V);
			}
		}
	}

	return;
}

const string analysis::asonly_transcripts() const
{
	ostringstream out;

	const map<string, map<string, gene> > *CHR = GTF.get_genes_ptr();

	map<string, map<string, gene> >::const_iterator ci = CHR->cbegin();

	while(ci != CHR->cend())
	{
		map<string, gene>::const_iterator gi = ci->second.cbegin();

		while(gi != ci->second.cend())
		{
			const map<string, transcript> *TRANSCRIPTS = gi->second.get_transcripts_ptr();

			map<string, transcript>::const_iterator ti = TRANSCRIPTS->cbegin();

			while(ti != TRANSCRIPTS->cend())
			{
				if(ti->second.get_count_var() > 0 && !ti->second.get_aspec_flag())
					out << ti->first << '\t' << ti->second.get_count_var()  << '\t' << ti->second.get_gene_id() << endl;	

				ti++;
			}

			gi++;
		}

		ci++;
	}

	return out.str();
}

void analysis::gene_report(const string out_file)
{
	ofstream out(out_file.c_str());
	
	if(!out)
        {
                cerr << "\nCan not open " << out_file << " for writing... ";
                exit(EXIT_FAILURE);
        }

	cerr << endl << "Writing report file: " << out_file << endl;

	const map<string, map<string, gene> > *CHR = GTF.get_genes_ptr();	
	
	map<string, map<string, gene> >::const_iterator ci = CHR->cbegin(); 

	while(ci != CHR->cend())
	{
		map<string, gene>::const_iterator gi = ci->second.cbegin();

		while(gi != ci->second.cend())
		{
			bool ot = false, mt = false;
			set<position *> USED_POS;

			const map<string, transcript> *TRANSCRIPTS = gi->second.get_transcripts_ptr();

			map<string, transcript>::const_iterator ti = TRANSCRIPTS->cbegin();

//			out << gi->first << endl;

			bool gflag = false;
	
			while(ti != TRANSCRIPTS->cend())
			{
				for(unsigned int i = 0; i < ti->second.get_pos_ptrv_size(); i++)
				{
					position *POS= RNA.get_position_ptr(ti->second.get_pos_ptrv_chr_ref_at(i), ti->second.get_pos_ptrv_pos_at(i));

					if(POS == NULL)
					{
						cerr << "\nWarning, something weird occurring (POS)..." << endl;
						continue;
					}

					const variant *V = POS->get_var();

					if(V == NULL)
					{
						cerr << "\nWarning, something weird occurring (VAR)..." << endl;
                                                continue;
					}

					if(USED_POS.find(POS) != USED_POS.end())
						continue;

					USED_POS.insert(POS);

					if(!gflag)
					{
	//					out << gi->first;
						out << gi->second.report() << endl;
						gflag = true;
					}

		//			out << ti->first << '\t' << i+1 << '/' << ti->second.get_pos_ptrv_size() << '\t' << position_full_report(V, POS) << endl;	 
					out << position_small_report(V, POS);

					if(V->get_zygosity() == variant::ETERO_Z)
						gi->second.add_etero_var_to_count();
					else if(V->get_zygosity() == variant::ALT_OMO_Z)
						gi->second.add_omo_var_to_count();
			
					if(V->get_zygosity() == variant::ETERO_Z && !POS->multiple_gene())	
					{
					//	out << POS->report_b_comb()  << endl;
						vector<double> PV, V_V1, G_V, CHI_RAW;
						double Att[2], gtest;
						unsigned int bvg;
	
						POS->compute_probability_binom(PV, V_V1, Att, V->get_ref_int(), V->get_alt_int());
						gtest = POS->compute_probability_gtest(V_V1, G_V, CHI_RAW, bvg, V->get_ref_int(), V->get_alt_int());
						POS->pass_prob_values_to_tr(CHI_RAW);

						if(POS->trv_size() == 1)
							ot = true;
						else if(POS->trv_size() > 1)
							mt = true;

			//			for(unsigned int m = 0; m < PV.size(); m++)
			//				out << '\t' << POS->report_b_comb_at(m) << ':' << PV.at(m) << '\t';
					
						double max = -1;
						unsigned int pmax = 0;

						for(unsigned int m = 0; m < PV.size(); m++)
						{
							if(PV.at(m) > max)
							{
								max = PV.at(m);
								pmax = m;
							}
						}

						out << '\t' << POS->report_b_comb_at(pmax) << '\t' << PV.at(pmax) << '\t' << "---" << '\t' << Att[0] << '\t' << Att[1]; 

						out << "\t---\t" << POS->report_b_comb_at(bvg) << '\t' << gtest;

						if(pmax != bvg)
							out << '\t' << "***";

						POS->add_var_to_tr(pmax);
					}
					else if(POS->multiple_gene())
						out << "\tMULTIPLE_GENE";

					out << endl;
							
				//	     << V->get_internal_id() << '\t' << ci->first << '\t' << V->get_pos() << '\t' << V->get_ref() << '\t' << V->get_alt() << '\t' << *V->zygosity_str() 
                                  //           << POS->report_s_atomic() << POS->report_s_transcripts() << endl;
				}

				ti++;
			}

			if(ot && mt)
				out << "LOOK HERE" << endl;

			gi++;
		}

		ci++;
	}

	out.close();

	
	return;
}


void analysis::simple_transcript_report(const string out_file)
{
	ofstream out(out_file.c_str());
	
	if(!out)
        {
                cerr << "\nCan not open " << out_file << " for writing... ";
                exit(EXIT_FAILURE);
        }

	cerr << endl << "Writing report file: " << out_file << endl;

	const map<string, map<string, gene> > *CHR = GTF.get_genes_ptr();	
	
	map<string, map<string, gene> >::const_iterator ci = CHR->cbegin(); 

	while(ci != CHR->cend())
	{
		map<string, gene>::const_iterator gi = ci->second.cbegin();

		while(gi != ci->second.cend())
		{
			const map<string, transcript> *TRANSCRIPTS = gi->second.get_transcripts_ptr();

			map<string, transcript>::const_iterator ti = TRANSCRIPTS->cbegin();

			while(ti != TRANSCRIPTS->cend())
			{
				for(unsigned int i = 0; i < ti->second.get_pos_ptrv_size(); i++)
				{
					position *POS= RNA.get_position_ptr(ti->second.get_pos_ptrv_chr_ref_at(i), ti->second.get_pos_ptrv_pos_at(i));

					if(POS == NULL)
					{
						cerr << "\nWarning, something weird occurring (POS)..." << endl;
						continue;
					}

					const variant *V = POS->get_var();

					if(V == NULL)
					{
						cerr << "\nWarning, something weird occurring (VAR)..." << endl;
                                                continue;
					}


					out << gi->first << '\t' << ti->first << '\t' << ti->second.get_expr() /*<< '\t' << ti->second.get_expr_weight()*/ << '\t' << i+1 << '/' << ti->second.get_pos_ptrv_size() << '\t' << position_full_report(V, POS) << endl;	 
			
					if(V->get_zygosity() == variant::ETERO_Z && !POS->multiple_gene())	
					{
					//	out << POS->report_b_comb()  << endl;
						vector<double> PV, V_V1, G_V, CHI_RAW;
						double Att[2];
						unsigned int foo;
	
						POS->compute_probability_binom(PV, V_V1, Att, V->get_ref_int(), V->get_alt_int());
						POS->compute_probability_gtest(V_V1, G_V, CHI_RAW,  foo, V->get_ref_int(), V->get_alt_int());
							
						for(unsigned int m = 0; m < PV.size(); m++)
						{
							out << POS->report_b_comb_at(m) << ':' << PV.at(m) << '/' << G_V.at(m); 
					/*		
							if(!POS->perf_omo())
							{
								if(!G_V.empty())
									out << gsl_cdf_chisq_Q(2 * G_V.at(m), 1);
							}
							else
							{
								if(!G_V.empty())
									out << G_V.at(m);
							}
*/
							out << '\t';
						}
						
						out << endl;

						

						
							
					}
					else if(POS->multiple_gene())
						out << "MULTIPLE_GENE" << endl;
				}

				ti++;
			}

			gi++;
		}

		ci++;
	}

	out.close();

	
	return;
}
	
void analysis::marker_var_transcript_report(const string out_file)
{
	ofstream out(out_file.c_str());
	set<string> TRANSCRIPTS_WITH_MARKER;
	
	if(!out)
        {
                cerr << "\nCan not open " << out_file << " for writing... ";
                exit(EXIT_FAILURE);
        }

	cerr << endl << "Writing report file: " << out_file << endl;

	const map<string, map<string, gene> > *CHR = GTF.get_genes_ptr();	
	
	map<string, map<string, gene> >::const_iterator ci = CHR->cbegin(); 

	while(ci != CHR->cend())
	{
		map<string, gene>::const_iterator gi = ci->second.cbegin();

		while(gi != ci->second.cend())
		{
			const map<string, transcript> *TRANSCRIPTS = gi->second.get_transcripts_ptr();

			map<string, transcript>::const_iterator ti = TRANSCRIPTS->cbegin();

			while(ti != TRANSCRIPTS->cend())
			{
				for(unsigned int i = 0; i < ti->second.get_pos_ptrv_size(); i++)
				{
					position *POS= RNA.get_position_ptr(ti->second.get_pos_ptrv_chr_ref_at(i), ti->second.get_pos_ptrv_pos_at(i));

					if(POS == NULL)
					{
						cerr << "\nWarning, something weird occurring (POS)..." << endl;
						continue;
					}

					if(!POS->tr_marker())
						continue;

					const variant *V = POS->get_var();

					if(V == NULL)
					{
						cerr << "\nWarning, something weird occurring (VAR)..." << endl;
                                                continue;
					}


					out << gi->first << '\t' << ti->first << '\t' << i+1 << '/' << ti->second.get_pos_ptrv_size() << '\t' << position_full_report(V, POS) /*<< '\t' << ti->second.get_expr()*/ << endl;
					  //   << V->get_internal_id() << '\t' << ci->first << '\t' << V->get_pos() << '\t' << V->get_ref() << '\t' << V->get_alt() << '\t' << *V->zygosity_str() 
                                          //  << POS->report_s_atomic() << '\t' << *POS->zygosity_str() << '\t' << *var_diagnostic(V->get_zygosity(), POS->get_zygosity()) << POS->report_s_transcripts() << endl;
				
					TRANSCRIPTS_WITH_MARKER.insert(ti->first);
				}

				ti++;
			}

			gi++;
		}

		ci++;
	}

	cerr << endl << "NUMBER OF TRANSCRIPTS WITH AT LEAST ONE MARKER: " << TRANSCRIPTS_WITH_MARKER.size() << endl;

	out.close();
}

const string analysis::position_full_report(const variant *V, position *POS) 
{
	ostringstream out;

	out << V->get_internal_id() << '\t' << V->get_chr() << '\t' << V->get_pos() << '\t' << V->get_ref() << '\t' << V->get_alt() << '\t' << *V->zygosity_str()
            << POS->report_s_atomic() << '\t' << *POS->zygosity_str() << '\t' << *POS->var_diagnostic_s(V->get_zygosity()) /**var_diagnostic(V->get_zygosity(), POS->get_zygosity())*/ << POS->report_s_transcripts() << '\t';
	
	double rec_p, mpv;

	mpv = POS->marker_p(V->get_ref_int(), V->get_alt_int(), &rec_p);	

	if(V->get_zygosity() == variant::ETERO_Z)
		out << mpv << '\t' << rec_p;
	else
		out << "-\t-";
	
	return out.str();
}

const string analysis::position_small_report(const variant *V, position *POS) 
{
	ostringstream out;

	out << V->get_internal_id() << '\t' << V->get_chr() << '\t' << V->get_pos() << '\t' << *V->zygosity_str() << '\t' << V->get_ref() << '\t' << V->get_alt() << '\t' << POS->get_count_at(V->get_ref_int()) << '\t' << POS->get_count_at(V->get_alt_int()) << '\t' << POS->report_s_transcripts();

	return out.str();
}
	
void analysis::simple_snp_report(const string out_file)
{
	ofstream out(out_file.c_str());

	if(!out)
	{
		cerr << "\nCan not open " << out_file << " for writing... ";
		exit(EXIT_FAILURE);
	}

	cerr << "\nWriting report file: " << out_file << endl;

	for(unsigned int i = 0; i < VCF.VV_MAP_size(); i++)
	{
		for(unsigned int j = 0; j < VCF.V_MAP_size(i); j++)
		{
			const variant *V = VCF.get_variant_at(i,j);

			if(V->get_type() == variant::SNP && V->get_chr_ref() != variant::INVALID_CHR)
			{
				const position *POS = RNA.get_position_ptr((unsigned int)V->get_chr_ref(), V->get_pos());

				if(POS != NULL)
				{
					out << V->get_internal_id() << '\t' << V->get_chr() << '\t' << V->get_pos() << '\t' << V->get_ref() << '\t' << V->get_alt() << '\t' << *V->zygosity_str() 
					     << POS->report_s_atomic() << '\t' << POS->report_s_transcripts();

					if(PHASE_ANALYSIS && V->get_zygosity() == variant::ETERO_Z)
						out << "\tPHASE = " << V->get_phase();	

					out << endl;
				}
			}	
		}
	}

	out.close();

	return;
}

void analysis::assign_transcript_to_pos()
{
	cerr << endl << "Assigning transcripts to genomic positions..." << endl;

	map<string, map<string, gene> > *GENES = GTF.get_genes_ptr();

	map<string, map<string, gene> >::iterator mi = GENES->begin(); 

//	set<position*> int_pos_set;

	while(mi != GENES->end())
	{
		map<string, gene>::iterator gi = mi->second.begin();	

		while(gi != mi->second.end())	
		{
			if(gi->second.get_chr_ref() == gene::INVALID_CHR)
			{
				gi++;
				continue;
			}

			map<string, transcript> *TRANSCRIPTS = gi->second.get_transcripts_ptr();
			map<string, transcript>::iterator ti = TRANSCRIPTS->begin();		

			while(ti != TRANSCRIPTS->end())
			{
				const vector<exon> *EXONS = ti->second.get_exons_ptr();

				for(unsigned int i = 0; i < EXONS->size(); i++)
				{
					for(unsigned int p = EXONS->at(i).get_start(); p < EXONS->at(i).get_end(); p++)
					{
						position *POS = RNA.get_position_ptr(gi->second.get_chr_ref(), p);	

						if(POS != NULL)
						{
							if(POS->get_var() != NULL)
							{
								if(EXPR_FILE.empty() || ti->second.get_expr() > transcript::NO_EXPRESSION)
								{
									POS->add_transcript(&ti->second);
									ti->second.add_pos(gi->second.get_chr_ref(), p, POS->get_var()->get_zygosity());
			//						int_pos_set.insert(POS);
									
									if(POS->get_var()->get_zygosity() == variant::ETERO_Z)
									{
										if(POS->get_ref_count() || POS->get_alt_count())
										{
											auto gpi = GENE_POS.find(&gi->second);

											if(gpi == GENE_POS.end())
												gpi = GENE_POS.insert(make_pair(&gi->second, set<position*>())).first;

											gpi->second.insert(POS);
										}

				//						cerr << endl << gi->first << '\t' << POS->get_var()->get_internal_id() << '\t' << gpi->second.size();
									}
								}
							}
						}
					}
				}	

				ti++;
			}

			gi++;
		}

		mi++;
	}
/*
	for(auto pi = int_pos_set.begin(); pi != int_pos_set.end(); ++pi)
	{
		if((*pi)->get_var()->get_zygosity() == variant::ETERO_Z)
			(*pi)->compute_probability_tr();
	}
*/
	return;
}
/*
const string* analysis::var_diagnostic(unsigned int vcf_z, unsigned int rna_z) const
{
	static const vector<string> D = {"ETERO_ZYG_EXPR", "ALLELE_SPECIFIC_EXPR", "RNA_EDIT", "OMOZYG_EXPR", "DUBIOUS_EXPR"};

	if(vcf_z == variant::ETERO_Z && rna_z == position::ETERO_Z)
		return &D[0];
	else if(vcf_z == variant::ETERO_Z && rna_z == position::OMO_Z)
		return &D[1];
	else if(vcf_z == variant::ALT_OMO_Z && rna_z == position::ETERO_Z)
		return &D[2];
	else if(vcf_z == variant::ALT_OMO_Z && rna_z == position::OMO_Z)
                return &D[3];
	else
		return &D[4];
}
*/
// END OF ANALYSIS CLASS
//
// BEGIN OF VARIANT CLASS

variant::variant(const string *Line, const RefVector *ref_ptr, const unsigned int var_count):internal_id(var_count)
{
	refv = ref_ptr;

	PHASED = false;	
	VAR_CHR = -1;

	istringstream str(*Line), sgtfields, sgt;

//	cerr << *Line << endl;

	str >> chr >> pos >> id >> ref >> alt >> qual >> filter >> info >> gtfields;

	if(FIX_VCF)
		chr.insert(0, "chr");

	pos--;  //CONVERTING POSITION TO 0-BASED, LIKE SAM/BAM

	while(gtfields.find("GT") == string::npos)
		str >> gtfields;	

//	cerr << endl << "INFO = " << gtfields << endl;

	str >> gt;	

	sgtfields.str(gtfields);
	sgt.str(gt);

	string buf1, buf2;

	while(getline(sgtfields, buf1, ':'))
	{
		getline(sgt, buf2, ':');	

		GT[buf1] = buf2;
	}

	set_chr_ref();
	set_type();
	set_zygosity();

//	gt_report();
//	report();

//	cerr << endl << "TYPE = " << type << endl;

	return;
}

const short int variant::get_phase() const
{
	return VAR_CHR;
}

const bool variant::is_phased() const
{
	return PHASED;
}

const unsigned int variant::get_internal_id() const
{
	return internal_id;
}

const string variant::get_chr() const
{
	return chr;
}

const string variant::get_ref() const
{
	return ref;
}

const unsigned int variant::get_char_int(const string *s) const
{
	char cref = 'N';

	if(!s->empty())
		cref = s->at(0);

	switch (cref)
        {
                case 'A':
                return A;
                break;
                case 'C':
                return C;
                break;
                case 'G':
                return G;
                break;
                case 'T':
                return T;
                break;
                default:
		return N; 
		break;
        }
}

const unsigned int variant::get_ref_int() const
{
	unsigned int i = get_char_int(&ref);

	if(i == N)
	{
		cerr << endl << endl << "Invalid reference nucleotide found: " << ref  << " in variant at pos " << chr << '\t' << pos << endl;
                exit(EXIT_FAILURE);
	}

	return i;
}

const unsigned int variant::get_alt_int() const
{
	unsigned int i = get_char_int(&alt);

        if(i == N)
        {
                cerr << endl << endl << "Invalid alternative nucleotide found: " << ref  << " in variant at pos " << chr << '\t' << pos << endl;
                exit(EXIT_FAILURE);
        }

        return i;
}


const string variant::get_alt() const
{
	return alt;
}

void variant::set_chr_ref()
{
	static set<string> UNUSED_CHR;
	chr_ref = INVALID_CHR;

	for(unsigned int i = 0; i < refv->size(); i++)
	{
		if(refv->at(i).RefName == chr)
		{
			chr_ref = (int)i;
			return;
		}
	}

	if(UNUSED_CHR.find(chr) == UNUSED_CHR.end())	
	{
		cerr << "\nWarning, this reference sequence found in vcf file does not exist in the rnaseq file: " << chr;
		cerr << "\n-> All occurrences of " << chr << " will be ignored." << endl;
		UNUSED_CHR.insert(chr);
	}

	return;
}

const int variant::get_chr_ref() const
{
	return chr_ref;
}

void variant::set_zygosity()
{
	map<string, string>::iterator mi = GT.find("GT");

	if(mi == GT.end())
		zygosity = UNKNOWN_Z;

	for(unsigned int i = 0; i < mi->second.size(); i++)
	{
		if(mi->second.at(i) == '|')  
		{
			mi->second.at(i) = ' ';
			PHASED = true;
		}
		else if( mi->second.at(i) == '/')
		{
			mi->second.at(i) = ' ';
			PHASED = false;
		}

/*		else if(mi->second.at(i) == '/')
		{
			zygosity = UNPHASED;
			return;
		}*/
	}

	istringstream str(mi->second);	
	string v1, v2;

	str >> v1 >> v2;

	if(v1 == "." || v2 == ".")
		zygosity = UNKNOWN_Z;
	else if(v1 == "0" && v1 != v2)
	{
		zygosity = ETERO_Z;

		if(PHASED)
			VAR_CHR = 1;
	}
	else if(v1 == "1" && v1 != v2)
	{
		zygosity = ETERO_Z;

		if(PHASED)	
			VAR_CHR = 0;
	}
	else if(v1 == "0" && v1 == v2)
		zygosity = REF_OMO_Z;
	else if(v1 == "1" && v1 == v2)
		zygosity = ALT_OMO_Z;


	if(PHASE_ANALYSIS && !PHASED)
		zygosity = UNPHASED;

	return;
}

unsigned int variant::get_zygosity() const
{
	return zygosity;
}

void variant::set_type()
{
	if(!VCF_FILTER.empty() && filter == VCF_FILTER)
	{
		type = FILTERED;
		return;
	}

	for(unsigned int i = 0; i < ref.size(); i++)
	{
		bool flag = false;

		for(unsigned int j = 0; j < ALPHABET_SIZE; j++)
			if(ref.at(i) == ALPHABET[j])
			{
				flag = true;
				break;
			}

		if(!flag)
		{
			type = OTHER;
			return;
		}
	}

	if(ref.size() == 1 && alt.size() == 1)
		type = SNP;
	else if(ref.size() == 1 && alt.size() > 1)
		type = INS;
	else if(ref.size() > 1 && alt.size() == 1)
		type = DEL;
	else
		type = OTHER;

	return;
}

unsigned int variant::get_type() const
{
	return type;
}

unsigned int variant::get_pos() const
{
	return pos;
}

void variant::gt_report()
{
	map<string, string>::iterator mi = GT.begin();

	while(mi != GT.end())
	{
		cerr << endl << mi->first << '\t' << mi->second;

		if(mi->first == "GT")	
			cerr << "\tZYG = " << *zygosity_str();

		mi++;
	}

	cerr << endl;

	return;
}

void variant::report()
{
	cerr << endl << chr << '\t' << pos << '\t' << "CHR_REF = " << chr_ref;

	return;
}

const string* variant::zygosity_str() const
{
	static const string zig_s[zyg_types] = {"VCF_REF_OMOZYG", "VCF_ALT_OMOZYG", "VCF_ETEROZYG", "VCF_DOUBLE_ETEROZYG", "VCF_UNKNOWN_ZYG", "VCF_UNPHASED"};

	return &zig_s[zygosity];	
}

// END OF VARIANT CLASS

// GTF_FILE CLASS
gtf_file::gtf_file(const string *File_gtf, const string *File_Names, const RefVector *R):refv(R)
{
	const static set<string> valid_features = {"exon"};
	const static set<string> gene_keys = {"gene_id"};
	const static set<string> transcript_keys = {"transcript_id"};

	read_tr_to_gene_file(File_Names);

	ifstream in(File_gtf->c_str());	

	if(!in)
	{
		cerr << "\nCan't find GTF file: " << *File_gtf << endl;
		exit(EXIT_FAILURE);
	}
	else
		cerr << endl << "Reading GTF file: " << *File_gtf << endl;

	string line;

	while(getline(in,line))
	{
		if(line.empty())
			continue;
		if(line[0] == '#')
			continue;
		
		istringstream str(line);
		string chr, trash, feature_type, gene_id, tr_id;
		unsigned int start, end;
		char strand;

		str >> chr >> trash >> feature_type;

		if(valid_features.find(feature_type) == valid_features.end())
			continue;

		str >> start >> end >> trash >> strand;

		str.seekg(ios_base::beg);

		while(str >> gene_id)
                {
                        if(gene_keys.find(gene_id) == gene_keys.end())
                                continue;
                        else
                        {
                                str >> gene_id;
                                break;
                        }
                }

		if(!gene_id.empty())         
                	clean_ids(&gene_id);

		str.seekg(ios_base::beg);

		while(str >> tr_id)	
		{
			if(transcript_keys.find(tr_id) == transcript_keys.end())
				continue;
			else
			{
				str >> tr_id;
				break;
			}
		}

		if(!tr_id.empty())       
			clean_ids(&tr_id);
		else 				// DISCARD ENTRY IF NOT ABLE TO IDENTIFY TRANSCRIPT IDENTIFIER
		{
			cerr << "\nWarning: can not determine transcript id for this gtf entry: " << endl << line << endl << endl;
			continue;
		}

		map<string, string>::const_iterator name_i = TR_TO_GENE.find(tr_id);

		if(name_i != TR_TO_GENE.cend())   //IF TRANSCRIPT HAS CORRESPONDING GENE NAME IN TR_TO_GENE THEN SET gene_id TO THAT VALUE
			gene_id = name_i->second;

		auto mi = GENES.find(chr);

		if(mi == GENES.end())
			mi = GENES.insert(make_pair(chr, map<string, gene>())).first;

		auto gi = mi->second.find(gene_id);

		if(gi == mi->second.end())
			gi = mi->second.emplace(make_pair(gene_id, gene(&chr, &gene_id, strand, refv))).first;

		auto direct_tr_it = gi->second.add_transcript_or_exon(&tr_id, start, end);

		auto id_to_tr_i  = ID_TO_TR_MAP.find(tr_id);

		if(id_to_tr_i == ID_TO_TR_MAP.end())
			ID_TO_TR_MAP.insert(make_pair(tr_id, direct_tr_it));
	}

	in.close();

/*

	read_expression(EXPR_FILE, EXPRESSION);
	read_expression(GENE_EXPR_FILE, GENE_EXPRESSION);

	set_expression();
	set_expression_gene();

	set_expression_weight();

//	report();
*/
	return;
}

void gtf_file::set_expression_weight()
{

	for(auto ci = GENES.begin(); ci != GENES.end(); ++ci)
		for(auto gi = ci->second.begin(); gi != ci->second.end(); ++gi)
			gi->second.set_expression_weight();
	return;
}

void gtf_file::set_expression()
{
	if(EXPRESSION.empty())
	{
		cerr << "\n\nWarning: NO EXPRESSION DATA AVAILABLE!" << endl;
		return;
	}

	for(auto mi = ID_TO_TR_MAP.begin(); mi != ID_TO_TR_MAP.end(); ++mi)
	{
		auto ei = EXPRESSION.find(mi->first);

		if(ei == EXPRESSION.end())
		{
			cerr << "\nWarning: no expression data available for: " << mi->first;
			continue;
		}
		
		mi->second->second.set_expr(ei->second);
	}

	return;
}

void gtf_file::set_expression_gene()
{
	for(auto ci = GENES.begin(); ci != GENES.end(); ++ci)
	{
		for(auto gi = ci->second.begin(); gi != ci->second.end(); ++gi)
		{
			auto xi = GENE_EXPRESSION.find(gi->first);

			if(xi != GENE_EXPRESSION.end())
				gi->second.set_expr(xi->second);
		}
	}

	return;
}

void gtf_file::read_expression(const string &file, map<string, double> &expression)
{
	ifstream in(file.c_str());

	if(!in)
	{
		cerr << "\n\nWarning: Can't find expression values file: " << file << endl;
		return;
	}

	cerr << endl << "Reading expression values file: " << file << endl;

	string line;

	while(getline(in,line))
	{
		if(line.empty())
			continue;
		if(line[0] == '#')
			continue;

		istringstream str(line);
		string id;
		double v;

		str >> id >> v;

		expression.insert(make_pair(id,v));
	}

	in.close();

	return;
}

void gtf_file::read_tr_to_gene_file(const string *File)
{
	if(File->empty())
		return;	

	ifstream in(File->c_str());

	if(!in)
	{
		cerr << "\nWarning: can't find name table file: " << *File << endl;
		return;
	}
	else
		cerr << "\nReading transcript/gene correspondence file: " << *File << endl;

	string line;

	while(getline(in,line))
	{
		if(line.empty())
			continue;
		if(line[0] == '#')
			continue;

		istringstream str(line);
		string s1,s2;

		str >> s1 >> s2;

		TR_TO_GENE[s1] = s2;
		
	}

	in.close();

	return;
}

const map<string, map<string, gene> >* gtf_file::get_genes_ptr() const
{
	return &GENES;
}

map<string, map<string, gene> >* gtf_file::get_genes_ptr() 
{
        return &GENES;
}

void gtf_file::clean_ids(string *s)
{
	s->erase(0,1);	 // REMOVE DOUBLE QUOTES  

	s->erase(s->size() - 2, 2);

	return;
}

const string * gtf_file::get_gene_name(const string *tr_id) const
{
	map<string,string>::const_iterator mi = TR_TO_GENE.find(*tr_id);

	if(mi != TR_TO_GENE.end())
		return &mi->second;
	else
		return NULL;
}

void gtf_file::report() const
{
	map<string, map<string, gene> >::const_iterator mi = GENES.cbegin();

	while(mi != GENES.cend())
	{
		map<string, gene>::const_iterator gi = mi->second.cbegin();

		while(gi != mi->second.cend())
		{
			gi->second.report();
			gi++;
		}

		mi++;
	}	

	return;
}

// END OF GTF_FILE CLASS

// GENE CLASS
gene::gene(const string *chr, const string *gene_id, const char strand, const RefVector *refv): refv(refv), gene_id(*gene_id), chr(*chr), strand(strand), tot_expr(0), expression(0), min_pv(1), max_pv(1), tot_pv(1), omo_count(0), etero_count(0)
{
	set_chr_ref();

	return;
}

const double gene::get_min_pv() const
{
	return min_pv;
}

const double gene::get_max_pv() const
{
        return max_pv;
}

const double gene::get_tot_pv() const
{
        return tot_pv;
}

const string gene::id() const
{
	return gene_id;
}

const double gene::get_expr() const
{
	return expression;
}

void gene::set_expr(double xpr)
{
	expression = xpr;	

	return;
}

const double gene::gtest_1() const
{
	double att, big, small, tmp;

	big = this->get_big_c();
	small = this->get_small_c();

	if(small > big)     //NECESSARY FOR PHASED ANALYSIS
		swap(big, small);

	att = (big + small) / 2.0;

        if(small == 0)
        	tmp = (big * log(big / att));
        else
        	tmp = (big * log(big / att)) + (small * log(small / att));

        tmp *= 2;

        return gsl_cdf_chisq_Q(tmp, 1);

}

void gene::gtest_2() 
{

	double tmp_sum = 0;
	vector<double> v_pv;

	for(unsigned int i = 0; i < SNP_C.size(); i++)
	{
        	double att, big, small, tmp;

        	big = SNP_C.at(i).big_c;
        	small = SNP_C.at(i).small_c;

		if(small > big)    //NECESSARY FOR PHASED ANALYSIS
			swap(big,small);

        	att = (big + small) / 2.0;

        	if(small == 0)
                	tmp = (big * log(big / att));
        	else
                	tmp = (big * log(big / att)) + (small * log(small / att));

        	tmp *= 2;

		v_pv.push_back(gsl_cdf_chisq_Q(tmp, 1));

		tmp_sum += tmp;
	}

	tot_pv = gsl_cdf_chisq_Q(tmp_sum, SNP_C.size());
	min_pv = *min_element(v_pv.begin(), v_pv.end());
	max_pv = *max_element(v_pv.begin(), v_pv.end());

        return;

}

void gene::add_snp_c(unsigned int ref, unsigned int alt)
{
//	if(ref == 0 && alt == 0)
//		return;

	snp_counts snp_c;

	if(ref >= alt)
	{
		snp_c.big_c = ref;
		snp_c.small_c = alt;
	}
	else
	{       
                snp_c.big_c = alt;
                snp_c.small_c = ref;
        }

	SNP_C.push_back(snp_c);

	return;
}

void gene::add_snp_c(unsigned int ref, unsigned int alt, short int phase)
{
	snp_counts snp_c;

	if(phase < 0)
	{
		cerr << endl << "ERROR: counting UNPHASED snp as PHASED" << endl;
		exit(EXIT_FAILURE);
	}

	if(phase == 1)
	{
                snp_c.big_c = ref;   //big_c acts as first chr and small_c acts as second chr of the chr couple
                snp_c.small_c = alt;
        }
	else if(phase == 0)
	{
		snp_c.big_c = alt;  //big_c acts as first chr and small_c acts as second chr of the chr couple
                snp_c.small_c = ref;
	}
	else
	{
		cerr << endl << "ERROR: unknown PHASE" << endl;
		exit(EXIT_FAILURE);
	}

	SNP_C.push_back(snp_c);
	
	return;

}

const unsigned int gene::get_big_c() const
{
	unsigned int count = 0;

	for(unsigned int i = 0; i < SNP_C.size(); i++)
		count += SNP_C.at(i).big_c;

	return count;
}

const unsigned int gene::get_small_c() const
{
        unsigned int count = 0;

        for(unsigned int i = 0; i < SNP_C.size(); i++)
                count += SNP_C.at(i).small_c;

        return count;
}

void gene::add_etero_var_to_count() const
{
	etero_count++;
	return;
}

void gene::add_omo_var_to_count() const
{
        omo_count++;
        return;
}

const string gene::report_var_counts() const
{
	ostringstream str;

	str << omo_count << '\t' << etero_count;

	return str.str();
}

const bool gene::has_transcript_with_var() const
{
	for(auto vi = TRANSCRIPTS.cbegin(); vi != TRANSCRIPTS.cend(); ++vi)
		if(vi->second.get_count_var() > 0)
			return true;	

	return false;
}

const int gene::get_chr_ref() const
{
	return chr_ref;
}

const string gene::get_chr() const
{
	return chr;
}

const map<string, transcript> * gene::get_transcripts_ptr() const
{
	return &TRANSCRIPTS;
}

map<string, transcript> * gene::get_transcripts_ptr() 
{
	return &TRANSCRIPTS;
}

map<string, transcript>::iterator gene::add_transcript_or_exon(const string *tr_id, unsigned int start, unsigned int end)
{
	map<string, transcript>::iterator mi = TRANSCRIPTS.find(*tr_id);

	if(mi == TRANSCRIPTS.end())
		mi = TRANSCRIPTS.emplace(make_pair(*tr_id, transcript(&gene_id, tr_id))).first;

	mi->second.add_exon(start, end);

	return mi;
}

const string gene::report() const
{
	ostringstream out;

	out << gene_id << "\tN_TR = " << TRANSCRIPTS.size() << '\t';
	
	map<string, transcript>::const_iterator mi = TRANSCRIPTS.cbegin();

	while(mi != TRANSCRIPTS.cend())
	{
//		mi->second.report();
		out << mi->first << ',';
		mi++;
	}

	return out.str();
}

void gene::set_chr_ref()
{
        static set<string> UNUSED_CHR;
	chr_ref = INVALID_CHR;

        for(unsigned int i = 0; i < refv->size(); i++)
                if(refv->at(i).RefName == chr)
                {
                        chr_ref = (int)i;
                        return;
                }

        if(UNUSED_CHR.find(chr) == UNUSED_CHR.end())
        {
                cerr << "\nWarning, this reference sequence found in GTF file does not exist in the rnaseq file: " << chr;
                cerr << "\n-> All occurrences of " << chr << " will be ignored." << endl;
                UNUSED_CHR.insert(chr);
        }

        return;
}

void gene::set_expression_weight()
{
	for(auto ti = TRANSCRIPTS.begin(); ti != TRANSCRIPTS.end(); ++ti)
		tot_expr += ti->second.get_expr();

	for(auto ti = TRANSCRIPTS.begin(); ti != TRANSCRIPTS.end(); ++ti)
		ti->second.set_expr_weight(tot_expr);

	return;
}
// END OF GENE CLASS

//TRANSCRIPT CLASS

transcript::transcript(const string *gene_id, const string *tr_id): tr_id(*tr_id), gene_id(*gene_id), expr(NO_EXPRESSION), expr_w(NO_EXPRESSION), count_var(0), aspec_flag(false)
{
	return;
}

const string transcript::snp_count_report() const
{
	ostringstream str;

	unsigned int omo = 0, etero = 0;

	for(unsigned int i = 0; i < POS_VAR_V.size(); i++)
		if(POS_VAR_V.at(i).T == variant::ALT_OMO_Z)
			omo++;
		else if(POS_VAR_V.at(i).T == variant::ETERO_Z)
			etero++;

	str << omo << '\t' << etero;

	return str.str();
}

const string transcript::all_prob_report(unsigned int t) const
{
	ostringstream str;

//	str << "size = " << V_prob[t].size() << '\t';

	for(unsigned int i = 0; i < V_prob[t].size(); i++)
		str << V_prob[t].at(i) << ',';

	return str.str();
}

const double transcript::get_prod_prob(unsigned int t) const
{
	if(V_prob[t].empty())
		return 0;

	double P = 1;

	for(unsigned int i = 0; i < V_prob[t].size(); i++)
		P *= V_prob[t].at(i);

	return P;
}

void transcript::add_pos_prob(double *P) const
{
	for(unsigned int i = 0; i < 2; i++)
	{
		if(P[i] < 0 || P[i] > 1)
		{
			cerr << "\nSomething very weird occurring with probability values! Getting a prob of " << P[i] << " !!! Transcript: " << tr_id << '[' << i << ']'  << endl;
	//		exit(EXIT_FAILURE);
		}

		V_prob[i].push_back(P[i]);
	}

	return;
}

void transcript::add_var_count(bool as_flag) const
{
	count_var++;
	aspec_flag += as_flag;
	
	return;
}

const bool transcript::get_aspec_flag() const
{
	return aspec_flag;
}

const unsigned int transcript::get_count_var() const
{
	return count_var;
}

const double transcript::get_expr() const
{
	return expr;
}

void transcript::set_expr(double xp)
{
	expr = xp;

	return;
}

const vector<exon>* transcript::get_exons_ptr() const
{
	return &EXONS;
}

void transcript::add_exon(const unsigned int start, const unsigned int end)
{
	EXONS.emplace_back(exon(start, end));

	return;
}

const unsigned int transcript::get_pos_ptrv_size() const
{
	return POS_VAR_V.size();
}

const unsigned int transcript::get_pos_ptrv_chr_ref_at(unsigned int i) const
{
	return POS_VAR_V.at(i).chr_ref;
}

const unsigned int transcript::get_pos_ptrv_pos_at(unsigned int i) const
{
	return POS_VAR_V.at(i).pos;
}

const string transcript::report() const
{
	ostringstream out;

	out <<	endl << tr_id << "\nN_EX = " << EXONS.size() << endl; 

	for(unsigned int i = 0; i < EXONS.size(); i++)	
		out << EXONS.at(i).report();

	return out.str();
}

const string * transcript::get_id() const
{
	return &tr_id;
}

void transcript::add_pos(unsigned int c, unsigned int p, unsigned int T)
{
	POS_VAR_V.emplace_back(pos_pointer{c, p, T});

	return;
}

void transcript::set_expr_weight(double tot_expr)
{
	if(tot_expr > 0)
		expr_w = expr / tot_expr;

	return;
}

/*const double transcript::get_expr_weight() const
{
	return expr_w;
}*/

const string transcript::get_gene_id() const
{
	return gene_id;
}
//END OF TRANSCRIPT CLASS

//EXON CLASS
exon::exon(const unsigned int start, const unsigned int end): start(start - 1), end(end - 1)  //CONVERT TO 0 BASED COORD SYS
{
	return;
}

const string exon::report() const
{
	ostringstream out;

	out << ' ' << start << ',' << end;

	return out.str();
}

const unsigned int exon::get_start() const
{
	return start;
}

const unsigned int exon::get_end() const
{
        return end;
}

//END OF EXON CLASS
