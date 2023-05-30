#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

using namespace std;
using namespace std::tr1;

/*
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * Date: Dec 2019
*/

typedef struct {
	unsigned int lineNum;
	unsigned int score;
} fraghit;

typedef struct {
	unsigned int lineNum;
	unsigned int unique;
	unsigned int dup;
	unsigned int discard;
} stat;

const unsigned int MAX_INSERT_SIZE = 1024;
const unsigned int MAX_FILE_NAME   = 256;
const uint64_t FLAG_CRICK          = 0xffULL << 32;

int get_readLen_from_cigar( const string &cigar );
void rmdup_pe(const char *ifile, const char *ofile, unordered_map<string, unordered_map<uint64_t, fraghit> *> &samRecord, stat &s );
void rmdup_se(const char *ifile, const char *ofile, unordered_map<string, unordered_map<uint64_t, fraghit> *> &samRecord, stat &s );

int main( int argc, char *argv[] ) {
	if( argc != 5 ) {
		cerr << "\nUsage: " << argv[0] << " <chr.info> <mode=PE|SE> <in.sam> <out.prefix>\n\n"
			 << "Remove duplicates that have the same start (and end for PE data),\n"
			 << "read with the BEST sequencing quality will be kept.\n"
			 << "Paired-end reads MUST be stored in two adjacent lines in SAM file.\n"
			 << "Note that I will read the 'in.sam' twice, therefore it could not be stdin.\n\n";
		return 2;
	}
	
	bool pe = false;
	if( strcmp(argv[2], "PE")==0 || strcmp(argv[2], "pe")==0 ) {
		pe = true;
	} else if( strcmp(argv[2], "SE")!=0 && strcmp(argv[2], "se")!=0 ) {
		cerr << "ERROR: Only PE or SE modes are acceptable!\n";
		exit(10);
	}

	// loading info file
	cerr << "Loading genome.info ...\n";
	unordered_map<string, unordered_map<uint64_t, fraghit> *> samRecord;
	string line, chr;
	stringstream ss;
	ifstream fin( argv[1] );
	if( fin.fail() ) {
		cerr << "Error: could not open file '" << argv[1] << "'!\n";
		exit( 1 );
	}
	while( true ) {
		getline( fin, line );
		if( fin.eof() )break;
		if( line[0] == '#' )continue;

		ss.str( line );
		ss.clear();
		ss >> chr;
//		cerr << "Adding " << chr << " ...\n";
		samRecord.insert( pair<string, unordered_map<uint64_t, fraghit>*>(chr, new unordered_map<uint64_t, fraghit>()) );
	}
	fin.close();
	
	stat samStat;
	
//	sprintf( outfile, "%s.sam", argv[4] );
	if( pe ) {
		rmdup_pe( argv[3], "/dev/stdout", samRecord, samStat );
	} else {
		rmdup_se( argv[3], "/dev/stdout", samRecord, samStat );
	}
	
	cerr << "Writing log file ...\n";
	// write log
	char * outfile = new char [ MAX_FILE_NAME ];
	sprintf( outfile, "%s.log", argv[4] );
	ofstream fout( outfile );
	if( fout.fail() ) {
		cerr << "Error: could not write output log file!\n";
		exit( 1 );
	}
	fout << "All\t"       << samStat.lineNum	<<  "\t100\n"
		 << "Unique\t"    << samStat.unique		<< '\t' << samStat.unique*100.0/samStat.lineNum	<< '\n'
		 << "Duplicate\t" << samStat.dup		<< '\t' << samStat.dup*100.0/samStat.lineNum	<< '\n'
		 << "Discard\t"   << samStat.discard	<< '\t' << samStat.discard*100.0/samStat.lineNum<< '\n';
	fout.close();

//	cerr << "Freeing memory ...\n";
	delete [] outfile;
}

// deal PE data
void rmdup_pe( const char *ifile, const char *ofile,
				unordered_map<string, unordered_map<uint64_t, fraghit> *> &samRecord, stat &samStat ) {
	// load sam file
	ifstream fin( ifile );
	if( fin.fail() ) {
		cerr << "Error: could not open file '" << ifile << "'!\n";
		exit( 11 );
	}
	unordered_map<string, unordered_map<uint64_t, fraghit> *> :: iterator sam_it;
	unordered_map<string, unordered_map<uint64_t, fraghit> *> :: iterator no_such_chr = samRecord.end();
	unordered_set<unsigned int> dup;
	unordered_set<unsigned int> discard;
	unordered_map<uint64_t, fraghit> :: iterator hit_it;
	fraghit hit;

	string line, line2, chr, cigar1, cigar2, qual1, qual2, unk;
	unsigned int pos1, pos2, fragSize;
	uint64_t key;
	stringstream ss;

	unsigned int lineNum = 0;
	cerr << "Loading sam file ...\n";
	while( true ) {
		getline( fin, line );
		if( fin.eof() )break;
		getline( fin, line2 );
		++ lineNum;

//		if( ! (lineNum & 0x3fffff) ) cerr << '\r' << lineNum << " reads loaded.";
		//499780R1	83	chrX	14710827	42	67M	*	0	0	TCCCAATTCTAAATAGTT	HHHHHHHHHHHHHHHHH XG:Z:GA
		//499780R2	163	chrX	14710378	42	67M	*	0	0	TAACATTTCTTTAATCAC	HHHHHHHH:;:987665 XG:Z:GA

		ss.str( line );
		ss.clear();
		ss >> unk >> unk >> chr >> pos1 >> unk >> cigar1 >> unk >> unk >> unk >> unk >> qual1;
		ss.str( line2 );
		ss.clear();
		ss >> unk >> unk >> chr >> pos2 >> unk >> cigar2 >> unk >> unk >> unk >> unk >> qual2;

		sam_it = samRecord.find( chr );
		if( sam_it == no_such_chr ) {
			discard.insert( lineNum );
			continue;
		}

		if( pos1 < pos2 ) {	// locate the end using pos2 and cigar2
			key = pos1;
			key <<= 32;
			int readLen = get_readLen_from_cigar( cigar2 );
			fragSize = pos2 + readLen - pos1;
			key |= fragSize;
		} else {	//pos1 > pos2, then locate the end using pos1 and cigar1
			key = pos2;
			key <<= 32;
			int readLen = get_readLen_from_cigar( cigar1 );
			fragSize = pos1 + readLen - pos2;
			key |= fragSize;
		}

		// calculate the read quality scores
		register unsigned int score = 0;
		const char *p = qual1.c_str();
		register unsigned int i = 0;
		register unsigned int readLen = qual1.size();
		for( i=0; i != readLen; ++i ) {
			score += p[i];
		}
		p = qual2.c_str();
		readLen = qual2.size();
		for( i=0; i != readLen; ++i ) {
			score += p[i];
		}
//		fprintf( stderr, "Line %d, chr=%s, pos1=%u (0x%x), cigar1=%s, pos2=%u (0x%x), cigar2=%s, score=%d\n",
//					lineNum, chr.c_str(), pos1, pos1, cigar1.c_str(), pos2, pos2, cigar2.c_str(), score );

		hit_it = sam_it->second->find( key );
		if( hit_it != sam_it->second->end() ) {	// there must be a duplicate
			if( hit_it->second.score >= score ) {	// the previous one is better, mark this one as duplicate
				dup.insert( lineNum );
//				cerr << "Line " << lineNum << " meet a duplicate at line " << hit_it->second.lineNum
//					 << " and the other one is better.\n";
			} else {	// this one is better, then mark the previous one as duplicate and update the record
//				cerr << "Line " << lineNum << " meet a duplicate at line " << hit_it->second.lineNum
//					 << " but this one is better.\n";
				dup.insert( hit_it->second.lineNum );
				hit_it->second.lineNum = lineNum;
				hit_it->second.score = score;
			}
		} else {	// no such record, add this one
			hit.lineNum = lineNum;
			hit.score = score;
			sam_it->second->insert( pair<uint64_t, fraghit>( key, hit ) );
//			cerr << "Line " << lineNum << " is new.\n";
		}
	}
	cerr << "\rDone: " << lineNum << " reads loaded, dup=" << dup.size() << ", discard=" << discard.size() << ".\n";

	cerr << "Writing output ...\n";
	// prepare output file
	ofstream fout( ofile );
	if( fout.fail() ) {
		cerr << "Error: could not write output SAM file!\n";
		exit( 1 );
	}

	// rewind sam file
	fin.clear();
	fin.seekg( ios_base::beg );
	lineNum = 0;
	unordered_set<unsigned int> :: iterator non_dup = dup.end();
	unordered_set<unsigned int> :: iterator non_discard = discard.end();
	unsigned int unique=0;
	while( true ) {
		getline( fin, line );
		if( fin.eof() )break;
		getline( fin, line2 );
		++ lineNum;
//		if( ! (lineNum & 0x3fffff) ) cerr << '\r' << lineNum << " reads loaded.";

		if( discard.find(lineNum)==non_discard && dup.find(lineNum)==non_dup ) {
			fout << line << '\n' << line2 << '\n';
			++ unique;
		} else {
//			cerr << "Line " << lineNum << " is marked as duplicate.\n";
		}
	}
	fin.close();
	fout.close();
	cerr << "\rDone: " << unique << " reads written.\n";

	samStat.lineNum = lineNum;
	samStat.unique  = unique;
	samStat.dup		= dup.size();
	samStat.discard = discard.size();

	for( sam_it=samRecord.begin(); sam_it!=no_such_chr; ++sam_it ) {
		delete sam_it->second;
	}
}

// deal SE data
void rmdup_se( const char *ifile, const char *ofile,
				unordered_map<string, unordered_map<uint64_t, fraghit> *> &samRecord, stat &samStat ) {
	// load sam file
	ifstream fin( ifile );
	if( fin.fail() ) {
		cerr << "Error: could not open file '" << ifile << "'!\n";
		exit( 1 );
	}
	unordered_map<string, unordered_map<uint64_t, fraghit> *> :: iterator sam_it;
	unordered_map<string, unordered_map<uint64_t, fraghit> *> :: iterator no_such_chr = samRecord.end();
	unordered_set<unsigned int> dup;
	unordered_set<unsigned int> discard;
	unordered_map<uint64_t, fraghit> :: iterator hit_it;
	fraghit hit;

	string line, chr, cigar, qual, unk;
	unsigned int samflag, pos;
	uint64_t key;
	stringstream ss;

	unsigned int lineNum = 0;
	while( true ) {
		getline( fin, line );
		if( fin.eof() ) break;
		++ lineNum;
		if( ! (lineNum & 0x3fffff) ) {
			cerr << '\r' << lineNum << " reads loaded.";
		}

		//499780R1	83	chrX	14710827	42	67M	*	0	0	TCCCAATTCTAAATAGTT	HHHHHHHHHHHHHHHHH XG:Z:GA
		ss.str( line );
		ss.clear();
		ss >> unk >> samflag >> chr >> pos >> unk >> cigar >> unk >> unk >> unk >> unk >> qual;

//		fprintf( stderr, "Line %d, chr=%s, pos1=%u (0x%x), cigar1=%s, pos2=%u (0x%x), cigar2=%s\n",
//					lineNum, chr.c_str(), pos1, pos1, cigar1.c_str(), pos2, pos2, cigar2.c_str() );

		sam_it = samRecord.find( chr );
		if( sam_it == no_such_chr ) {
			discard.insert( lineNum );
			continue;
		}

		if( samflag & 0x10 ) {	//read unordered_mapped to CRICK strand
			pos += get_readLen_from_cigar( cigar );
			key = pos;
			key |= FLAG_CRICK;
		} else {	// read unordered_mapped to WATSON strand
			key = pos;
		}

		register unsigned int score = 0;
		const char * p = qual.c_str();
		register unsigned int i = 0;
		register unsigned int readLen = qual.size();
		for( i=0; i != readLen; ++i ) {
			score += p[i];
		}
//		fprintf( stderr, "line %d => key=0x%llx, score=%u\n",
//					lineNum, key, score );

		hit_it = sam_it->second->find( key );
		if( hit_it != sam_it->second->end() ) {	// there must be a duplicate
			if( hit_it->second.score >= score ) {	// the previous one is better, mark this one as duplicate
				dup.insert( lineNum );
//				cerr << "Line " << lineNum << " meet a duplicate at line " << hit_it->second.lineNum
//					 << " and the other one is better.\n";
			} else {	// this one is better, then mark the previous one as duplicate and update the record
//				cerr << "Line " << lineNum << " meet a duplicate at line " << hit_it->second.lineNum
//					 << " but this one is better.\n";
				dup.insert( hit_it->second.lineNum );
				hit_it->second.lineNum = lineNum;
				hit_it->second.score = score;
			}
		} else {	// no such record, add this one
			hit.lineNum = lineNum;
			hit.score = score;
			sam_it->second->insert( pair<uint64_t, fraghit>( key, hit ) );
//			cerr << "Line " << lineNum << " is new.\n";
		}
	}
	cerr << "\rDone: " << lineNum << " reads loaded, dup=" << dup.size() << ", discard=" << discard.size() << ".\n";

	// prepare output file
	cerr << "Writing output ...\n";
	ofstream fout( ofile );
	if( fout.fail() ) {
		cerr << "Error: could not write output SAM ile!\n";
		exit( 1 );
	}
	// rewind sam file
	fin.clear();
	fin.seekg( ios_base::beg );
	lineNum = 0;
	unordered_set<unsigned int> :: iterator non_dup = dup.end();
	unordered_set<unsigned int> :: iterator non_discard = discard.end();
	unsigned int unique=0;
	while( true ) {
		getline( fin, line );
		if( fin.eof() )break;
		++ lineNum;
		if( ! (lineNum & 0x3fffff) ) {
			cerr << '\r' << lineNum << " reads loaded.";
		}

		//499780R1	83	chrX	14710827	42	67M	*	0	0	TCCCAATTCTAAATAGTT	HHHHHHHHHHHHHHHHH ???
		if( discard.find(lineNum)==non_discard && dup.find(lineNum)==non_dup ) {
			fout << line << '\n';
			++ unique;
		} else {
//			cerr << "Line " << lineNum << " is marked as duplicate.\n";
		}
	}
	fin.close();
	fout.close();
	
	cerr << "\rDone: " << unique << " reads written.\n";
	
	samStat.lineNum = lineNum;
	samStat.unique  = unique;
	samStat.dup		= dup.size();
	samStat.discard = discard.size();

	for( sam_it=samRecord.begin(); sam_it!=no_such_chr; ++sam_it ) {
		delete sam_it->second;
	}
}

// CIGARs could contain I for insertion and D for deletions
int get_readLen_from_cigar( const string &cigar ) {
	register int i, j;
	register int size = 0;
	register int cs = cigar.size();
	register const char * p = cigar.c_str();

	for(i=0, j=0; i!=cs; ++i) {
		if( p[i] <= '9' ) {   // digital
			j *= 10;
			j += p[i] - '0';
		} else {        // MUST be M, I, or D
			if( p[i] == 'M' ) { // match or mismatch, keep
				size += j;
			} else if ( p[i] == 'I' ) { // insertion, ignore
				// do nothing
			} else if ( p[i] == 'D' ) { // deletion, add place holders
				size += j;
			} else {	// unsupported CIGAR element
				return 0;
			}
			j = 0;
		}
	}

	return size;
}


