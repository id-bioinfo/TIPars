/////g++ test.o -o test -fopenmp -lpthread
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <sstream>
#include <vector>
#include <map>
#include <math.h> 
#include <cstring>
#include <time.h>
#include "omp.h"

using namespace std;

unordered_map<string, string> readMSA(string msafile, int minLength = 1, int maxLength = 1000000)
{
    unordered_map<string, string> msaNameHashMap;
    ifstream in_msafile(msafile);
    if (!in_msafile.is_open())  cout << "Unable to open file" << msafile << endl;
    string seqName, seqline, line;
	int lowerflag = 0;
    while ( getline (in_msafile, line) )
    {
        if(line[0] == '>') 
        {
            if(seqline.length() >= minLength && seqline.length() <= maxLength) 
            {
				97-122
				if(lowerflag == 0)
				{
					if(seqline[0] >= 97 && seqline[0] <= 122) lowerflag = 1;
					else lowerflag = -1;
				}
                if(lowerflag == 1) transform(seqline.begin(),seqline.end(),seqline.begin(),::toupper);
                msaNameHashMap.insert({seqName, seqline});
            }
            else if (seqline.length() > 0)
            {
                cout << "error: " << seqline.length() << "," << seqName << endl;
            }
            seqName = line;
            seqline = "";
        }
        else
        {
            seqline += line;
        }
    }
    if(seqline.length() >= minLength && seqline.length() <= maxLength) 
    {
        if(lowerflag == 1) transform(seqline.begin(),seqline.end(),seqline.begin(),::toupper);
        msaNameHashMap.insert({seqName, seqline});
    }
    else
    {
        cout << "error: " << seqline.length() << "," << seqName << endl;
    }

    in_msafile.close();

    return msaNameHashMap;
}

unordered_map<char, unordered_set<char>> generate_nucleotide_nomenclature()
{
    unordered_map<char, unordered_set<char>> nucleotide_nomenclature;
	//-
    unordered_set<char> tempGap; tempGap.insert('X');
    nucleotide_nomenclature.insert({'-',tempGap});
    //A
    unordered_set<char> tempA; tempA.insert('A');
    nucleotide_nomenclature.insert({'A',tempA});
    //C
    unordered_set<char> tempC; tempC.insert('C');
    nucleotide_nomenclature.insert({'C',tempC});
    //G
    unordered_set<char> tempG; tempG.insert('G');
    nucleotide_nomenclature.insert({'G', tempG});
    //T
    unordered_set<char> tempT; tempT.insert('T');
    nucleotide_nomenclature.insert({'T', tempT});
    //R
    unordered_set<char> tempR; tempR.insert('A'); tempR.insert('G');
    nucleotide_nomenclature.insert({'R', tempR});
    //R
    unordered_set<char> tempY; tempY.insert('C'); tempY.insert('T');
    nucleotide_nomenclature.insert({'Y', tempY});
    //M
    unordered_set<char> tempM; tempM.insert('A'); tempM.insert('C');
    nucleotide_nomenclature.insert({'M', tempM});
    //K
    unordered_set<char> tempK; tempK.insert('G'); tempK.insert('T');
    nucleotide_nomenclature.insert({'K', tempK});
    //S
    unordered_set<char> tempS; tempS.insert('C'); tempS.insert('G');
    nucleotide_nomenclature.insert({'S', tempS});
    //W
    unordered_set<char> tempW; tempW.insert('A'); tempW.insert('T');
    nucleotide_nomenclature.insert({'W', tempW});
    //H
    unordered_set<char> tempH; tempH.insert('A'); tempH.insert('C'); tempH.insert('T');
    nucleotide_nomenclature.insert({'H', tempH});
    //B
    unordered_set<char> tempB; tempB.insert('C'); tempB.insert('G'); tempB.insert('T');
    nucleotide_nomenclature.insert({'B', tempB});
    //V
    unordered_set<char> tempV; tempV.insert('A'); tempV.insert('C'); tempV.insert('G'); 
    nucleotide_nomenclature.insert({'V',tempV});
    //D
    unordered_set<char> tempD; tempD.insert('A'); tempD.insert('G'); tempD.insert('T'); 
    nucleotide_nomenclature.insert({'D',tempD});
    //N
    unordered_set<char> tempN; tempN.insert('A'); tempN.insert('C'); tempN.insert('G'); tempN.insert('T'); 
    nucleotide_nomenclature.insert({'N',tempN});
    
    return nucleotide_nomenclature;   	 
}

int main(int argc,char *argv[])
{
	if(argc != 4)
	{
		cout << "usege: splitEachColumn msafile outdir numberofthreads" << endl;
		return -1;
	}
	
    //int num_threads = omp_get_num_procs();
	int num_threads = atoi(argv[3]);
    omp_set_num_threads(num_threads);
    cout << num_threads << " threads" << endl;

    string msafile = argv[1];

    // Start measuring time
    clock_t start = clock();
    unordered_map<string, string> msa = readMSA(msafile);
    // Stop measuring time and calculate the elapsed time
    clock_t end = clock();
    double elapsed = double(end - start)/CLOCKS_PER_SEC;
    printf("Reading MSA Time: %.3f seconds.\n", elapsed);
    start = end;

	int alnLength = (msa.begin())->second.length();
	cout << "msa size = " << msa.size() << " alignmentLength = " << alnLength << endl;
	
	string outdir = argv[2];
    
    unordered_map<char, unordered_set<char>> nomenclatur_table = generate_nucleotide_nomenclature();
    
    #pragma omp parallel for
    for (int i=0; i < alnLength; ++i)
    { 
        string outfilepath = outdir + "/position_"  + to_string(i) + ".txt"; //Output file creation
		ofstream outfile;
		outfile.open(outfilepath, ios::out); 
		if(!outfile.is_open()) 
		{
			cout << "Unable to open file" << outfilepath << endl;
		}
		else
		{
			outfile << "TAXA\tNUCL" << endl;
			for (auto it =  msa.begin(); it != msa.end(); ++it)
			{
				string seqname = (it->first).substr(1);
				char base = it->second[i];
				if (nomenclatur_table.find(base) != nomenclatur_table.end())
				{
					for(auto it2 = nomenclatur_table[base].begin(); it2 != nomenclatur_table[base].end(); ++it2)
					{
						outfile << seqname << "\t" << *it2 << endl;
					}
				}
			}
			if(i % 100 == 99)
				cout << "processing at columns " << i << endl;
		}
		outfile.close();
    }
    
    end = clock();
    elapsed = double(end - start)/CLOCKS_PER_SEC;
    cout << "Spliting columns Time: " << elapsed << " seconds" << endl;
 
    return 0;
}
