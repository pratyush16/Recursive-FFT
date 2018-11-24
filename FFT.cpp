//Pratyush Potu
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<string>
#include<math.h>
#include<complex>
using namespace std;

#define BYTE unsigned char

typedef struct  WAV_HEADER{
    char       RIFF[4];        // RIFF Header
    BYTE       FileLength[4];   // Total Length of File - 8  
    char       WAVE[4];        // WAVE Header      
    char       fmt[4];         // FMT header       
    BYTE       ST[4];          // Sixteen                                
    BYTE       AudioFormat[2];    // Audio format  
    BYTE       Stereo[2];      // Stereo                   
    BYTE       SamplesPerSec[4];  // Data points per second                             
    BYTE       BytesPerSec[4];    // Bytes per second 
    BYTE       BytesPS[2];     // Number of bytes per sample 
    BYTE       BitsPS[2];  // Number of bits per sample      
    char       data[4]; // "data"  string   
    BYTE       dataSize[4];  // Sampled data length    
}wav_hdr;

wav_hdr wavHeader;

void function(ifstream& in, ofstream& out);
void FFT(complex<double> In[], complex<double> Out[], int N, bool fInv);

int main()
{
    //stream declarations
    ifstream in_stream;
    ofstream out_stream;

    in_stream.open("file5.wav");
    
    if (in_stream.fail())
    {
        cout << "Input file opening failed. \n";
        exit(EXIT_FAILURE); //exit if cannot open file
    }
    
    out_stream.open("file5_edited.wav"); //connect to the output file and test
    if (out_stream.fail())
    {
        cout << "Output file opening failed. \n";
        exit(EXIT_FAILURE); //exit if cannot open file
    }

    function(in_stream, out_stream);

    in_stream.close();
    out_stream.close();
    return 0;
}

void function(ifstream& in, ofstream& out)
{
    for(int i=0;i<4;i++)
    {
        in >> wavHeader.RIFF[i];
        out << wavHeader.RIFF[i];
    }
    for(int i=0;i<4;i++)
    {
        in >> wavHeader.FileLength[i];
        out << wavHeader.FileLength[i];
    }
    for(int i=0;i<4;i++)
    {
        in >> wavHeader.WAVE[i];
        out << wavHeader.WAVE[i];
    }
    for(int i=0;i<4;i++)
    {
	in.get(wavHeader.fmt[i]);
	out.put(wavHeader.fmt[i]);
    }
    for(int i=0;i<4;i++)
    {
        in >> wavHeader.ST[i];
        out << wavHeader.ST[i];
    }
    for(int i=0;i<2;i++)
    {
        in >> wavHeader.AudioFormat[i];
        out << wavHeader.AudioFormat[i];
    }
    for(int i=0;i<2;i++)
    {
        in >> wavHeader.Stereo[i];
        out << wavHeader.Stereo[i];
    }
    for(int i=0;i<4;i++)
    {
        in >> wavHeader.SamplesPerSec[i];
        out << wavHeader.SamplesPerSec[i];
    }
    for(int i=0;i<4;i++)
    {
        in >> wavHeader.BytesPerSec[i];
        out << wavHeader.BytesPerSec[i];
    }
    for(int i=0;i<2;i++)
    {
        in >> wavHeader.BytesPS[i];
        out << wavHeader.BytesPS[i];
    }
    for(int i=0;i<2;i++)
    {
        in >> wavHeader.BitsPS[i];
        out << wavHeader.BitsPS[i];
    }
    for(int i=0;i<4;i++)
    {
        in >> wavHeader.data[i];
        out << wavHeader.data[i];
    }
    for(int i=0;i<4;i++)
    {
        in >> wavHeader.dataSize[i];
        out << wavHeader.dataSize[i];
    }

    int n = 0;
    for (int i=0;i<4;i++)
    {
	n += wavHeader.dataSize[i] * pow(256, i); //size of data
    }

    if (ceil(log2(n)) != floor(log2(n)))
    {
	n = pow(2, ceil(log2(n)));
    }
    
    BYTE *dataArray = new BYTE[n];
    for (int i=0;i<n;i++)
    {
    	dataArray[i] = 0;
    }
    for (int i=0;i<n;i++)
    {
        in.get(reinterpret_cast <char&> (dataArray[i]));
    }

    complex<double> *dataArrayC=new complex<double>[n];
    for(int i=0;i<n;i++)
    {
        dataArrayC[i].real() = dataArray[i];
	dataArrayC[i].imag() = 0.0;
    }

    complex<double> *FFTData=new complex<double>[n] ;//FFT'd array
    complex<double> *finalData=new complex<double>[n];//FFT'd back array

    FFT(dataArrayC, FFTData, n, false);
    for(int i=0;i<n;i++)
    {
	FFTData[i] /= sqrt(n);
    }

    double max = 0;
    for(int i=1;i<n;i++)
    {
	if(abs(FFTData[i]) > max)
	{
	    max = abs(FFTData[i]);
	}
    }

    double val;

    cout << "Enter Cutoff Percentage: ";
    cin >> val;
    val /= 100;

    for(int i=0;i<n;i++)
    {
	if (abs(FFTData[i]) <= val*max)
	{
	    FFTData[i].real() = 0.0;
	    FFTData[i].imag() = 0.0;
	}
    }

    FFT(FFTData, finalData, n, true);
    for(int i=0;i<n;i++)
    {
        finalData[i] /= sqrt(n);
    }

    for(int i=0;i<n;i++)
    {
/*	if (finalData[i].real() < 0)
	{
	    finalData[i].real() = 0;
	}
	else if (finalData[i].real() > 255)
	{
	    finalData[i].real() = 255;
	}*/
	out.put(finalData[i].real());
    }

    delete []finalData;
    delete []FFTData;
    delete []dataArrayC;

    return;
}

void FFT(complex<double> In[], complex<double> Out[], int N, bool fInv)
{
    complex<double> im(0.0, 1.0);
    complex<double> *NewIn= new complex<double>[N/2];
    complex<double> *Even=new complex<double>[N/2];
    complex<double> *Odd= new complex<double>[N/2];

    if (N == 2)
    {
	Out[0] = In[0]+In[1];
	Out[1] = In[0]-In[1];
    }
    else
    {
	for (int i=0;i<(N/2);i++)
	{
	    NewIn[i] = In[2*i];
	}
	FFT(NewIn, Even, N/2, fInv);
	for (int i=0;i<(N/2);i++)
	{
	    NewIn[i] = In[2*i+1];
	}
	FFT(NewIn, Odd, N/2, fInv);
        if(fInv == true)
        {
            for (int n=0;n<(N/2);n++)
            {
                Out[n] = Even[n] + Odd[n] * pow(M_E, ((2*M_PI*im*static_cast<double>(n))/static_cast<double>(N)));
            }
            for (int n=(N/2);n<N;n++)
            {
                Out[n] = Even[n-(N/2)] + Odd[n-(N/2)] * pow(M_E, ((2*M_PI*im*static_cast<double>(n))/static_cast<double>(N)));
            }
        }
        else
        {
            for (int n=0;n<(N/2);n++)
            {
                Out[n] = Even[n] + Odd[n] * pow(M_E, (((-2)*M_PI*im*static_cast<double>(n))/static_cast<double>(N)));
            }
            for (int n=(N/2);n<N;n++)
            {
                Out[n] = Even[n-(N/2)] + Odd[n-(N/2)] * pow(M_E, (((-2)*M_PI*im*static_cast<double>(n))/static_cast<double>(N)));
            }
        }
    }

    delete []NewIn;
    delete []Even;
    delete []Odd;

    return;
}
