
#include <iostream>
#include <ctime>

#include "MigEditConfig.h"

#include "AudioFileIf.h"
#include "AudioInfo.h"
#include "CommandLineOptions.h"

#include "Preproc.h"
#include "NmfIf.h"
#include "MatrixRepresentation.h"

using std::cout;
using std::endl;

// local function declarations
void    showClInfo ();


enum ClOptionIdx_t
{
    kInputName,             //!< input file path
    //kOutputName,            //!< output file path
    kNormalize,             //!< processing: normalize level
    kNmf,                   //!< non-negative matrix factorization
    //kDownmix,
    //kResample,
    //kSpectrogram,           //!< compute spectrogram
    //kCorrelation,           //!< compute ACF coefficients
    //kPitchChroma,           //!< compute pitch chroma

    kNumClOptions
};

const CCommandLineOptions::COption MyOptions[kNumClOptions] = 
{
    CCommandLineOptions::COption(kInputName,    "-i",   "input audio file path",                        CCommandLineOptions::COption::kString),
    CCommandLineOptions::COption(kNormalize,    "-n",   "normalize the input file",                     CCommandLineOptions::COption::kBool),
    CCommandLineOptions::COption(kNmf,          "-nmf", "compute non-negative matrix factorization",    CCommandLineOptions::COption::kBool),
};

/////////////////////////////////////////////////////////////////////////////////
// main function
int main(int argc, char* argv[])
{
    std::string             sInputFilePath,                 //!< file paths
                            sOutputFilePath;
    CCommandLineOptions     *phClArgs           = 0;

    CAudioFileIf            *phOutputFile       = 0;        //!< output audio file
    CAudioFileIf            *phInputFile        = 0;        //!< input audio file

    CAudioFileIf::FileSpec_t stFileSpec;                    //!< input/output file spec

    long long               iInFileLength       = 0;        //!< length of input file

    clock_t                 time                = 0;

    CAudioInfo              *phAudioInfo        = 0;
    float                   **ppfAudioData      = 0;


    showClInfo ();

    //////////////////////////////////////////////////////////////////////////////
    // parse command line arguments
    CCommandLineOptions::create(phClArgs);
    phClArgs->init(MyOptions, kNumClOptions);
    phClArgs->process(argc, argv);

    //////////////////////////////////////////////////////////////////////////////
    // open the input wave file
    if (phClArgs->isOptionSet(kInputName))
    {
        phClArgs->getOption(kInputName, sInputFilePath);

        CAudioFileIf::create(phInputFile);
        phInputFile->openFile(sInputFilePath, CAudioFileIf::kFileRead);
        if (!phInputFile->isOpen())
        {
            cout << "Input wave file open error!";
            return -1;
        }
        phInputFile->getFileSpec(stFileSpec);
        phInputFile->getLength(iInFileLength);
    }
    else
    {
        cout << "input path not set!"  << endl;
        return -1;
    }

    //////////////////////////////////////////////////////////////////////////////
    // allocate buffer for the whole file
    if (iInFileLength > 0 && stFileSpec.iNumChannels > 0)
    {
        ppfAudioData  = new float* [stFileSpec.iNumChannels];
        for (int i = 0; i < stFileSpec.iNumChannels; i++)
        {
            ppfAudioData[i]   = new float [iInFileLength];
        }
    }

    // get audio data
    phInputFile->readData(ppfAudioData, iInFileLength);

    //////////////////////////////////////////////////////////////////////////////
    // get audio info and print it to stdout
    CAudioInfo::create(phAudioInfo);
    phAudioInfo->init(stFileSpec.fSampleRateInHz, stFileSpec.iNumChannels);
    phAudioInfo->process(ppfAudioData,iInFileLength);

    for (int i = 0; i < CAudioInfo::kNumInfoTypes; i++)
    {
        cout << phAudioInfo->getResultName(static_cast<CAudioInfo::InfoType_t>(i)) << ":\t";
        for (int c = 0; c < stFileSpec.iNumChannels; c++)
        {
            double dValue;
            phAudioInfo->getResult(dValue, static_cast<CAudioInfo::InfoType_t>(i), c);
            cout << dValue << "\t";
        }
        cout << endl;
    }

    //////////////////////////////////////////////////////////////////////////////
    // compute nmf
    if (phClArgs->isOptionSet(kNmf))
    {
        int iBlockLength    = 2048;
        int iHopLength      = 1024;
        int iRank           = 10;
        CPreproc                *phPreProc          = 0;
        CMatrixRepresentation   *phSpecGramComp     = 0;
        CNmfIf                  *phNmf              = 0;
        CMatrixRepresentationResult hSpecGram;
        CNmfParametrization     hNmfInit;
        CNmfResult              hNmfResult;

        // pre-processing
        CPreproc::create(phPreProc);
        phPreProc->init(phAudioInfo);
        phPreProc->setStepActive(CPreproc::kPpDownmix);
        phPreProc->setStepActive(CPreproc::kPpNormalize);
        phPreProc->process(ppfAudioData, ppfAudioData, static_cast<int>(iInFileLength)); // inplace

        // compute spectrogram
        CMatrixRepresentation::create(phSpecGramComp);
        phSpecGramComp->init(   static_cast<int>(iInFileLength), 
            phAudioInfo->getSampleRate(),
            1,//phAudioInfo->getNumChannels(), 
            iBlockLength, 
            iHopLength, 
            hSpecGram);
        phSpecGramComp->process(ppfAudioData, static_cast<int>(iInFileLength));

        // compute nmf
        CMatrix SpecGram = hSpecGram.getResultPtr()->transpose();
        CNmfIf::create(phNmf);
        hNmfInit.init(phSpecGramComp->getParam(CMatrixRepresentation::kBlockLength)/2+1, iRank);
        phNmf->init(hNmfInit);
        time = clock();
        phNmf->process(&SpecGram, hNmfResult);
        cout << endl << "Elapsed processing time: " << (clock() - time )*1.F/CLOCKS_PER_SEC << endl;
        cout << "#Iterations:\t" << hNmfResult.getNumIterations() << endl;
        cout << "Min. Error: \t" << hNmfResult.getError() << endl;
    }

    //////////////////////////////////////////////////////////////////////////////
    // clean-up
    for (int i = 0; i < stFileSpec.iNumChannels; i++)
    {
        delete [] ppfAudioData[i];
    }
    delete [] ppfAudioData ;

    // close the files
    CAudioFileIf::destroy(phInputFile);
    CCommandLineOptions::destroy(phClArgs);
    CAudioInfo::destroy(phAudioInfo);

    return 0;
    
}


void     showClInfo()
{
    cout << "GTCMT MigEdit" << endl;
    cout << "(c) 2015 by Alexander Lerch" << endl;
    //cout    << "output" << endl;
    cout  << endl;

    return;
}

