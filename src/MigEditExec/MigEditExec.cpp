
#include <iostream>
#include <ctime>

#include "MigEditConfig.h"

#include "AudioFileIf.h"
#include "AudioInfo.h"
#include "CommandLineOptions.h"

#define WITH_FLOATEXCEPTIONS
#define WITH_MEMORYCHECK

// include exception header
#if (defined(WITH_FLOATEXCEPTIONS) && !defined(NDEBUG) && defined (GTCMT_WIN32))
#include <float.h>
#endif // #ifndef WITHOUT_EXCEPTIONS

// include memory leak header
#if (defined(WITH_MEMORYCHECK) && !defined(NDEBUG) && defined (GTCMT_WIN32))
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#else
#include <stdlib.h>
#endif

using std::cout;
using std::endl;

// local function declarations
void    showClInfo ();


enum ClOptionIdx_t
{
    kInputName,             //!< input file path
    //kOutputName,            //!< output file path
    kNormalize,             //!< processing: normalize level
    //kDownmix,
    //kResample,
    //kSpectrogram,           //!< compute spectrogram
    //kCorrelation,           //!< compute ACF coefficients
    //kPitchChroma,           //!< compute pitch chroma

    kNumClOptions
};

const CCommandLineOptions::COption MyOptions[kNumClOptions] = 
{
    CCommandLineOptions::COption(kInputName, "-i", "input audio file path", CCommandLineOptions::COption::kString),
    CCommandLineOptions::COption(kNormalize, "-n", "normalize the input file", CCommandLineOptions::COption::kBool),
};

/////////////////////////////////////////////////////////////////////////////////
// main function
int main(int argc, char* argv[])
{
    std::string             sInputFilePath,                 //!< file paths
                            sOutputFilePath;

    CAudioFileIf            *phOutputFile       = 0;        //!< output audio file
    CAudioFileIf            *phInputFile        = 0;        //!< input audio file

    CAudioFileIf::FileSpec_t stFileSpec;                    //!< input/output file spec

    long long               iInFileLength       = 0;        //!< length of input file

    clock_t                 time                = 0;

    CAudioInfo              *phAudioInfo        = 0;
    float                   **ppfAudioData      = 0;

    CCommandLineOptions     *phClArgs           = 0;

    // detect memory leaks in win32
#if (defined(WITH_MEMORYCHECK) && !defined(NDEBUG) && defined (GTCMT_WIN32))
    // set memory checking flags
    int iDbgFlag = _CrtSetDbgFlag(_CRTDBG_REPORT_FLAG);
    iDbgFlag       |= _CRTDBG_CHECK_ALWAYS_DF;
    iDbgFlag       |= _CRTDBG_LEAK_CHECK_DF;
    _CrtSetDbgFlag( iDbgFlag );
#endif

    // enable floating point exceptions in win32
#if (defined(WITH_FLOATEXCEPTIONS) && !defined(NDEBUG) && defined (GTCMT_WIN32))
    // enable check for exceptions (don't forget to enable stop in MSVC!)
    _controlfp(~(_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW | _EM_UNDERFLOW | _EM_DENORMAL), _MCW_EM) ;
#endif // #ifndef WITHOUT_EXCEPTIONS

    CCommandLineOptions::create(phClArgs);
    CAudioFileIf::create(phInputFile);
    CAudioInfo::create(phAudioInfo);

    showClInfo ();

    //////////////////////////////////////////////////////////////////////////////
    // parse command line arguments
    phClArgs->init(MyOptions, kNumClOptions);
    phClArgs->process(argc, argv);

    if (phClArgs->isOptionSet(kInputName))
        phClArgs->getOption(kInputName, sInputFilePath);
    else
    {
        cout << "input path not set";
        return -1;
    }

    //////////////////////////////////////////////////////////////////////////////
    // open the input wave file
    phInputFile->openFile(sInputFilePath, CAudioFileIf::kFileRead);
    if (!phInputFile->isOpen())
    {
        cout << "Input wave file open error!";
        return -1;
    }
    phInputFile->getFileSpec(stFileSpec);
    phInputFile->getLength(iInFileLength);


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

    // get audio info
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


    //////////////////////////////////////
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

