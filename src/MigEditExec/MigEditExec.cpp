
#include <iostream>
#include <ctime>

#include "MigEditConfig.h"

#include "Util.h"
#include "AudioFileIf.h"

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

void    getClArgs (std::string &sInputFilePath, std::string &sOutputFilePath, int argc, char* argv[]);

/////////////////////////////////////////////////////////////////////////////////
// main function
int main(int argc, char* argv[])
{
    std::string             sInputFilePath,                 //!< file paths
                            sOutputFilePath;

    CAudioFileIf            *phOutputFile       = 0;        //!< output audio file
    CAudioFileIf            *phInputFile        = 0;        //!< input audio file

    CAudioFileIf::FileSpec_t stFileSpec;                    //!< input/output file spec

    int                     iLengthOfBlock      = 0;        //!< length of one input block to process

    long long               iInFileLength       = 0;        //!< length of input file

    clock_t                 time                = 0;

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

    showClInfo ();

    // parse command line arguments
    getClArgs (sInputFilePath, sOutputFilePath, argc, argv);

    //////////////////////////////////////////////////////////////////////////////
    // open files
    {

        // open the input wave file
        CAudioFileIf::createInstance(phInputFile);
        phInputFile->openFile(sInputFilePath, CAudioFileIf::kFileRead);
        if (!phInputFile->isOpen())
        {
            cout << "Input wave file open error!";
            return -1;
        }
        phInputFile->getFileSpec(stFileSpec);
        phInputFile->getLength(iInFileLength);

        // open the output wave file
        CAudioFileIf::createInstance(phOutputFile);
        phOutputFile->openFile(sOutputFilePath, CAudioFileIf::kFileWrite, &stFileSpec);
        if (!phOutputFile->isOpen())
        {
            cout << "Input wave file open error!";
            return -1;
        }
    }

    // write remaining samples to file
    //phOutputFile->writeData(ppfIrData, static_cast<int>(iLengthOfIr)-1);

    //////////////////////////////////////
    // clean-up
    // close the files
    CAudioFileIf::destroyInstance(phInputFile);
    CAudioFileIf::destroyInstance(phOutputFile);

    return 0;
    
}


void     showClInfo()
{
    cout << "GTCMT MigEdit" << endl;
    cout << "(c) 2015 by Alexander Lerch" << endl;
    cout    << "output" << endl;
    cout  << endl;

    return;
}

void getClArgs( std::string &sInputFilePath, std::string &sOutputFilePath, int argc, char* argv[] )
{
    if (argc > 1)
        sInputFilePath.assign (argv[1]);
    if (argc > 2)
        sOutputFilePath.assign (argv[3]);
}
