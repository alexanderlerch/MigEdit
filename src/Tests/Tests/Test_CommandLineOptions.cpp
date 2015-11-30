#include "MigEditConfig.h"

#ifdef WITH_TESTS
#include <cassert>
#include <cstdio>

#include "UnitTest++.h"

#include "SignalGen.h"

#include "CommandLineOptions.h"

enum ClOptionIdx_t
{
    kInputName,             //!< input file path
    kOutputName,            //!< output file path
    kNormalize,             //!< processing: normalize level
    kDownmix,

    kNumClOptions
};

const CCommandLineOptions::COption MyOptions[kNumClOptions] = 
{
    CCommandLineOptions::COption(kInputName, "-i", "input audio file path", CCommandLineOptions::COption::kString),
    CCommandLineOptions::COption(kOutputName, "-o", "output file path", CCommandLineOptions::COption::kString),
    CCommandLineOptions::COption(kNormalize, "-n", "normalize the input file", CCommandLineOptions::COption::kBool),
    CCommandLineOptions::COption(kDownmix, "-dm", "downmix to one channel", CCommandLineOptions::COption::kBool),
};

SUITE(CommandLineOptions)
{
    struct CommandLineOptionsData
    {
        CommandLineOptionsData() 
        {
            CCommandLineOptions::create(m_phClArgs);
            m_phClArgs->init(MyOptions, kNumClOptions);
            argc = 10;
        }

        ~CommandLineOptionsData() 
        {
            m_phClArgs->reset();
            CCommandLineOptions::destroy(m_phClArgs);

        }

        CCommandLineOptions *m_phClArgs;

        int argc;
        char *argv[10];

    };

    TEST_FIXTURE(CommandLineOptionsData, Nonsense)
    {
        argc = 3;

        argv[1] = "asd;lfkj";
        argv[2] = "d'apiwefnl";

        m_phClArgs->process (argc, argv);

        for (int i = 0; i < kNumClOptions; i++)
        {
            CHECK_EQUAL(false, m_phClArgs->isOptionSet (i));
        }
    }
     TEST_FIXTURE(CommandLineOptionsData, ApiCalls)
    {
        CHECK_EQUAL(kFunctionInvalidArgsError,m_phClArgs->process (argc, 0));
        CHECK_EQUAL(kFunctionInvalidArgsError,m_phClArgs->process (-1, argv));
        CHECK_EQUAL(false, m_phClArgs->isOptionSet (-1));
        CHECK_EQUAL(false, m_phClArgs->isOptionSet (kNumClOptions));
    }

     TEST_FIXTURE(CommandLineOptionsData, RealWorld)
     {
         argc = 7;

         argv[1] = "-i";
         argv[2] = "input";
         argv[3] = "-o";
         argv[4] = "output";
         argv[5] = "-n";
         argv[6] = "-dm";

         m_phClArgs->process (argc, argv);

         for (int i = 0; i < kNumClOptions; i++)
         {
             CHECK_EQUAL(true, m_phClArgs->isOptionSet (i));
         }
         std::string test;
         m_phClArgs->getOption (kInputName,test);
         CHECK_EQUAL("input", test);
         m_phClArgs->getOption (kOutputName,test);
         CHECK_EQUAL("output", test);
         m_phClArgs->getOption (kNormalize,test);
         CHECK_EQUAL("", test);
     }
}

#endif //WITH_TESTS