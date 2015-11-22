#include "MigEditConfig.h"

#ifdef WITH_TESTS
#include <cassert>
#include <cstdio>

#include "UnitTest++.h"

#include "SignalGen.h"

#include "AudioInfo.h"

SUITE(AudioInfo)
{
    struct AudioInfoData
    {
        AudioInfoData() : m_fSampleRate(48000)
        {
            m_pfBuff1   = new float [m_iAudioInfoLength];
            m_pfBuff2   = new float [m_iAudioInfoLength];

            CAudioInfo::createInstance(m_pCAudioInfoInstance);
            m_pCAudioInfoInstance->initInstance(m_fSampleRate, 1);
        }

        ~AudioInfoData() 
        {
            m_pCAudioInfoInstance->resetInstance();
            CAudioInfo::destroyInstance(m_pCAudioInfoInstance);

            delete [] m_pfBuff1;
            delete [] m_pfBuff2;
        }

        float *m_pfBuff1;
        float *m_pfBuff2;

        static const int m_iAudioInfoLength = 16384;
        const float m_fSampleRate;

        CAudioInfo *m_pCAudioInfoInstance;
    };

    TEST_FIXTURE(AudioInfoData, SimpleSine)
    {
        double dResult = 0;
        CSignalGen::generateSine(m_pfBuff1, 1000.F, m_fSampleRate, m_iAudioInfoLength, .5F, 0);
        
        m_pCAudioInfoInstance->process(&m_pfBuff1, m_iAudioInfoLength);

        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kMax);
        CHECK_CLOSE(.5, dResult, 1e-3); 
        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kMin);
        CHECK_CLOSE(-.5, dResult, 1e-3); 
        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kAbsMax);
        CHECK_CLOSE(.5, dResult, 1e-3); 
        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kRange);
        CHECK_CLOSE(1.0, dResult, 1e-3); 
        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kMean);
        CHECK_CLOSE(0, dResult, 1e-3); 
        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kRms);
        CHECK_CLOSE(.5/sqrt(2.0), dResult, 1e-3); 
        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kStd);
        CHECK_CLOSE(.5/sqrt(2.0), dResult, 1e-3); 

        for (int i = 0; i < m_iAudioInfoLength; i++)
        {
            m_pfBuff1[i]   -= .6F;
        }

        m_pCAudioInfoInstance->initInstance(m_fSampleRate, 1);
        m_pCAudioInfoInstance->process(&m_pfBuff1, m_iAudioInfoLength);

        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kMax);
        CHECK_CLOSE(-.1, dResult, 1e-3); 
        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kMin);
        CHECK_CLOSE(-1.1, dResult, 1e-3); 
        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kAbsMax);
        CHECK_CLOSE(1.1, dResult, 1e-3); 
        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kRange);
        CHECK_CLOSE(1.0, dResult, 1e-3); 
        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kMean);
        CHECK_CLOSE(-.6, dResult, 1e-3); 
        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kStd);
        CHECK_CLOSE(.5/sqrt(2.0), dResult, 1e-3); 
    }
    TEST_FIXTURE(AudioInfoData, Dc)
    {
        double dResult = 0;
        CSignalGen::generateDc(m_pfBuff1, m_iAudioInfoLength, -.1F);

        m_pCAudioInfoInstance->process(&m_pfBuff1, m_iAudioInfoLength);

        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kMax);
        CHECK_CLOSE(-.1, dResult, 1e-3); 
        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kMin);
        CHECK_CLOSE(-.1, dResult, 1e-3); 
        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kAbsMax);
        CHECK_CLOSE(.1, dResult, 1e-3); 
        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kRange);
        CHECK_CLOSE(.0, dResult, 1e-3); 
        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kMean);
        CHECK_CLOSE(-.1, dResult, 1e-3); 
        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kRms);
        CHECK_CLOSE(.1, dResult, 1e-3); 
        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kStd);
        CHECK_CLOSE(.0, dResult, 1e-3); 
    }
    TEST_FIXTURE(AudioInfoData, ZeroInput)
    {
        double dResult = 0;
        CUtil::setZero(m_pfBuff1, m_iAudioInfoLength);

        m_pCAudioInfoInstance->process(&m_pfBuff1, m_iAudioInfoLength);

        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kMax);
        CHECK_CLOSE(.0, dResult, 1e-3); 
        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kMin);
        CHECK_CLOSE(.0, dResult, 1e-3); 
        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kAbsMax);
        CHECK_CLOSE(.0, dResult, 1e-3); 
        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kRange);
        CHECK_CLOSE(.0, dResult, 1e-3); 
        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kMean);
        CHECK_CLOSE(.0, dResult, 1e-3); 
        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kRms);
        CHECK_CLOSE(.0, dResult, 1e-3); 
        m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kStd);
        CHECK_CLOSE(.0, dResult, 1e-3); 
    }
    TEST_FIXTURE(AudioInfoData, ApiCalls)
    {
        double dResult = 0;
        CUtil::setZero(m_pfBuff1, m_iAudioInfoLength);

        m_pCAudioInfoInstance->resetInstance();

        CHECK_EQUAL(kNotInitializedError, m_pCAudioInfoInstance->process(&m_pfBuff1, m_iAudioInfoLength));
        CHECK_EQUAL(kNotInitializedError, m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kMax));
        CHECK_EQUAL(kFunctionInvalidArgsError, m_pCAudioInfoInstance->initInstance(-1, 1));

        m_pCAudioInfoInstance->initInstance(m_fSampleRate, 1);

        CHECK_EQUAL(kFunctionInvalidArgsError, m_pCAudioInfoInstance->process(0, m_iAudioInfoLength));
        CHECK_EQUAL(kFunctionInvalidArgsError, m_pCAudioInfoInstance->process(&m_pfBuff1, -1));
    }

    TEST_FIXTURE(AudioInfoData, MultiChannel)
    {
        int iNumChannels = 8;
        float ** ppfAudioData = new float *[iNumChannels];
        double dResult = 0;

        for (int c = 0; c < iNumChannels;c++)
        {
            if (c%2)
            {
                ppfAudioData[c] = m_pfBuff1;
            }
            else
            {
                ppfAudioData[c]= m_pfBuff2;
            }
        }

        CSignalGen::generateSine(m_pfBuff1, 1000.F, m_fSampleRate, m_iAudioInfoLength, .5F, 0);
        CSignalGen::generateSine(m_pfBuff2, 1000.F, m_fSampleRate, m_iAudioInfoLength, .25F, 0);

        for (int i = 0; i < m_iAudioInfoLength; i++)
        {
            m_pfBuff2[i]   -= .6F;
        }
        m_pCAudioInfoInstance->initInstance(m_fSampleRate, iNumChannels);

        m_pCAudioInfoInstance->process(ppfAudioData, m_iAudioInfoLength);

        for (int c = 0; c < iNumChannels;c++)
        {
            if (c%2)
            {
                m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kMax, c);
                CHECK_CLOSE(.5, dResult, 1e-3); 
                m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kMin, c);
                CHECK_CLOSE(-.5, dResult, 1e-3); 
                m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kAbsMax, c);
                CHECK_CLOSE(.5, dResult, 1e-3); 
                m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kRange, c);
                CHECK_CLOSE(1.0, dResult, 1e-3); 
                m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kMean, c);
                CHECK_CLOSE(0, dResult, 1e-3); 
                m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kRms, c);
                CHECK_CLOSE(.5/sqrt(2.0), dResult, 1e-3); 
                m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kStd, c);
                CHECK_CLOSE(.5/sqrt(2.0), dResult, 1e-3); 
            }
            else
            {
                m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kMax, c);
                CHECK_CLOSE(-.35, dResult, 1e-3); 
                m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kMin, c);
                CHECK_CLOSE(-.85, dResult, 1e-3); 
                m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kAbsMax, c);
                CHECK_CLOSE(.85, dResult, 1e-3); 
                m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kRange, c);
                CHECK_CLOSE(.5, dResult, 1e-3); 
                m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kMean, c);
                CHECK_CLOSE(-.6, dResult, 1e-3); 
                m_pCAudioInfoInstance->getResult (dResult, CAudioInfo::kStd, c);
                CHECK_CLOSE(.25/sqrt(2.0), dResult, 1e-3); 
            }
        }

        delete [] ppfAudioData;
    }
}

#endif //WITH_TESTS