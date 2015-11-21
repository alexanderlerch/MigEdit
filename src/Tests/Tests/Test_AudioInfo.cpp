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
}

#endif //WITH_TESTS