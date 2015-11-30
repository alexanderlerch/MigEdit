#include "MigEditConfig.h"

#ifdef WITH_TESTS
#include <cassert>
#include <cstdio>

#include "UnitTest++.h"

#include "SignalGen.h"

#include "AudioInfo.h"
#include "Preproc.h"

SUITE(Preproc)
{
    struct PreprocData
    {
        PreprocData() 
        {
            m_ppfBuff   = new float *[m_iNumChannels];
            m_ppfOut    = new float *[m_iNumChannels];
            for (int c=0;c<m_iNumChannels;c++)
            {
                m_ppfBuff[c]    = new float [m_iLength];
                m_ppfOut[c]     = new float [m_iLength];

                CSignalGen::generateSine(m_ppfBuff[c], 1000.F, 44100, m_iLength, .5F/(c+1), (float)(c*M_PI));
            }
            CAudioInfo::create(m_phAudioInfo);
            m_phAudioInfo->init(44100, m_iNumChannels);
            m_phAudioInfo->process(m_ppfBuff, m_iLength);
            CPreproc::create(m_phPreproc);
            m_phPreproc->init(m_phAudioInfo);
        }

        ~PreprocData() 
        {
            for (int c=0;c<m_iNumChannels;c++)
            {
                delete [] m_ppfOut[c];
                delete [] m_ppfBuff[c];
            }
            delete [] m_ppfOut;
            delete [] m_ppfBuff;

            m_phPreproc->reset();
            CPreproc::destroy(m_phPreproc);

            m_phAudioInfo->reset();
            CAudioInfo::destroy(m_phAudioInfo);
        }

        CPreproc *m_phPreproc;
        CAudioInfo *m_phAudioInfo;

        float **m_ppfBuff;
        float **m_ppfOut;
        static const int   m_iLength       = 1024;
        static const int   m_iNumChannels  = 2;
        

    };

    TEST_FIXTURE(PreprocData, Downmix)
    {
        m_phPreproc->setStepActive(CPreproc::kPpDownmix, true);

        m_phPreproc->process (m_ppfBuff, &m_ppfOut[0], m_iLength);

        for (int i = 0; i < m_iLength; i++)
            CHECK_CLOSE(m_ppfOut[0][i],  0.5F*(m_ppfBuff[0][i] + m_ppfBuff[1][i]), 1e-4F);
    }

    TEST_FIXTURE(PreprocData, Normalize)
    {
        m_phPreproc->setStepActive(CPreproc::kPpNormalize, true);
        m_phPreproc->setStepActive(CPreproc::kPpDownmix, true);

        m_phPreproc->process (m_ppfBuff, &m_ppfOut[0], m_iLength);

        CHECK_CLOSE(1.F,  CUtil::getMax(m_ppfOut[0], m_iLength, true), 1e-4F);
        
        m_phPreproc->setStepActive(CPreproc::kPpNormalize, true);
        m_phPreproc->setStepActive(CPreproc::kPpDownmix, false);

        m_phPreproc->process (m_ppfBuff, m_ppfOut, m_iLength);
        m_phAudioInfo->process(m_ppfOut, m_iLength);
        double dTmp;
        m_phAudioInfo->getResult(dTmp,CAudioInfo::kAbsMax,0);
        CHECK_CLOSE(1.F,  dTmp, 1e-4F);
        m_phAudioInfo->getResult(dTmp,CAudioInfo::kAbsMax,1);
        CHECK_CLOSE(.5F,  dTmp, 1e-4F);
    }
}

#endif //WITH_TESTS