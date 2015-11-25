#include "MigEditConfig.h"

#ifdef WITH_TESTS
#include <cassert>
#include <cstdio>

#include "UnitTest++.h"

#include "SignalGen.h"

#include "AudioInfo.h"
#include "Matrix.h"

SUITE(Matrix)
{
    struct MatrixData
    {
        MatrixData() 
        {
            //m_ppfBuff   = new float *[m_iNumChannels];
            //m_ppfOut    = new float *[m_iNumChannels];
            //for (int c=0;c<m_iNumChannels;c++)
            //{
            //    m_ppfBuff[c]    = new float [m_iLength];
            //    m_ppfOut[c]     = new float [m_iLength];

            //    CSignalGen::generateSine(m_ppfBuff[c], 1000.F, 44100, m_iLength, .5F/(c+1), (float)(c*M_PI));
            //}
            //CAudioInfo::createInstance(m_phAudioInfo);
            //m_phAudioInfo->initInstance(44100, m_iNumChannels);
            //m_phAudioInfo->process(m_ppfBuff, m_iLength);
            //CMatrix::createInstance(m_phMatrix);
            //m_phMatrix->initInstance(m_phAudioInfo);
        }

        ~MatrixData() 
        {
            //for (int c=0;c<m_iNumChannels;c++)
            //{
            //    delete [] m_ppfOut[c];
            //    delete [] m_ppfBuff[c];
            //}
            //delete [] m_ppfOut;
            //delete [] m_ppfBuff;

            //m_phMatrix->resetInstance();
            //CMatrix::destroyInstance(m_phMatrix);

            //m_phAudioInfo->resetInstance();
            //CAudioInfo::destroyInstance(m_phAudioInfo);
        }

        CMatrix *m_phMatrix;

        //float **m_ppfBuff;
        //float **m_ppfOut;
        //static const int   m_iLength       = 1024;
        //static const int   m_iNumChannels  = 2;
    };

    TEST_FIXTURE(MatrixData, Downmix)
    {
        CHECK_EQUAL(0,0);
    }
}

#endif //WITH_TESTS