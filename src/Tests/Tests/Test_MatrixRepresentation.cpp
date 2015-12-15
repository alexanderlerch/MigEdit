#include "MigEditConfig.h"

#ifdef WITH_TESTS
#include <cassert>
#include <cstdio>

#include "UnitTest++.h"

#include "SignalGen.h"

#include "AudioInfo.h"
#include "MatrixRepresentation.h"

SUITE(MatrixRepresentation)
{
    struct MatrixRepresentationData
    {
        MatrixRepresentationData() 
        {
            m_fSampleRate   = 44100;
            m_ppfBuff   = new float *[m_iNumChannels];
            for (int c=0;c<m_iNumChannels;c++)
            {
                m_ppfBuff[c]    = new float [m_iLength];
            }

            CMatrixRepresentation::create(m_pCMatrixRepresentation);
            m_pCMatrixRepresentation->init(m_iLength, m_fSampleRate, m_iNumChannels,m_iBlockLength, m_iHopLength, m_Result);
        }

        ~MatrixRepresentationData() 
        {
            for (int c=0;c<m_iNumChannels;c++)
            {
                delete [] m_ppfBuff[c];
            }
            delete [] m_ppfBuff;

            m_pCMatrixRepresentation->reset();
            
        }

        CMatrixRepresentation *m_pCMatrixRepresentation;

        static const int m_iLength      = 161717;
        static const int m_iNumChannels = 4;
        float m_fSampleRate;
        static const int m_iBlockLength = 1024;
        static const int m_iHopLength   = 256;

        CMatrixRepresentationResult m_Result;

        float **m_ppfBuff;
    };

    TEST_FIXTURE(MatrixRepresentationData, Process)
    {
        for (int c=0; c<m_iNumChannels; c++)
            CSignalGen::generateSine(m_ppfBuff[c], 1000.F + 1000*c, m_fSampleRate, m_iLength, .5F, static_cast<float>(M_PI*.5 * c));

        m_pCMatrixRepresentation->process(m_ppfBuff, m_iLength);

         // check dimensions
        CHECK_EQUAL(513, m_Result.getNumCols());
        CHECK_EQUAL(635, m_Result.getNumRows());
        CHECK_EQUAL(m_iNumChannels, m_Result.getNumChannels());

        for (int c=0; c<m_iNumChannels; c++)
        {
            float fMax;
            long long int iMax;
            CMatrix *pMatrixResult = m_Result.getResultPtr (c);

            pMatrixResult->getRow(17, m_ppfBuff[c], pMatrixResult->getNumCols());
            CUtil::findMax(m_ppfBuff[c],fMax, iMax, pMatrixResult->getNumCols());

            CHECK_CLOSE(1000.F + 1000*c, iMax*m_fSampleRate/m_iBlockLength,m_fSampleRate/m_iBlockLength);
        }

    }
    TEST_FIXTURE(MatrixRepresentationData, Api)
    {
        //TODO: implement me
    }
}

#endif //WITH_TESTS