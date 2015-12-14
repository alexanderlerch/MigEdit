#include "MigEditConfig.h"

#ifdef WITH_TESTS
#include <cassert>
#include <cstdio>
#include <ctime>

#include "UnitTest++.h"

#include "SignalGen.h"

#include "Util.h"
#include "Ccf.h"

SUITE(Ccf)
{
    struct CcfData
    {
        CcfData() 
        {
            m_pfSrc1    = new float [m_iMaxInpLength];
            m_pfSrc2    = new float [m_iMaxInpLength];
            m_pfRes     = new float [2*m_iMaxInpLength];

            CCcf::create(m_pCCcfInstance);
        }

        ~CcfData() 
        {
            m_pCCcfInstance->reset();
            CCcf::destroy(m_pCCcfInstance);

            delete [] m_pfSrc2;
            delete [] m_pfRes;
            delete [] m_pfSrc1;
        }

        float *m_pfSrc1;
        float *m_pfSrc2;
        float *m_pfRes;

        static const int m_iMaxInpLength  = 1024;

        CCcf *m_pCCcfInstance;
    };

    TEST_FIXTURE(CcfData, SimpleAcf)
    {
        int iLength = m_iMaxInpLength - 2;
        // zerotest
        int aiLengths[2] = { iLength, iLength};
        CUtil::setZero(m_pfSrc1, iLength);
        
        m_pCCcfInstance->init(aiLengths, true);
        CHECK_EQUAL(2* iLength -1, m_pCCcfInstance->getResultLength(aiLengths));

        m_pCCcfInstance->process (m_pfSrc1, m_pfSrc1);
        m_pCCcfInstance->getResult (m_pfRes);
        
        for (int i = 0; i < 2* iLength -1; i++)
            CHECK_CLOSE(0.F,m_pfRes[i], 1e-4);
        
        // impulse test
        m_pCCcfInstance->reset();
        m_pCCcfInstance->init(aiLengths, false);
        m_pfSrc1[m_iMaxInpLength>>1]   = 2.F;
        
        m_pCCcfInstance->process (m_pfSrc1, m_pfSrc1);
        m_pCCcfInstance->getResult (m_pfRes);

        for (int i = 0; i < 2* iLength -1; i++)
        {
            if (i != iLength-1)
                CHECK_CLOSE(0.F,m_pfRes[i], 1e-4);
            else
                CHECK_CLOSE(4.F,m_pfRes[i], 1e-4);
        }

        m_pCCcfInstance->init(aiLengths, true);
        m_pCCcfInstance->process(m_pfSrc1, m_pfSrc1);
        m_pCCcfInstance->getResult(m_pfRes);

        for (int i = 0; i < 2 * iLength - 1; i++)
        {
            if (i != iLength - 1)
                CHECK_CLOSE(0.F, m_pfRes[i], 1e-4);
            else
                CHECK_CLOSE(1.F, m_pfRes[i], 1e-4);
        }

        // rect test
        CSignalGen::generateRect(m_pfSrc1, .5F, static_cast<float>(m_iMaxInpLength), m_iMaxInpLength, .5F);

        m_pCCcfInstance->init(aiLengths, false);

        m_pCCcfInstance->process (m_pfSrc1, m_pfSrc1);
        m_pCCcfInstance->getResult (m_pfRes);

        CHECK_CLOSE(.25F*iLength, m_pfRes[iLength-1], 1e-4F);
        for (int i = 1; i < iLength; i++)
            CHECK_EQUAL(true, m_pfRes[i]>m_pfRes[i-1]);
        for (int i = iLength; i < 2* iLength -2; i++)
            CHECK_EQUAL(true, m_pfRes[i]>m_pfRes[i+1]);

    }
    TEST_FIXTURE(CcfData, Identity)
    {
        int iLength = 1;
        // zerotest
        int aiLengths[2] = { m_iMaxInpLength,iLength  };
        CUtil::setZero(m_pfSrc1, m_iMaxInpLength);
        CUtil::setZero(m_pfSrc2, m_iMaxInpLength);
        m_pfSrc2[0] = 1.F;
        CSignalGen::generateSine(m_pfSrc1, 4.25F, static_cast<float>(m_iMaxInpLength), m_iMaxInpLength);

        m_pCCcfInstance->init(aiLengths, true);
        CHECK_EQUAL(m_iMaxInpLength, m_pCCcfInstance->getResultLength(aiLengths));

        m_pCCcfInstance->process(m_pfSrc1, m_pfSrc2);
        m_pCCcfInstance->getResult(m_pfRes);

        CHECK_ARRAY_CLOSE(m_pfSrc1, m_pfRes, m_iMaxInpLength, 1e-4F);

        aiLengths[0] = iLength;
        aiLengths[1] = m_iMaxInpLength;
    
        CUtil::setZero(m_pfSrc1, m_iMaxInpLength);
        CUtil::setZero(m_pfSrc2, m_iMaxInpLength);
        m_pfSrc1[0] = 1.F;
        CSignalGen::generateSine(m_pfSrc2, 4.25F, static_cast<float>(m_iMaxInpLength), m_iMaxInpLength);
        
        m_pCCcfInstance->init(aiLengths, true);
        CHECK_EQUAL(m_iMaxInpLength, m_pCCcfInstance->getResultLength(aiLengths));

        m_pCCcfInstance->process(m_pfSrc1, m_pfSrc2);
        m_pCCcfInstance->getResult(m_pfRes);

        for (int i = 0; i < m_iMaxInpLength; i++)
        {
            CHECK_CLOSE(m_pfSrc2[i], m_pfRes[m_iMaxInpLength - 1 - i], 1e-4F);
        }

    }
    TEST_FIXTURE(CcfData, Performance)
    {
        int aiLengths[2] = { 30000,30000 };
        float *pfBuffer = new float[aiLengths[0]];
        clock_t	clStartTime = 0,
            clTotalTime[2] = { 0,0 };

        CSignalGen::generateSine(pfBuffer, 117.3F, static_cast<float>(aiLengths[0]), aiLengths[0]);
        m_pCCcfInstance->init(aiLengths, true);
        clStartTime = clock();
        m_pCCcfInstance->process(pfBuffer, pfBuffer);
        clTotalTime[0] = clock() - clStartTime;

        clStartTime = clock();
        for (int i = 0; i< aiLengths[0]; i++)
        {
            float temp = 0.F;
            for (int j = 0; j < aiLengths[0] - i - 1; j++)
            {
                temp += pfBuffer[i] * pfBuffer[i + j];
            }
        }
        clTotalTime[1] = clock() - clStartTime;

        CHECK  (clTotalTime[0] < clTotalTime[1]);

        delete pfBuffer;
    }

}

#endif //WITH_TESTS