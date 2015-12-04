#include "MigEditConfig.h"

#ifdef WITH_TESTS
#include <cassert>
#include <cstdio>

#include "UnitTest++.h"

#include "SignalGen.h"

#include "AudioInfo.h"
#include "Matrix.h"
#include "NmfIf.h"

SUITE(Nmf)
{
    struct NmfData
    {
        NmfData() 
        {
            CNmfIf::create(m_pCNmf);
            m_pCInitAndResult   = new CNmfSharedData(8,3);

            m_pCNmf->init(*m_pCInitAndResult);
        }

        ~NmfData() 
        {
            m_pCNmf->reset();
            CNmfIf::destroy(m_pCNmf);
            m_pCInitAndResult->reset();
            delete m_pCInitAndResult;
        }

        CNmfIf *m_pCNmf;
        CNmfSharedData *m_pCInitAndResult;
    };

    TEST_FIXTURE(NmfData, Basic)
    {
        CMatrix myMatrix = CMatrix(8, 13);
        CMatrix myFixed  = CMatrix(8,3);
        CMatrix myFixedAct = CMatrix(3,13);
        myFixed.setOnes();
        myFixedAct.setOnes();
        CMatrix myNonFixed = myFixed;
        CMatrix myNonFixedAct = myFixedAct;

        for (int j = 0; j < 4; j++)
            myMatrix.setElement(2, j, 1.F);
        for (int j = 4; j < 9; j++)
            myMatrix.setElement(5, j, 1.F);
        for (int j = 9; j < 13; j++)
            myMatrix.setElement(7, j, 1.F);

        m_pCInitAndResult->setMatrix(kDict,kSplit1, &myNonFixed);
        m_pCInitAndResult->setMatrix(kAct,kSplit1, &myNonFixedAct);
        m_pCInitAndResult->setMatrix(kDict,kSplit2, &myFixed);
        m_pCInitAndResult->setMatrix(kAct,kSplit2, &myFixedAct);
        m_pCInitAndResult->finalizeSettings();
        m_pCNmf->process (&myMatrix);
    }
}

#endif //WITH_TESTS