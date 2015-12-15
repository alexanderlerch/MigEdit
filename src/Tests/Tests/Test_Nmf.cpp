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
        }

        ~NmfData() 
        {
            m_pCNmf->reset();
            CNmfIf::destroy(m_pCNmf);
            m_CInit.reset();
        }

        CNmfIf *m_pCNmf;
        CNmfParametrization m_CInit;
        CNmfResult m_CResult;
    };

    TEST_FIXTURE(NmfData, BasicNmf)
    {
        CMatrix Result[kNumMatrices][kNumSplits];
        int iNmfRank = 3;
        CMatrix myMatrix = CMatrix(8, 13);

        for (int j = 0; j < 4; j++)
            myMatrix.setElement(2, j, 1.F);
        for (int j = 4; j < 9; j++)
            myMatrix.setElement(5, j, 1.F);
        for (int j = 9; j < 13; j++)
            myMatrix.setElement(7, j, 1.F);

        m_CInit.init(8,iNmfRank);

        m_pCNmf->init(m_CInit);

        m_pCNmf->process (&myMatrix, m_CResult);
        //for (int i = 0; i < kNumMatrices; i++)
        //{
        //    for (int j = 0; j < kNumSplits; j++)
        //    {
        //        Result[i][j] = m_CResult.getMatrix((Matrices_t)i,(MatrixSplit_t)j);
        //        Result[i][j].dbgPrint2StdOut();
        //    }
        //}
        
        CHECK_CLOSE(3,Result[kDict][kSplit1].getSum(), 1e-4F);
        for (int j = 0; j < 3; j++)
            CHECK_CLOSE(1,Result[kDict][kSplit1].getColNorm(j,1),1e-4F);
        CHECK_CLOSE(13,Result[kAct][kSplit1].getSum(), 1e-4F);
        for (int j = 0; j < 13; j++)
            CHECK_CLOSE(1,Result[kAct][kSplit1].getColNorm(j,1),1e-4F);

    }
}

#endif //WITH_TESTS