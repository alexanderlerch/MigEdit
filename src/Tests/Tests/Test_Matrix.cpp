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
            m_ppfBuff   = new float *[m_iNumRows];
            for (int c=0;c<m_iNumRows;c++)
            {
                m_ppfBuff[c]    = new float [m_iNumCols];
            }

            Matrix.init(m_iNumRows, m_iNumCols);
        }

        ~MatrixData() 
        {
            for (int c=0;c<m_iNumRows;c++)
            {
                delete [] m_ppfBuff[c];
            }
            delete [] m_ppfBuff;

            Matrix.reset();
            
        }

        CMatrix Matrix;

        static const int m_iNumRows = 3;
        static const int m_iNumCols = 4;

        float **m_ppfBuff;
    };

    TEST_FIXTURE(MatrixData, SetGet)
    {
        float afArray[4];

        CHECK_EQUAL(m_iNumRows, Matrix.getNumRows());
        CHECK_EQUAL(m_iNumCols, Matrix.getNumCols());
        for (int i = 0; i < m_iNumRows; i++)
        {
            for (int j = 0; j < m_iNumCols; j++)
            {
                m_ppfBuff[i][j] = static_cast<float>(i*j);
                Matrix.setElement(i,j, static_cast<float>(i*j));
            }

        }
        
        Matrix.getCol(2, afArray, m_iNumRows);
        for (int i = 0; i < m_iNumRows; i++)
            CHECK_EQUAL(m_ppfBuff[i][2],afArray[i]);

        for (int i = 0; i < m_iNumRows; i++)
        {
            for (int j = 0; j < m_iNumCols; j++)
            {
                CHECK_EQUAL(m_ppfBuff[i][j],Matrix.getElement(i,j));   
            }
            Matrix.getRow(i, afArray, m_iNumCols);
            CHECK_ARRAY_EQUAL(m_ppfBuff[i],afArray, m_iNumCols);
        }

        Matrix.setZero();
        for (int i = 0; i < m_iNumRows; i++)
        {
            for (int j = 0; j < m_iNumCols; j++)
            {
                CHECK_EQUAL(0,Matrix.getElement(i,j));   
            }
        }
    }
    TEST_FIXTURE(MatrixData, Api)
    {
        float afArray[4];
        CMatrix TestMatrix;

        //no init
        TestMatrix.reset();

        CHECK_EQUAL(0, TestMatrix.getNumRows());
        CHECK_EQUAL(0, TestMatrix.getNumCols());

        CHECK_EQUAL(CMatrix::kMatrixIllegalFunctionParam, TestMatrix.getRow(2,afArray, m_iNumCols));


        for (int i = 0; i < m_iNumRows; i++)
        {
            for (int j = 0; j < m_iNumCols; j++)
            {
                m_ppfBuff[i][j] = static_cast<float>(i*j);
                Matrix.setElement(i,j, static_cast<float>(i*j));
            }
        }
        CHECK_EQUAL(CMatrix::kMatrixIllegalFunctionParam, Matrix.setElement(m_iNumRows,1, static_cast<float>(m_iNumRows)));
    }
    TEST_FIXTURE(MatrixData, Operators)
    {
        CMatrix TestMatrix;
        CMatrix Result;
        float afArray[4] = {5,4,3,1};
        Result.init(1,1);
        TestMatrix.init(4,1);

        TestMatrix.setCol (0,afArray,4);
        //[ 0 0 0 0
        //  0 1 2 3
        //  0 2 4 6]
        for (int i = 0; i < m_iNumRows; i++)
        {
            for (int j = 0; j < m_iNumCols; j++)
            {
                m_ppfBuff[i][j] = static_cast<float>(i*j);
                Matrix.setElement(i,j, static_cast<float>(i*j));
            }
        }

        //Result = TestMatrix * Matrix;

        //CHECK_EQUAL(0, TestMatrix.getNumRows());
        //CHECK_EQUAL(0, TestMatrix.getNumCols());

        Result =  Matrix * TestMatrix;
        CHECK_EQUAL(m_iNumRows, Result.getNumRows());
        CHECK_EQUAL(1, Result.getNumCols());
        CHECK_EQUAL(0, Result.getElement(0,0));
        CHECK_EQUAL(13, Result.getElement(1,0));
        CHECK_EQUAL(26, Result.getElement(2,0));

        TestMatrix = Matrix;
        for (int i = 0; i < m_iNumRows; i++)
        {
            TestMatrix.getRow(i, afArray, m_iNumCols);
            CHECK_ARRAY_EQUAL(m_ppfBuff[i],afArray, m_iNumCols);
        }
        TestMatrix = Matrix + Matrix;
        for (int i = 0; i < m_iNumRows; i++)
        {
            for (int j = 0; j < m_iNumCols; j++)
            {
                CHECK_EQUAL(m_ppfBuff[i][j]*2,TestMatrix.getElement(i,j));   
            }
        }
        TestMatrix = Matrix * 2;
        for (int i = 0; i < m_iNumRows; i++)
        {
            for (int j = 0; j < m_iNumCols; j++)
            {
                CHECK_EQUAL(m_ppfBuff[i][j]*2,TestMatrix.getElement(i,j));   
            }
        }
        TestMatrix = Matrix + 2.5F;
        for (int i = 0; i < m_iNumRows; i++)
        {
            for (int j = 0; j < m_iNumCols; j++)
            {
                CHECK_EQUAL(m_ppfBuff[i][j]+2.5F,TestMatrix.getElement(i,j));   
            }
        }
        TestMatrix = TestMatrix - 2.5F;
        for (int i = 0; i < m_iNumRows; i++)
        {
            for (int j = 0; j < m_iNumCols; j++)
            {
                CHECK_EQUAL(m_ppfBuff[i][j],TestMatrix.getElement(i,j));   
            }
        }

        TestMatrix = TestMatrix;
    }
    TEST_FIXTURE(MatrixData, Transpose)
    {
        for (int i = 0; i < m_iNumRows; i++)
        {
            for (int j = 0; j < m_iNumCols; j++)
            {
                m_ppfBuff[i][j] = static_cast<float>(i*j);
                Matrix.setElement(i,j, static_cast<float>(m_ppfBuff[i][j]));
            }
        }
        Matrix.transpose();

        for (int i = 0; i < m_iNumRows; i++)
        {
            for (int j = 0; j < m_iNumCols; j++)
            {
                CHECK_EQUAL(m_ppfBuff[i][j],Matrix.getElement(j,i));
            }
        }
    }
    TEST_FIXTURE(MatrixData, Norm)
    {
        for (int i = 0; i < m_iNumRows; i++)
        {
            for (int j = 0; j < m_iNumCols; j++)
            {
                m_ppfBuff[i][j] = static_cast<float>(i*j);
                Matrix.setElement(i,j, static_cast<float>(m_ppfBuff[i][j]));
            }
        }
        CHECK_CLOSE(9, Matrix.getNorm(1), 1e-4);
        CHECK_CLOSE(8.36660, Matrix.getNorm(2), 1e-4);
        Matrix.transpose();
        CHECK_CLOSE(9, Matrix.getNorm(1), 1e-4);
        CHECK_CLOSE(8.36660, Matrix.getNorm(2), 1e-4);
    }
}

#endif //WITH_TESTS