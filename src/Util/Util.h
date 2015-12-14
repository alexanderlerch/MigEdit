#if !defined(__Util_hdr__)
#define __Util_hdr__

#include <cassert>
#include <cstring>
#include <limits>

class CUtil
{
public:
    template<typename T>
    static T float2int (float fInput)
    {
        if (fInput >= 0.F)
            return static_cast<T>(fInput + .5F);
        else
            return static_cast<T>(fInput - .5F);
    }
    template<typename T>
    static T double2int (double fInput)
    {
        if (fInput >= 0)
            return static_cast<T>(fInput + .5);
        else
            return static_cast<T>(fInput - .5);
    }

    static bool isPowOf2 (int n) 
    {
        return !(n & (n-1));
    }

    static int nextPowOf2(int n)
    {
        int iOrder = 0;

        if (n == 0)
            return 0;

        while (n>>iOrder)
            iOrder++;

        if (!(n%(1<<(iOrder-1))))
            iOrder--;

        return (1<<(iOrder));
    }

    template<typename T>
    static void setZero (T *pInput, int iLength)
    {
        assert (iLength >= 0);
        assert (pInput);

        if (iLength > 0)
            memset (pInput, 0, sizeof(T)*iLength);
    }
    template<typename T>
    static void copyBuff(T *pDest, const T *pSource, int iLength)
    {
        assert(iLength >= 0);

        if (iLength > 0)
        {
            assert(pDest);
            assert(pSource);
            memcpy(pDest, pSource, sizeof(T)*iLength);
        }
    }
    template<typename T>
    static void flipBuff(T *pSrcDest, int iLength)
    {
        assert(iLength >= 0);

        if (iLength > 0)
        {
            assert(pSrcDest);

            int iLoopLength = iLength / 2; // integer division!
            for (int i = 0; i < iLoopLength; i++)
            {
                T Tmp                       = pSrcDest[i];
                pSrcDest[i]                 = pSrcDest[iLength - 1 - i];
                pSrcDest[iLength - 1 - i]   = Tmp;
            }
        }
    }
    template<typename T>
    static void moveBuff (T *pSrcDest, int iDestIdx, int iSrcIdx, int iLength)
    {
        assert (iLength >= 0);
        assert (pSrcDest);

        if (iLength > 0)
            memmove (&pSrcDest[iDestIdx], &pSrcDest[iSrcIdx], sizeof(T)*iLength);
    }

    static void setValue (float *pDest, float fValue, int iLength)
    {
        assert (iLength >= 0);
        assert (pDest);

        for (int i = 0; i < iLength; i++)
            pDest[i] = fValue;
    }
    static void mulBuffC (float *pSrcDest, float fScale, int iLength)
    {
        assert (iLength >= 0);
        assert (pSrcDest);
        
        for (int i = 0; i < iLength; i++)
            pSrcDest[i] *= fScale;
    }

    static void mulBuff (float *pSrcDest, const float *pSrc, int iLength)
    {
        assert (iLength >= 0);
        assert (pSrcDest);
        assert (pSrc);

        for (int i = 0; i < iLength; i++)
            pSrcDest[i] *= pSrc[i];
    }

    static float mulBuffScalar (const float *pSrc1, const float *pSrc2, int iLength)
    {
        assert (iLength >= 0);
        assert (pSrc1);
        assert (pSrc2);
        float  fResult = 0;

        for (int i = 0; i < iLength; i++)
            fResult += pSrc1[i] * pSrc2[i];

        return fResult;
    }

    static void divBuff (float *pSrcDest, const float *pSrc, int iLength)
    {
        assert (iLength >= 0);
        assert (pSrcDest);
        assert (pSrc);

        for (int i = 0; i < iLength; i++)
        {
            assert(pSrc[i] != 0);
            pSrcDest[i] /= pSrc[i];
        }
    }

    static void addBuff (float *pSrcDest, const float *pSrc, int iLength)
    {
        assert (iLength >= 0);
        assert (pSrcDest);
        assert (pSrc);

        for (int i = 0; i < iLength; i++)
            pSrcDest[i] += pSrc[i];
    }

    static void subBuff (float *pSrcDest, const float *pSrc, int iLength)
    {
        assert (iLength >= 0);
        assert (pSrcDest);
        assert (pSrc);

        for (int i = 0; i < iLength; i++)
            pSrcDest[i] -= pSrc[i];
    }

    static float sum (float *pSrc, int iLength, bool bAbs = false)
    {
        assert (iLength >= 0);
        assert (pSrc);

        float fResult = 0;
        if (bAbs)
        {
            for (int i = 0; i < iLength; i++)
                fResult += std::abs(pSrc[i]);
        }
        else
        {
            for (int i = 0; i < iLength; i++)
                fResult += pSrc[i];
        }
        return fResult;
    }

    static bool isBuffEqual (float *pSrc1, const float *pSrc2, int iLength)
    {
        assert (iLength >= 0);
        assert (pSrc1);
        assert (pSrc2);

        return (memcmp (pSrc1, pSrc2, iLength * sizeof(float)) == 0);
    }

    template<typename T>
    static double getMax (T *pSrc, long long int iLength, bool bAbs = false)
    {
        T fMax;
        long long iMax;

        findMax<T>(pSrc, fMax, iMax, iLength, bAbs);

        return fMax;
    }
    template<typename T>
    static double getMin (T *pSrc, long long int iLength, bool bAbs = false)
    {
        T fMin;
        long long iMin;

        findMin<T>(pSrc, fMin, iMin, iLength, bAbs);

        return fMin;
    }
    template<typename T>
    static double getMean (T *pSrc, long long int iLength)
    {
        assert (iLength >= 0);

        double dMean = 0;

        for (int i=0; i < iLength; i++)
        {
            dMean  += pSrc[i];
        }

        if (iLength > 0)
        {
            dMean  /= iLength;
        }

        return dMean;
    }
    template<typename T>
    static double getStd (T *pSrc, long long int iLength, double dMean = std::numeric_limits<double>::max())
    {
        assert (iLength >= 0);

        double dStd = 0;

        if (dMean == std::numeric_limits<double>::max())
        {
            dMean   = getMean(pSrc, iLength);
        }

        for (int i=0; i < iLength; i++)
        {
            dStd   += (pSrc[i] - dMean) * (pSrc[i] - dMean);
        }

        if (iLength > 1)
        {
            dStd   /= (iLength - 1);
        }

        return std::sqrt(dStd);
    }
    template<typename T>
    static double getRms (T *pSrc, long long int iLength)
    {
        assert (iLength >= 0);

        double dRms = 0;


        for (int i=0; i < iLength; i++)
        {
            dRms   += pSrc[i] * pSrc[i];
        }

        if (iLength > 0)
        {
            dRms   /= iLength;
        }

        return std::sqrt(dRms);
    }
    template<typename T>
    static void findMax (T *pSrc, T &fMax, long long &iMax, long long int iLength, bool bAbs = false)
    {
        assert (iLength >= 0);
        assert (pSrc);

        fMax    = -std::numeric_limits<T>::max();
        iMax    = -1;

        for (int i = 0; i < iLength; i++)
        {
            T fCurr   = (bAbs)? std::abs(pSrc[i]) : pSrc[i];

            if (fCurr > fMax)
            {
                fMax = fCurr;
                iMax = i;
            }
        }
    }
    template<typename T>
    static void findMin (T *pSrc, T &fMin, long long &iMin, long long int iLength, bool bAbs = false)
    {
        assert (iLength >= 0);
        assert (pSrc);

        fMin    = std::numeric_limits<T>::max();
        iMin    = -1;

        for (int i = 0; i < iLength; i++)
        {
            T fCurr   = (bAbs)? std::abs(pSrc[i]) : pSrc[i];

            if (fCurr < fMin)
            {
                fMin    = fCurr;
                iMin    = i;
            }
        }
    }
};
#endif // __Util_hdr__