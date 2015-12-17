#if !defined(__SignalGen_hdr__)
#define __SignalGen_hdr__

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdlib>

#include "ErrorDef.h"
#include "Util.h"
class CSignalGen
{
public:
    static Error_t generateSine (float *pfOutBuf, float fFreqInHz, float fSampleFreqInHz, int iLength, float fAmplitude = 1.F, float fStartPhaseInRad = 0.F)
    {
        if (!pfOutBuf)
            return kFunctionInvalidArgsError;

        for (int i = 0; i < iLength; i++)
        {
            pfOutBuf[i] = fAmplitude * static_cast<float>(sin (2*M_PI*fFreqInHz * i/fSampleFreqInHz + fStartPhaseInRad));
        }

        return kNoError;
    }
    static Error_t generateRect (float *pfOutBuf, float fFreqInHz, float fSampleFreqInHz, int iLength, float fAmplitude = 1.F)
    {
        if (!pfOutBuf)
            return kFunctionInvalidArgsError;

        float fPeriodLength = fSampleFreqInHz / fFreqInHz;
        for (int i = 0; i < iLength; i++)
        {
            if (i%CUtil::float2int<int>(fPeriodLength) <= .5*fPeriodLength)
            {            
                pfOutBuf[i] = fAmplitude;
            }
            else
            {
                pfOutBuf[i] = -fAmplitude;
            }
        }

        return kNoError;
    }
    static Error_t generateSaw (float *pfOutBuf, float fFreqInHz, float fSampleFreqInHz, int iLength, float fAmplitude = 1.F)
    {
        if (!pfOutBuf)
            return kFunctionInvalidArgsError;

        float fIncr = 2*fAmplitude / fSampleFreqInHz * fFreqInHz;
        pfOutBuf[0] = 0;
        for (int i = 1; i < iLength; i++)
        {
            pfOutBuf[i] = fmodf(pfOutBuf[i-1] + fIncr + fAmplitude, 2*fAmplitude) - fAmplitude;
        }

        return kNoError;
    }
    static Error_t generateDc (float *pfOutBuf, int iLength, float fAmplitude = 1.F)
    {
        if (!pfOutBuf)
            return kFunctionInvalidArgsError;

        for (int i = 0; i < iLength; i++)
        {
            pfOutBuf[i] = fAmplitude;
        }

        return kNoError;
    }
    static Error_t generateNoise (float *pfOutBuf, int iLength, float fAmplitude = 1.F)
    {
        if (!pfOutBuf)
            return kFunctionInvalidArgsError;

        for (int i = 0; i < iLength; i++)
        {
            pfOutBuf[i] = rand()*fAmplitude/RAND_MAX;
        }

        return kNoError;
    }
};
#endif // __SignalGen_hdr__