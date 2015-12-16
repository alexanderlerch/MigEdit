#include <string>
#include <cassert>

#include "Util.h"
#include "Fft.h"
#include "MatrixRepresentation.h"
#include "InputBuffSrc.h"



CMatrixRepresentationResult::CMatrixRepresentationResult() :
    m_bIsInitialized(false),
    m_ppCResult(0),
    m_iNumChannels(0)
{
    CUtil::setZero(m_apfTickLabels,kNumAxes);

    reset();
}

CMatrixRepresentationResult::~CMatrixRepresentationResult()
{
    reset();
}

Error_t CMatrixRepresentationResult::init( int iNumRows, int iNumCols, int iNumChannels /*= 1*/ )
{
    if (iNumRows <= 0 || iNumCols <= 0 || iNumChannels <=0)
        return kFunctionInvalidArgsError;

    reset();

    m_ppCResult = new CMatrix*[iNumChannels];
    for (int c = 0; c< iNumChannels; c++)
    {
        m_ppCResult[c]      = new CMatrix();
        if (m_ppCResult[c]->init(iNumRows,iNumCols) != CMatrix::kMatrixNoError)
            return kMemError;
    }

    m_apfTickLabels[kRow]   = new float [iNumRows];        
    m_apfTickLabels[kCol]   = new float [iNumCols];

    m_iNumChannels          = iNumChannels;

    m_bIsInitialized        = true;

    return kNoError;
}

Error_t CMatrixRepresentationResult::reset()
{
    m_bIsInitialized    = false;

    for (int c = 0; c< m_iNumChannels; c++)
    {
        delete m_ppCResult[c];
        m_ppCResult[c]  = 0;
    }
    delete [] m_ppCResult;
    m_ppCResult = 0;

    for (int i = 0; i < kNumAxes; i++)
    {
        delete [] m_apfTickLabels[i];
        m_apfTickLabels[i]  = 0;

        m_asLabels[i].clear();
        m_asUnits[i].clear();
    }

    return kNoError;
}

Error_t CMatrixRepresentationResult::setAxisTickLabels( enum Axes_t eAxis, const float *pfAxisTickLabels, int iNumOfAxisTickLabels )
{
    if (!m_bIsInitialized)
        return kNotInitializedError;
    int  iLength = (eAxis == kRow)? m_ppCResult[0]->getNumRows() : m_ppCResult[0]->getNumCols();

    if (!pfAxisTickLabels || iNumOfAxisTickLabels != iLength)
        return kFunctionInvalidArgsError;

    assert(m_apfTickLabels[eAxis] != 0);

    CUtil::copyBuff(m_apfTickLabels[eAxis], pfAxisTickLabels, iNumOfAxisTickLabels);

    return kNoError;
}

Error_t CMatrixRepresentationResult::getAxisTickLabels( enum Axes_t eAxis, float *pfAxisTickLabels, int iNumOfAxisTickLabels ) const
{
    if (!m_bIsInitialized)
        return kNotInitializedError;
    int  iLength = (eAxis == kRow)? m_ppCResult[0]->getNumRows() : m_ppCResult[0]->getNumCols();

    if (!pfAxisTickLabels || iNumOfAxisTickLabels != iLength)
        return kFunctionInvalidArgsError;

    assert(m_apfTickLabels[eAxis] != 0);

    CUtil::copyBuff(pfAxisTickLabels, m_apfTickLabels[eAxis], iNumOfAxisTickLabels);

    return kNoError;
}

Error_t CMatrixRepresentationResult::setAxisLabel( enum Axes_t eAxis, std::string sAxisLabel )
{
    if (!m_bIsInitialized)
        return kNotInitializedError;

    m_asLabels[eAxis]   = sAxisLabel;

    return kNoError;
}

std::string CMatrixRepresentationResult::getAxisLabel( enum Axes_t eAxis ) const
{
    return m_asLabels[eAxis];
}

Error_t CMatrixRepresentationResult::setAxisUnit( enum Axes_t eAxis, std::string sAxisLabel )
{
    if (!m_bIsInitialized)
        return kNotInitializedError;

    m_asUnits[eAxis]   = sAxisLabel;

    return kNoError;
}

std::string CMatrixRepresentationResult::getAxisUnit( enum Axes_t eAxis ) const
{
    return m_asUnits[eAxis];
}

Error_t CMatrixRepresentationResult::setResultCol( int iCol, const float *pfValues, int iNumValues, int iChannelIdx /*= 0*/ )
{
    CMatrix::MatrixError_t rErr;
    if (!m_bIsInitialized)
        return kNotInitializedError;

    if (iChannelIdx < 0 || iChannelIdx >= m_iNumChannels)
        return kFunctionInvalidArgsError;

    rErr = m_ppCResult[iChannelIdx]->setCol(iCol, pfValues, iNumValues);

    if (rErr != CMatrix::kMatrixNoError)
        return kUnknownError;
    else
        return kNoError;

}

Error_t CMatrixRepresentationResult::setResultRow( int iRow, const float *pfValues, int iNumValues, int iChannelIdx /*= 0*/ )
{
    CMatrix::MatrixError_t rErr;
    if (!m_bIsInitialized)
        return kNotInitializedError;

    if (iChannelIdx < 0 || iChannelIdx >= m_iNumChannels)
        return kFunctionInvalidArgsError;

    rErr = m_ppCResult[iChannelIdx]->setRow(iRow, pfValues, iNumValues);

    if (rErr != CMatrix::kMatrixNoError)
        return kUnknownError;
    else
        return kNoError;
}

Error_t CMatrixRepresentationResult::getCol( int iCol, float *pfValues, int iNumValues, int iChannelIdx /*= 0*/ ) const
{
    CMatrix::MatrixError_t rErr;
    if (!m_bIsInitialized)
        return kNotInitializedError;

    if (iChannelIdx < 0 || iChannelIdx >= m_iNumChannels)
        return kFunctionInvalidArgsError;

    rErr = m_ppCResult[iChannelIdx]->getCol(iCol, pfValues, iNumValues);

    if (rErr != CMatrix::kMatrixNoError)
        return kUnknownError;
    else
        return kNoError;

}

Error_t CMatrixRepresentationResult::getRow( int iRow, float *pfValues, int iNumValues, int iChannelIdx /*= 0*/ ) const
{
    CMatrix::MatrixError_t rErr;
    if (!m_bIsInitialized)
        return kNotInitializedError;

    if (iChannelIdx < 0 || iChannelIdx >= m_iNumChannels)
        return kFunctionInvalidArgsError;

    rErr = m_ppCResult[iChannelIdx]->getRow(iRow, pfValues, iNumValues);

    if (rErr != CMatrix::kMatrixNoError)
        return kUnknownError;
    else
        return kNoError;

}

CMatrix* CMatrixRepresentationResult::getResultPtr( int iChannelIdx /*= 0*/ )
{
    if (iChannelIdx < 0 || iChannelIdx >= m_iNumChannels || !m_bIsInitialized)
        return 0;

    return m_ppCResult[iChannelIdx];
}

int CMatrixRepresentationResult::getNumRows() const
{
    if (!m_bIsInitialized)
        return 0;

    assert(m_ppCResult != 0 && m_ppCResult[0] != 0);

    return m_ppCResult[0]->getNumRows();
}

int CMatrixRepresentationResult::getNumCols() const
{
    if (!m_bIsInitialized)
        return 0;

    assert(m_ppCResult != 0 && m_ppCResult[0] != 0);

    return m_ppCResult[0]->getNumCols();
}

int CMatrixRepresentationResult::getNumChannels() const
{
    if (!m_bIsInitialized)
        return 0;

    return m_iNumChannels;
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////


CMatrixRepresentation::CMatrixRepresentation() :
    m_bIsInitialized(false),
    m_pCInputBuffer(0),
    m_pCFft(0),
    m_iNumChannels(0),
    m_ppfProcBuffer(0),
    m_pCResult(0)
{
    reset ();
}

CMatrixRepresentation::~CMatrixRepresentation() 
{
    reset ();
}

Error_t CMatrixRepresentation::create( CMatrixRepresentation*& pCMatrixRepresentation )
{
    pCMatrixRepresentation = new CMatrixRepresentation ();

    if (!pCMatrixRepresentation)
        return kMemError;

    return kNoError;
}

Error_t CMatrixRepresentation::destroy( CMatrixRepresentation*& pCMatrixRepresentation )
{
    if (!pCMatrixRepresentation)
        return kNoError;

    delete pCMatrixRepresentation;
    pCMatrixRepresentation   = 0;

    return kNoError;
}


Error_t CMatrixRepresentation::reset()
{
    for (int c = 0; c < m_iNumChannels; c++)
        delete [] m_ppfProcBuffer[c];
    delete [] m_ppfProcBuffer;
    m_ppfProcBuffer = 0;

    CFft::destroy(m_pCFft);
    delete m_pCInputBuffer;
    m_pCInputBuffer = 0;

    m_pCResult      = 0;

    CUtil::setZero(m_aiParams,kNumParams);
    m_iNumChannels  = 0;

    return kNoError;
}

Error_t CMatrixRepresentation::init( int iLengthInSamples, float fSampleRateInHz, int iNumChannels, int iBlockLength, int iHopLength, CMatrixRepresentationResult &Result )
{
    if (iLengthInSamples <= 0 || fSampleRateInHz <+ 0 || iNumChannels <= 0 || iBlockLength <= 0 || iHopLength <= 0 || !CUtil::isPowOf2(iBlockLength))
        return kFunctionInvalidArgsError;

    reset();

    int iNumObservations    = 0;
    int iNumFreqBins        = iBlockLength/2 + 1;
    float *pfTmpTickLabels  = 0;

    // compute number of blocks
    iNumObservations = calcNumProcessingBlocks(iLengthInSamples, iBlockLength, iHopLength);

    // alloc local memory
    pfTmpTickLabels = new float[std::max(iNumObservations,iNumFreqBins)];
    
    // allocate internal memory and create private instances
    CFft::create(m_pCFft);
    m_pCFft->init(iBlockLength);
    m_pCInputBuffer     = new CInputBuffSrc<float>(iNumChannels,iBlockLength, iBlockLength-iHopLength);
    m_ppfProcBuffer     = new float*[iNumChannels];
    for (int c = 0; c < iNumChannels; c++)
    {
        m_ppfProcBuffer[c]  = new float [iBlockLength];
    }

    // init result class
    Result.init(iNumObservations, iNumFreqBins, iNumChannels);
    calcFreqTickLabels(pfTmpTickLabels, iNumFreqBins, fSampleRateInHz);
    Result.setAxisTickLabels(CMatrixRepresentationResult::kCol, pfTmpTickLabels, iNumFreqBins);
    calcTimeTickLabels(pfTmpTickLabels, iNumObservations, iBlockLength, iHopLength, fSampleRateInHz);
    Result.setAxisTickLabels(CMatrixRepresentationResult::kRow, pfTmpTickLabels, iNumObservations);

    m_pCResult      = &Result;
   
    
    m_aiParams[kBlockLength]  = iBlockLength;
    m_aiParams[kHopLength]    = iHopLength;
    m_iNumChannels  = iNumChannels;

    // free local memory
    delete [] pfTmpTickLabels;
    
    m_bIsInitialized    = true;

    return kNoError;
}

int CMatrixRepresentation::calcNumProcessingBlocks( int iLengthInSamples, int iBlockLength, int iHopLength )
{
    return static_cast<int>(ceil((iLengthInSamples + iBlockLength - iHopLength) * 1.F/iHopLength));
}

void CMatrixRepresentation::calcFreqTickLabels( float * pfTickLabels, int iNumFreqBins, float fSampleRateInHz )
{
    assert(pfTickLabels);

    for (int k = 0; k < iNumFreqBins; k++)
        pfTickLabels[k]  = k * fSampleRateInHz * .5F / (iNumFreqBins-1);
}

void CMatrixRepresentation::calcTimeTickLabels( float * pfTickLabels, int iNumObservations, int iBlockLength, int iHopLength, float fSampleRateInHz )
{
    assert(pfTickLabels);
    
    pfTickLabels[0] = fSampleRateInHz * (-iBlockLength*.5F + iHopLength);

    for (int i = 1; i < iNumObservations; i++)
    {
        pfTickLabels[i] = i * iHopLength / fSampleRateInHz + pfTickLabels[0];
    }
}

Error_t CMatrixRepresentation::process( float **ppfInputBuffer, int iNumOfFrames )
{
    int iBlockIdx = 0;

    // add new input data of any length
    m_pCInputBuffer->setDataPtr2Hold(ppfInputBuffer, iNumOfFrames);
    while (m_pCInputBuffer->getBlock(m_ppfProcBuffer, m_aiParams[kBlockLength], m_aiParams[kHopLength]))
    {
        for (int c=0; c < m_iNumChannels; c++)
        {
            Error_t rErr;
            // inplace fft
            m_pCFft->doFft(m_ppfProcBuffer[c], m_ppfProcBuffer[c]);
            m_pCFft->getMagnitude(m_ppfProcBuffer[c], m_ppfProcBuffer[c]);
            rErr = m_pCResult->setResultRow(iBlockIdx, m_ppfProcBuffer[c], m_aiParams[kBlockLength]/2+1, c);

            if (rErr != kNoError)
                return rErr;
        }
        iBlockIdx++;
    }
    // store remaining data
    m_pCInputBuffer->releaseDataPtr();

    for (int c=0; c < m_iNumChannels; c++)
    {
       CUtil::setZero(m_ppfProcBuffer[c], m_aiParams[kBlockLength]);
    }
    // at the end of the processing - get remaining frames from the internal buffer if needed
    m_pCInputBuffer->flush(m_ppfProcBuffer);
    for (int c=0; c < m_iNumChannels; c++)
    {
        Error_t rErr;
        // inplace fft
        m_pCFft->doFft(m_ppfProcBuffer[c], m_ppfProcBuffer[c]);
        m_pCFft->getMagnitude(m_ppfProcBuffer[c], m_ppfProcBuffer[c]);
        m_pCResult->setResultRow(iBlockIdx, m_ppfProcBuffer[c], m_aiParams[kBlockLength]/2+1, c);

        rErr = m_pCResult->setResultRow(iBlockIdx, m_ppfProcBuffer[c], m_aiParams[kBlockLength]/2+1, c);

        if (rErr != kNoError)
            return rErr;
    }

    return kNoError;

}

int CMatrixRepresentation::getParam( Params_t eParam ) const
{
    return m_aiParams[eParam];
}


