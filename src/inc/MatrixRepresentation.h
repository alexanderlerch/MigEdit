#if !defined(__MatrixRepresentation_hdr__)
#define __MatrixRepresentation_hdr__

#include "ErrorDef.h"
#include "Matrix.h"
#include "InputBuffSrc.h"

class CFft;

class CMatrixRepresentationResult
{
    friend class CMatrixRepresentation;
public:
    enum Axes_t
    {
        kRow,
        kCol,

        kNumAxes
    };
    CMatrixRepresentationResult();;
    virtual ~CMatrixRepresentationResult();
    Error_t reset ();

    int         getNumRows() const;
    int         getNumCols() const;
    int         getNumChannels() const;
    Error_t     getAxisTickLabels (enum Axes_t eAxis, float *pfAxisTickLabels, int iNumOfAxisTickLabels) const;
    std::string getAxisLabel (enum Axes_t eAxis) const;
    std::string getAxisUnit (enum Axes_t eAxis) const;

    Error_t     getCol (int iCol, float *pfValues, int iNumValues, int iChannelIdx = 0) const;
    Error_t     getRow (int iRow, float *pfValues, int iNumValues, int iChannelIdx = 0) const;
    CMatrix*    getResultPtr (int iChannelIdx = 0) ;

protected:
    Error_t init(int iNumRows, int iNumCols, int iNumChannels = 1);
    Error_t setResultCol (int iCol, const float *pfValues, int iNumValues, int iChannelIdx = 0);
    Error_t setResultRow (int iRow, const float *pfValues, int iNumValues, int iChannelIdx = 0);
    Error_t setResult (const CMatrix & pCValues, int iChannelIdx = 0);
    Error_t setAxisTickLabels (enum Axes_t eAxis, const float *pfAxisTickLabels, int iNumOfAxisTickLabels);
    Error_t setAxisLabel (enum Axes_t eAxis, std::string sAxisLabel);
    Error_t setAxisUnit (enum Axes_t eAxis, std::string sAxisLabel);
private:
    std::string m_asLabels[kNumAxes];
    std::string m_asUnits[kNumAxes];

    CMatrix  **m_ppCResult;

    bool    m_bIsInitialized;

    float   *m_apfTickLabels[kNumAxes];
    int     m_iNumChannels;
};

class CMatrixRepresentation
{
public:
    /*! version number */
    enum MatrixRepresentations_t
    {
        kSpectrogram,                         
        //kAcf,                       
        //kCepstrum,                       

        kNumMatrixRepresentations
    };

    enum Params_t
    {
        kBlockLength,
        kHopLength,

        kNumParams
    };

    /*! creates a new CMatrixRepresentation instance
    \param CMatrixRepresentation * & pCInstance: pointer to the new instance
    \return Error_t
    */
    static Error_t create (CMatrixRepresentation*& pCInstance);
    /*! destroys a CMatrixRepresentation instance
    \param CMatrixRepresentation * & pCInstance: pointer to instance to be destroyed
    \return Error_t
    */
    static Error_t destroy (CMatrixRepresentation*& pCInstance);
    
    Error_t init (int iLengthInSamples, float fSampleRateInHz, int iNumChannels, int iBlockLength, int iHopLength, CMatrixRepresentationResult &Result);

    /*! resets the instance
    \return Error_t
    */
    Error_t reset ();

    int     getParam (Params_t eParam) const;

    Error_t process (float **ppfInputBuffer, int iNumOfFrames);

protected:
    CMatrixRepresentation ();
    virtual ~CMatrixRepresentation ();
private:

    int calcNumProcessingBlocks( int iLengthInSamples, int iBlockLength, int iHopLength );
    void calcFreqTickLabels( float * pfTickLabels, int iNumFreqBins, float fSampleRateInHz );
    void calcTimeTickLabels( float * pfTickLabels, int iNumObservations, int iBlockLength, int iHopLength, float fSampleRateInHz );


    bool    m_bIsInitialized;

    int     m_aiParams[kNumParams];
    int     m_iNumChannels;

    CFft            *m_pCFft;
    CInputBuffSrc<float>   *m_pCInputBuffer;

    float   **m_ppfProcBuffer;
    CMatrixRepresentationResult *m_pCResult;

};

#endif // #if !defined(__MatrixRepresentation_hdr__)



