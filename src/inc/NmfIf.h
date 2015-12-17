#if !defined(__NmfIf_hdr__)
#define __NmfIf_hdr__

#include "ErrorDef.h"

enum MatrixSplit_t
{
    kSplit1,
    kSplit2,

    kNumSplits
};
enum Matrices_t
{
    kDict,
    kAct,

    kNumMatrices
};

class CMatrix;

class CNmfParametrization
{
    friend class CNmfResult;
public:

    enum AdaptationMethod_t
    {
        kNoAdaptation,
        kAdaptation01,
        kAdaptation02,

        kNumAdaptationMethods
    };

    CNmfParametrization();
    virtual ~CNmfParametrization();

    Error_t init(int iTemplateLength, int iRankSplit1, int iRankSplit2 = 0);
    Error_t reset ();

    Error_t setMatrixInit(Matrices_t eMatrix, MatrixSplit_t eSplit, const CMatrix *pCMatrix);
    Error_t setIsUpdated(Matrices_t eMatrix, MatrixSplit_t eSplit, bool bIsUpdate = true);
    Error_t setNmfTerminationCriteria (int iMaxIterations = 300, float fMinError = 1e-3F);
    Error_t setAdaptTerminationCriteria (int iMaxIterations = 20, float fMinError = 1e-2F);
    Error_t setAdaptRho (float fRho = .5F);
    Error_t setSparsityLambda(MatrixSplit_t eSplit, float fValue = 0);
    Error_t setAdaptationMethod(AdaptationMethod_t eMethod = kNoAdaptation);

    // get settings
    int     getNmfMaxIterations() const;
    int     getAdaptMaxIterations() const;
    float   getNmfMinError() const;
    float   getAdaptMinError() const;
    float   getAdaptRhoThresh() const;
    float   getSparsityLambda(MatrixSplit_t eSplit) const;
    bool    getIsUpdated(Matrices_t eMatrix, MatrixSplit_t eSplit) const;
    int     getRank(MatrixSplit_t eSplit) const;
    int     getTemplateLength() const;
    AdaptationMethod_t getAdaptationMethod() const;

protected:
    const CMatrix* getMatrixPtr(Matrices_t eMatrix, MatrixSplit_t eSplit);

private:
    CNmfParametrization (const CNmfParametrization &/*other*/) {};
    CMatrix *m_aapCMatrices[kNumMatrices][kNumSplits];

    AdaptationMethod_t  m_eMethod;
    int m_iTemplateLength;
    bool m_aabIsUpdated[kNumMatrices][kNumSplits];
    float m_afSparsity[kNumSplits];
    int m_aiRank[kNumSplits];

    int m_iNmfMaxIter;
    float m_fNmfMinError;
    int m_iAdaptMaxIter;
    float m_fAdaptMinError;
    float m_fAdaptRho;
};

class CNmfResult
{
    friend class CNmf;
public:
    CNmfResult();
    virtual ~CNmfResult();


    // get result
    int     getNumNmfIterations() const;
    int getNumAdaptIterations() const;
    CMatrix getMatrix(Matrices_t eMatrix, MatrixSplit_t eSplit) const;
    float   getNmfError() const;
    float getAdaptError() const;
protected:

    bool    isInitialized() const {return m_bIsInitialized;}
    Error_t init(CNmfParametrization& ParamsAndInit, int iNumObservations);
    Error_t reset();
    CMatrix* getMatrixPtr(Matrices_t eMatrix, MatrixSplit_t eSplit);
    float*  getNmfErrorPtr();
    float* getAdaptErrorPtr();

    int m_iNumNmfIter;
    int m_iNumAdaptIter;
private:
    CNmfResult (const CNmfResult &/*other*/) {};
    CMatrix *m_aapCMatrices[kNumMatrices][kNumSplits];

    float *m_pfAdaptError;
    float *m_pfNmfError;
    int m_iMaxNmfIter;
    int m_iMaxAdaptIter;

    bool m_bIsInitialized;
};


class CNmfIf
{
public:
    /*! version number */
    enum Version_t
    {
        kMajor,                         //!< major version number
        kMinor,                         //!< minor version number
        kPatch,                         //!< patch version number

        kNumVersionInts
    };

    static const int  getVersion (const Version_t eVersionIdx);
    static const char* getBuildDate ();

    static Error_t create (CNmfIf*& pCInstance);
    static Error_t destroy (CNmfIf*& pCInstance);
    
    virtual Error_t init (CNmfParametrization &NmfInit) = 0;
    virtual Error_t reset () = 0;
    
    virtual Error_t process (const CMatrix *pCInput, CNmfResult& NmfResult) = 0;

protected:
    CNmfIf () {};
    virtual ~CNmfIf () {};
};

#endif // #if !defined(__NmfIf_hdr__)



