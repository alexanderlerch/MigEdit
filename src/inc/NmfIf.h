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

class CNmfSharedData
{
    friend class CNmf;
public:
    CNmfSharedData(int iTemplateLength, int iRankSplit1, int iRankSplit2 = 0);
    virtual ~CNmfSharedData();

    Error_t setMatrix(Matrices_t eMatrix, MatrixSplit_t eSplit, const CMatrix *pCMatrix);
    Error_t setIsUpdated(Matrices_t eMatrix, MatrixSplit_t eSplit, bool bIsUpdate = true);
    Error_t setTerminationCriteria (int iMaxIterations = 300, float fMinError = 1e-4F);
    Error_t setSparsityLambda(MatrixSplit_t eSplit, float fValue = 0);
    Error_t finalizeSettings();
    Error_t reset ();

    // get settings
    int     getMaxIterations() const;
    float   getMinError() const;
    float   getSparsityLambda(MatrixSplit_t eSplit) const;
    bool    getIsUpdated(Matrices_t eMatrix, MatrixSplit_t eSplit) const;
    int     getRank(MatrixSplit_t eSplit) const;
    int     getTemplateLength() const;

    bool    isInitialized() const {return m_bIsInitialized;}
 
    // get result
    int     getResultNumIterations() const;
protected:
    CMatrix* getMatrixPtr(Matrices_t eMatrix, MatrixSplit_t eSplit);
    float*  getErrorPtr();

    int m_iNumIterations;
private:
    CMatrix *m_aapCMatrices[kNumMatrices][kNumSplits];

    float *m_pfError;

    bool m_bIsInitialized;

    int m_iTemplateLength;
    bool m_aabIsUpdated[kNumMatrices][kNumSplits];
    float m_fMinError;
    float m_afSparsity[kNumSplits];
    int m_iMaxIter;
    int m_aiRank[kNumSplits];
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
    
    virtual Error_t init (CNmfSharedData &NmfSharedData) = 0;
    virtual Error_t reset () = 0;
    
    virtual Error_t process (const CMatrix *pCInput) = 0;

protected:
    CNmfIf () {};
    virtual ~CNmfIf () {};
};

#endif // #if !defined(__NmfIf_hdr__)



