#if !defined(__Nmf_hdr__)
#define __Nmf_hdr__

#include "NmfIf.h"


class CNmf : public CNmfIf
{
public:
    Error_t init (CNmfParametrization &NmfSharedData) override;
    Error_t reset () override;

    Error_t process (const CMatrix *pCInput, CNmfResult& NmfResult) override;

    void adaptTemplateXCorr( CMatrix * paaCMatrix[kNumMatrices][kNumSplits], const int aiRank[], float fAdaptCoeff );

    CMatrix calcCrossCorr( CMatrix * paCMatrix[kNumSplits] );


    CNmf ();
    virtual ~CNmf ();

private:
    void runNmf( const CMatrix * pCInput, CNmfResult &NmfResult );
    float calcKlDivergence(const CMatrix &Mat1, const CMatrix &Mat2) const;
    bool isRelativeErrorBelowThresh(int iCurrIteration, const float *pfErr) const;

    static const float m_kMinOffset;

    CNmfParametrization *m_phCNmfConfig;
};

#endif // #if !defined(__Nmf_hdr__)



