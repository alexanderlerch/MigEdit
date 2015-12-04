#if !defined(__Nmf_hdr__)
#define __Nmf_hdr__

#include "NmfIf.h"


class CNmf : public CNmfIf
{
public:
    Error_t init (CNmfSharedData &NmfSharedData) override;
    Error_t reset () override;

    Error_t process (const CMatrix *pCInput) override;

    CNmf ();
    virtual ~CNmf ();

private:
    float calcKlDivergence(const CMatrix &Mat1, const CMatrix &Mat2) const;
    bool isRelativeErrorBelowThresh(int iCurrIteration) const;

    static const float m_kMinOffset;

    float *m_pfErr;
    CNmfSharedData *m_pCConfigAndResults;
};

#endif // #if !defined(__Nmf_hdr__)



